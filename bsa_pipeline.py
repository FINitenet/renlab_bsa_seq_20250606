#!/usr/bin/env python3
import os
import json
import subprocess
import pandas as pd
from datetime import datetime

# =============== 配置参数 ===============
threads = 24
genome_fasta = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bwa_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
output_base = "/bios-store1/chenyc/Project/YHR/Project_yhr_20250418BSA-seq"
fastq_dir = os.path.join(output_base, "0_fastq")
meta_tsv = os.path.join(output_base, "fastq_meta.tsv")
meta_json = os.path.join(output_base, "fastq_meta.json")

# =============== 工具函数 ===============
def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}")

def run(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

# =============== Step 0: 生成 fastq_meta 文件 ===============
def generate_fastq_meta(fastq_dir, output_meta_tsv, output_meta_json):
    samples = {}
    for fname in os.listdir(fastq_dir):
        if not fname.endswith(".fastq.gz"):
            continue
        parts = fname.split('.cleaned_')
        if len(parts) != 2:
            continue
        sample, read = parts
        read = read[0]
        if sample not in samples:
            samples[sample] = {"read1": "", "read2": ""}
        full_path = os.path.abspath(os.path.join(fastq_dir, fname))
        if read == '1':
            samples[sample]["read1"] = full_path
        elif read == '2':
            samples[sample]["read2"] = full_path

    df = pd.DataFrame([
        {"sample": sample, "read1": paths["read1"], "read2": paths["read2"]}
        for sample, paths in samples.items()
    ])
    df.to_csv(output_meta_tsv, sep='\t', index=False)
    with open(output_meta_json, 'w') as f:
        json.dump(samples, f, indent=4)
    print(f"[✓] Meta files saved:\n  TSV: {output_meta_tsv}\n  JSON: {output_meta_json}")

# =============== Step 1: FastQC + MultiQC ===============
def fastqc_step(df):
    log("Step 1: FastQC")
    ensure_dir("1_fastqc")
    files = df["read1"].tolist() + df["read2"].tolist()
    run(f"source activate seq && fastqc -t {threads} -o 1_fastqc {' '.join(files)}")
    run("source activate seq && multiqc -f --fullnames -o 1_fastqc 1_fastqc")

# =============== Step 2: BWA 比对 ===============
def mapping_step(df):
    log("Step 2: Mapping with BWA")
    ensure_dir("2_mapping")
    for _, row in df.iterrows():
        prefix = row["sample"]
        read1, read2 = row["read1"], row["read2"]
        raw_bam = f"{prefix}.bam"
        sorted_bam = f"2_mapping/{prefix}_sorted.bam"
        run(f"bwa mem -t {threads} -R '@RG\\tID:{prefix}\\tSM:{prefix}\\tPL:illumina' {genome_fasta} {read1} {read2} | "
            f"samtools view -@ {threads} -bhS -o {raw_bam} -")
        run(f"samtools sort -@ {threads} -o {sorted_bam} {raw_bam}")
        run(f"samtools index -@ {threads} {sorted_bam}")
        os.remove(raw_bam)

# =============== Step 3: 去重 ===============
def dedup_step(df):
    log("Step 3: Deduplication with GATK4")
    ensure_dir("3_dedup")
    for _, row in df.iterrows():
        prefix = row["sample"]
        sorted_bam = f"2_mapping/{prefix}_sorted.bam"
        dedup_bam = f"3_dedup/{prefix}_dedup.bam"
        metrics = f"3_dedup/{prefix}_metrics.txt"
        run(f"gatk --java-options '-Xmx10g' MarkDuplicates "
            f"-I {sorted_bam} -O {dedup_bam} -M {metrics} --REMOVE_DUPLICATES true")
        run(f"samtools index -@ {threads} {dedup_bam}")
        os.remove(metrics)

# =============== Step 4: 变异检测 ===============
def call_variants_bcftools(parent_bam, mutant_bam, output_prefix, output_dir, reference_fasta):
    log("Step 4: Variant calling using bcftools")
    ensure_dir(output_dir)

    raw_vcf = os.path.join(output_dir, f"{output_prefix}_variants.vcf")
    filtered_vcf = os.path.join(output_dir, f"{output_prefix}_variants_filtered.vcf.gz")

    mpileup_cmd = (
        f"bcftools mpileup -C 50 --ignore-RG -a DP,AD,INFO/AD -O u --threads {threads} "
        f"-f {reference_fasta} {parent_bam} {mutant_bam} | "
        f"bcftools call -O v --threads {threads} -f GQ -mv -o {raw_vcf}"
    )
    filter_cmd = (
        f"bcftools filter --threads {threads} -g 5 -G 10 "
        f"-e '%QUAL < 20 || INFO/MQ < 20 || %MIN(FORMAT/DP) < 4' "
        f"-s LowQual -O z -o {filtered_vcf} {raw_vcf}"
    )
    run(mpileup_cmd)
    run(filter_cmd)
    log(f"[✓] Output VCF files:\n  - {raw_vcf}\n  - {filtered_vcf}")

# =============== 主程序入口 ===============
if __name__ == "__main__":
    # Step 0
    generate_fastq_meta(fastq_dir, meta_tsv, meta_json)
    df_meta = pd.read_csv(meta_tsv, sep='\t')

    # Step 1
    fastqc_step(df_meta)

    # Step 2
    mapping_step(df_meta)

    # Step 3
    dedup_step(df_meta)

    # Step 4（你可以修改样本名）
    parent_sample = "YE_L1_UDI164"
    mutant_sample = "Ya10_L7_UDI089"
    parent_bam = f"3_dedup/{parent_sample}_dedup.bam"
    mutant_bam = f"3_dedup/{mutant_sample}_dedup.bam"

    call_variants_bcftools(
        parent_bam=parent_bam,
        mutant_bam=mutant_bam,
        output_prefix="YA10",
        output_dir=os.path.join(output_base, "4_variant"),
        reference_fasta=genome_fasta
    )

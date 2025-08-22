#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :bsa_pipeline.py
# @Author    :Yuchen@rlab
# @Description: ESM BSA-Seq ΔIndex variant analysis pipeline

import os
import json
import argparse
import subprocess
import logging
from datetime import datetime
from pathlib import Path
import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def log(msg):
    logging.info(msg)

def run(cmd):
    log(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def match_sample_names(reader_samples, parent_key, mutant_key):
    def match(key):
        matches = [s for s in reader_samples if key in s]
        if not matches:
            raise ValueError(f"No sample matched for {key}")
        elif len(matches) > 1:
            raise ValueError(f"Multiple samples matched for {key}: {matches}")
        return matches[0]

    return match(parent_key), match(mutant_key)

def plot_facet_by_chr(df, parent, mutant, output_dir):
    sns.set(style="whitegrid")

    # 设置绘图尺寸和分面
    g = sns.FacetGrid(df, col="Chr", col_wrap=5, sharey=True, height=3, aspect=1.5)
    g.map_dataframe(sns.scatterplot, x="Pos", y="Delta_Index", s=5, linewidth=0)

    # 设置每个子图的标签
    g.set_axis_labels("Position", f"{mutant} - {parent}")
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Delta Index per Chromosome", fontsize=9)

    # 输出保存
    pdf_out = os.path.join(output_dir, cfg["mutant_sample"] + "_Mapping.pdf")
    g.savefig(pdf_out)
    log(f"[✓] Faceted plot saved to {pdf_out}")

def generate_fastq_meta(fastq_dir, output_meta_tsv, output_meta_json):
    log("Step 1: Generating FASTQ meta")
    samples = {}
    for fname in os.listdir(fastq_dir):
        if not fname.endswith(".fastq.gz"):
            continue
        parts = fname.split('.R')
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
    return df

def fastqc_step(df, threads, output_dir):
    log("Step 2: FastQC & MultiQC")
    fastqc_dir = os.path.join(output_dir, "1_fastqc")
    ensure_dir(fastqc_dir)
    files = df["read1"].tolist() + df["read2"].tolist()
    run(f"fastqc -t {threads} -o {fastqc_dir} {' '.join(files)}")
    run(f"multiqc -f --fullnames -o {fastqc_dir} {fastqc_dir}")

def mapping_step(df, genome_fasta, threads, output_dir):
    log("Step 3: Mapping with BWA")
    map_dir = os.path.join(output_dir, "2_mapping")
    ensure_dir(map_dir)
    for _, row in df.iterrows():
        prefix = row["sample"]
        read1, read2 = row["read1"], row["read2"]
        raw_bam = f"{prefix}.bam"
        sorted_bam = f"{map_dir}/{prefix}_sorted.bam"
        run(f"bwa mem -t {threads} -R '@RG\\tID:{prefix}\\tSM:{prefix}\\tPL:illumina' {genome_fasta} {read1} {read2} | "
            f"samtools view -@ {threads} -bhS -o {raw_bam} -")
        run(f"samtools sort -@ {threads} -o {sorted_bam} {raw_bam}")
        run(f"samtools index -@ {threads} {sorted_bam}")
        os.remove(raw_bam)

def dedup_step(df, threads, output_dir):
    log("Step 4: MarkDuplicates")
    dedup_dir = os.path.join(output_dir, "3_dedup")
    ensure_dir(dedup_dir)
    for _, row in df.iterrows():
        prefix = row["sample"]
        sorted_bam = os.path.join(output_dir, "2_mapping", f"{prefix}_sorted.bam")
        dedup_bam = os.path.join(dedup_dir, f"{prefix}_dedup.bam")
        metrics = os.path.join(dedup_dir, f"{prefix}_metrics.txt")
        run(f"gatk --java-options '-Xmx10g' MarkDuplicates -I {sorted_bam} -O {dedup_bam} -M {metrics} --REMOVE_DUPLICATES true")
        run(f"samtools index -@ {threads} {dedup_bam}")
        os.remove(metrics)

def call_variants(parent_bam, mutant_bam, output_prefix, output_dir, reference_fasta):
    log("Step 5: Variant calling")
    ensure_dir(output_dir)
    raw_vcf = os.path.join(output_dir, f"{output_prefix}_variants.vcf")
    filtered_vcf = os.path.join(output_dir, f"{output_prefix}_variants_filtered.vcf.gz")
    run(
        f"bcftools mpileup --threads 24 -C 50 --ignore-RG -a DP,AD,INFO/AD -O u --threads 16 "
        f"-f {reference_fasta} {parent_bam} {mutant_bam} | "
        f"bcftools call -O v --threads 24 -f GQ -mv -o {raw_vcf}"
    )
    run(
        f"bcftools filter --threads 24 -g 5 -G 10 "
        f"-e 'QUAL < 20 || MQ < 20 || INFO/DP < 4' "
        f"-s LowQual -O z -o {filtered_vcf} {raw_vcf}"
    )
    return filtered_vcf

def parse_and_plot_vcf(vcf_path, parent, mutant, output_dir):
    log("Step 6: Plotting ΔIndex Manhattan plot")
    reader = vcf.Reader(filename=vcf_path)

    parent_real, mutant_real = match_sample_names(reader.samples, parent, mutant)
    log(f"[✓] Matched parent: {parent} → {parent_real}")
    log(f"[✓] Matched mutant: {mutant} → {mutant_real}")

    records = []
    for record in reader:
        if record.FILTER or record.is_indel: continue
        if len(record.ALT) != 1: continue
        ref, alt = str(record.REF), str(record.ALT[0])
        # if not ((ref == "G" and alt == "A") or (ref == "C" and alt == "T")): continue

        ad_parent = record.genotype(parent_real).data.AD
        ad_mutant = record.genotype(mutant_real).data.AD
        dp_parent = sum(ad_parent) if ad_parent else 0
        dp_mutant = sum(ad_mutant) if ad_mutant else 0
        idx_parent = ad_parent[1]/dp_parent if dp_parent else 0
        idx_mutant = ad_mutant[1]/dp_mutant if dp_mutant else 0
        delta = idx_mutant - idx_parent
        records.append({
            "Chr": record.CHROM, "Pos": record.POS,
            "Ref": ref, "Alt": alt,
            "Delta_Index": delta
        })

    df = pd.DataFrame(records)
    if df.empty:
        log("[!] Warning: No valid SNPs found. Check filtering conditions.")
        return
    df["Chr"] = pd.Categorical(df["Chr"], ordered=True)
    excel_out = os.path.join(output_dir, cfg["mutant_sample"] + "_Mapping.xlsx")
    df.to_excel(excel_out, index=False)

    log(f"[✓] Excel saved to {excel_out}")

    # Plotting
    plot_facet_by_chr(df, parent, mutant, output_dir)

# === Main ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BSA-Seq ΔIndex pipeline")
    parser.add_argument("-c", "--config", type=str, default="config.json", help="Path to config JSON")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    ensure_dir(cfg["output_dir"])
    setup_logging(cfg["output_dir"])
    log("=== BSA-Seq Pipeline Start ===")

    df_meta = generate_fastq_meta(cfg["fastq_dir"],
                                   os.path.join(cfg["output_dir"], "fastq_meta.tsv"),
                                   os.path.join(cfg["output_dir"], "fastq_meta.json"))
    fastqc_step(df_meta, cfg["threads"], cfg["output_dir"])
    mapping_step(df_meta, cfg["genome_fasta"], cfg["threads"], cfg["output_dir"])
    dedup_step(df_meta, cfg["threads"], cfg["output_dir"])

    parent_bam = os.path.join(cfg["output_dir"], "3_dedup", f"{cfg['parent_sample']}_dedup.bam")
    mutant_bam = os.path.join(cfg["output_dir"], "3_dedup", f"{cfg['mutant_sample']}_dedup.bam")
    vcf_out = call_variants(parent_bam, mutant_bam, cfg["mutant_sample"],
                            os.path.join(cfg["output_dir"], "4_variant"), cfg["genome_fasta"])
    parse_and_plot_vcf(vcf_out, cfg["parent_sample"], cfg["mutant_sample"], cfg["output_dir"])
    log("[✓] Pipeline finished.")

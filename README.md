## Purpose and Scope
This document provides a comprehensive overview of the BSA-seq (Bulk Segregant Analysis sequencing) analysis pipeline repository. The system implements a complete genomic variant analysis workflow for identifying genetic variants associated with specific traits through comparative analysis of parent and mutant sample populations.
The pipeline processes raw FASTQ sequencing data through quality control, read mapping, duplicate removal, variant calling, and specialized BSA analysis tools. For detailed information about individual analysis components, see Specialized Analysis Tools. For system configuration and external dependencies, see System Configuration and External Dependencies.

## Data Processing Workflow
The pipeline implements a linear data processing workflow with the following stages:

| Stage | Function            | Input                | Output                              | External Tools                     |
| ----: | ------------------- | -------------------- | ----------------------------------- | ---------------------------------- |
|     1 | Metadata Generation | FASTQ directory      | `fastq_meta.tsv`, `fastq_meta.json` | –                                  |
|     2 | Quality Control     | FASTQ files          | FastQC/MultiQC reports              | `fastqc`, `multiqc`                |
|     3 | Read Mapping        | FASTQ + Reference    | Sorted BAM files                    | `bwa mem`, `samtools`              |
|     4 | Duplicate Removal   | Sorted BAMs          | Deduplicated BAMs                   | `gatk MarkDuplicates`              |
|     5 | Variant Calling     | Parent + Mutant BAMs | Filtered VCF                        | `bcftools`                         |
|     6 | Visualization       | VCF + Sample names   | Delta Index plots                   | `matplotlib`, `seaborn` *(Python)* |


## Specialized Analysis Tools
The pipeline integrates three specialized analysis tools for post-variant calling analysis:

### 1) Delta Index Analysis — bsa_delta_index.R
Purpose: Calculate and visualize Delta Index (Sample − Parent) for BSA mapping, focusing on G→A / C→T transitions.
Input: Filtered VCF, parent/sample name prefixes (after BAM name cleanup).
Output: PDF plot and XLSX table, named as <sample>vs<parent>_Mapping.*.

Usage:

Rscript scripts/bsa_delta_index.R \
  -p <PARENT_PREFIX> \
  -s <SAMPLE_PREFIX> \
  -v <FILTERED_VCF.gz> \
  -o <OUTDIR>
Notes:
列名前缀按规则自动清洗（移除路径、_dedup*、_sorted*、.bam 等后缀）。
亲本粗过滤：GT == "0/0" 的位点保留。
自动防零除；输出包含 Parent_Index、Sample_Index、Delta_Index。

### 2) QTL-seq Analysis — method_qtlseq.R
Purpose: Statistical mapping of QTL using QTLseqr.
Input: GATK-formatted variant tables (or converted from VCF).
Output: CSV results and QTL mapping plots.

Usage:

Rscript scripts/method_qtlseq.R \
  --input table_for_qtlseq.csv \
  --outdir 5_bsa/qtlseq \
  --win-size 1e6 \
  --step-size 1e5

### 3) VCF Comparison — compare_vcf_gatk_bcftools.py
Purpose: Cross-validate variants between GATK 与 BCFtools 调用结果。
Input: 成对的 VCF 文件。
Output: Excel 对比报告（匹配/不匹配位点、统计汇总）与 Δ 值直方图。

Usage:

python scripts/compare_vcf_gatk_bcftools.py \
  --gatk 4_variant/gatk_filtered.vcf.gz \
  --bcf  4_variant/bcftools_filtered.vcf.gz \
  --out  5_bsa/vcf_compare

### Configuration

参考序列、索引路径等在 config.json里维护：
reference: "/path/to/reference.fa"
bwa_index: "/path/to/reference.fa"
dbsnp: "/path/to/dbsnp.vcf.gz"

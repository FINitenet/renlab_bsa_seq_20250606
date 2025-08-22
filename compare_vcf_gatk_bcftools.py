#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :compare_vcf_gatk_bcftools.py
# @Time      :2025/08/22 17:39:23
# @Author    :Yuchen@rlab
# @Description: Compare GATK and BCFtools VCF files for ΔIndex variant analysis

import argparse
import pandas as pd
import vcf
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_vcf(path, tag):
    reader = vcf.Reader(filename=path)
    records = []
    for record in reader:
        if record.is_snp and len(record.ALT) == 1:
            for sample in record.samples:
                gt = sample['GT'] if sample['GT'] else './.'
                records.append({
                    "Chr": record.CHROM,
                    "Pos": record.POS,
                    "Ref": str(record.REF),
                    "Alt": str(record.ALT[0]),
                    "Sample": sample.sample,
                    f"{tag}_GT": gt,
                    f"{tag}_DP": sample.get('DP'),
                    f"{tag}_GQ": sample.get('GQ'),
                    f"{tag}_AD": sample.get('AD')
                })
    return pd.DataFrame(records)

def compare_variants(df1, df2):
    # Merge on Chr:Pos:Ref:Alt + Sample
    merged = pd.merge(df1, df2, on=["Chr", "Pos", "Ref", "Alt", "Sample"], how="outer", indicator=True)
    return merged

def plot_delta_index(df, out_pdf):
    df = df.dropna(subset=["bcf_AD", "gatk_AD"])
    def get_index(ad):
        return ad[1] / sum(ad) if ad and isinstance(ad, (list, tuple)) and sum(ad) > 0 else None
    df["bcf_idx"] = df["bcf_AD"].apply(get_index)
    df["gatk_idx"] = df["gatk_AD"].apply(get_index)
    df["Delta"] = df["gatk_idx"] - df["bcf_idx"]

    plt.figure(figsize=(10,4))
    sns.histplot(df["Delta"].dropna(), bins=100, kde=True, color="skyblue")
    plt.title("Delta Index (GATK - BCFtools)")
    plt.xlabel("GATK - BCFtools Index")
    plt.tight_layout()
    plt.savefig(out_pdf)
    print(f"[✓] Delta_Index histogram saved: {out_pdf}")

def main(args):
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    print(f"🔍 Loading GATK VCF: {args.gatk_vcf}")
    df_gatk = load_vcf(args.gatk_vcf, "gatk")
    print(f"🔍 Loading bcftools VCF: {args.bcf_vcf}")
    df_bcf = load_vcf(args.bcf_vcf, "bcf")

    print("🔍 Comparing variants...")
    df_merged = compare_variants(df_gatk, df_bcf)

    # Save to Excel
    out_xlsx = args.out_prefix + ".xlsx"
    with pd.ExcelWriter(out_xlsx) as writer:
        df_gatk.to_excel(writer, sheet_name="GATK_VCF", index=False)
        df_bcf.to_excel(writer, sheet_name="BCFtools_VCF", index=False)
        df_merged.to_excel(writer, sheet_name="Merged_Comparison", index=False)
    print(f"[✓] Excel saved: {out_xlsx}")

    # Plot delta index
    out_pdf = args.out_prefix + "_delta.pdf"
    plot_delta_index(df_merged, out_pdf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gatk_vcf", required=True, help="GATK VCF path (.vcf.gz)")
    parser.add_argument("--bcf_vcf", required=True, help="bcftools VCF path (.vcf.gz)")
    parser.add_argument("--out_prefix", default="compare_result", help="Output prefix (no extension)")
    args = parser.parse_args()
    main(args)

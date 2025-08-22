#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(VariantAnnotation)
  library(tidyverse)
  library(ragg)
  library(writexl)
  library(GenomicRanges)
})

# -------------------- CLI --------------------
option_list <- list(
  make_option(c("-p", "--parent"), type = "character",
              help = "亲本样本名前缀（清洗后），如：YE_L1_UDI164"),
  make_option(c("-s", "--sample"), type = "character",
              help = "目标样本名前缀（清洗后），如：Ya10_L7_UDI089"),
  make_option(c("-v", "--vcf"), type = "character",
              help = "输入 VCF(.vcf.gz) 文件路径"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "输出目录（不存在会自动创建）")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("parent", "sample", "vcf", "outdir")
missing_args <- required[!nzchar(trimws(unlist(opt[required])))]
if (length(missing_args) > 0) {
  cat("缺少参数：", paste(missing_args, collapse = ", "), "\n", sep = "")
  cat("示例：\n",
      "Rscript bsa_delta_index.R ",
      "-p YE_L1_UDI164 -s Ya10_L7_UDI089 ",
      "-v 4_variant/YA10_variants_filtered.vcf.gz -o results\n", sep = "")
  quit(status = 1)
}

if (!file.exists(opt$vcf)) {
  stop("找不到 VCF 文件：", opt$vcf)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

parent_key <- opt$parent
sample_key <- opt$sample
vcf_path   <- opt$vcf
outdir     <- opt$outdir
tag        <- paste0(sample_key, "vs", parent_key)

pdf_file  <- file.path(outdir, paste0(tag, "_Mapping.pdf"))
xlsx_file <- file.path(outdir, paste0(tag, "_Mapping.xlsx"))

# -------------------- Helpers --------------------
clean_sample <- function(x) {
  x %>%
    gsub(".*/", "", .) %>%                    # 去路径
    gsub("_dedup.*|_sorted.*|\\.bam", "", .)  # 去冗余后缀
}

mb_fmt <- scales::unit_format(accuracy = 1, scale = 1e-6, unit = "")

find_col <- function(df, key, suffix) {
  hits <- grep(paste0("^", key, ".*_", suffix, "$"), names(df), value = TRUE)
  if (length(hits) != 1) {
    stop("无法唯一定位列：前缀=", key, " 后缀=", suffix,
         "；匹配到(", length(hits), "): ", paste(hits, collapse = ", "))
  }
  hits
}

extract_AD_mats <- function(ADx) {
  # 情况 A：3D array: [variants x samples x alleles]
  if (is.array(ADx) && length(dim(ADx)) == 3 && dim(ADx)[3] >= 2) {
    AD_REF <- ADx[, , 1, drop = FALSE]
    AD_ALT <- ADx[, , 2, drop = FALSE]
    `dim<-`(AD_REF, dim(AD_REF)[1:2])
    `dim<-`(AD_ALT, dim(AD_ALT)[1:2])
    rownames(AD_REF) <- rownames(ADx); colnames(AD_REF) <- colnames(ADx)
    rownames(AD_ALT) <- rownames(ADx); colnames(AD_ALT) <- colnames(ADx)
    return(list(REF = AD_REF, ALT = AD_ALT))
  }
  # 情况 B：矩阵/数据框，每个单元是长度为2的向量
  if (is.matrix(ADx) || inherits(ADx, "DataFrame") || is.data.frame(ADx)) {
    nvar <- nrow(ADx); nsmp <- ncol(ADx)
    AD_REF <- matrix(NA_integer_, nrow = nvar, ncol = nsmp,
                     dimnames = list(rownames(ADx), colnames(ADx)))
    AD_ALT <- AD_REF
    for (j in seq_len(nsmp)) {
      col_list <- ADx[, j]
      AD_REF[, j] <- vapply(col_list, function(x) {
        x <- as.integer(x); if (length(x) >= 1) x[1] else NA_integer_
      }, 1L)
      AD_ALT[, j] <- vapply(col_list, function(x) {
        x <- as.integer(x); if (length(x) >= 2) x[2] else NA_integer_
      }, 1L)
    }
    return(list(REF = AD_REF, ALT = AD_ALT))
  }
  stop("无法识别的 AD 数据结构：", paste(class(ADx), collapse = ", "))
}

# -------------------- 主流程 --------------------
message("[1/6] 读取 VCF：", vcf_path)
vcf <- readVcf(
  vcf_path,
  param = ScanVcfParam(info = c("INDEL", "MQ"), geno = c("GT", "PL", "DP", "AD", "GQ"))
)

message("[2/6] 过滤 PASS + 二等位 SNP（排除 INDEL）")
keep <- rowRanges(vcf)$FILTER == "PASS" &
  VariantAnnotation::isSNV(vcf) &
  elementNROWS(alt(vcf)) == 1 &
  !isTRUE(as.logical(info(vcf)$INDEL))
vcf <- vcf[keep]
if (length(vcf) == 0) stop("过滤后无变异！")

# 位点信息：Chr/Pos 来自 rowRanges；等位基因来自 ref/alt(vcf)
rr  <- rowRanges(vcf)
REF <- as.character(VariantAnnotation::ref(vcf))
ALT <- as.character(unlist(VariantAnnotation::alt(vcf)))

site_df <- tibble(
  Chr = as.character(GenomicRanges::seqnames(rr)),
  Pos = as.integer(GenomicRanges::start(rr)),
  Ref = REF,
  Alt = ALT
)

# 基因型（用于亲本粗过滤）
GT <- geno(vcf)$GT %>% as.data.frame(check.names = FALSE)
colnames(GT) <- paste0(clean_sample(colnames(GT)), "_GT")

parent_gt_col <- grep(paste0("^", parent_key, ".*_GT$"), colnames(GT), value = TRUE)
if (length(parent_gt_col) != 1) {
  stop("无法唯一定位亲本 GT 列（前缀=", parent_key, "），匹配：",
       paste(parent_gt_col, collapse = ", "))
}

message("[3/6] 亲本粗过滤：", parent_gt_col, " == '0/0'")
mask_parent <- GT[[parent_gt_col]] == "0/0"
vcf     <- vcf[mask_parent]
site_df <- site_df[mask_parent, , drop = FALSE]
GT      <- GT[mask_parent, , drop = FALSE]
if (length(vcf) == 0) stop("亲本粗过滤后无变异！")

# AD/DP/GQ
message("[4/6] 提取 AD/DP/GQ 并整理列名")
AD_raw <- geno(vcf)$AD
if (is.null(AD_raw)) stop("VCF 中不含 AD 字段")

samples_raw <- if (!is.null(colnames(AD_raw))) colnames(AD_raw) else colnames(geno(vcf)$GT)
samples_cln <- clean_sample(samples_raw)

adm <- extract_AD_mats(AD_raw)
AD_REF <- as.data.frame(adm$REF, check.names = FALSE)
AD_ALT <- as.data.frame(adm$ALT, check.names = FALSE)
colnames(AD_REF) <- paste0(samples_cln, "_REF")
colnames(AD_ALT) <- paste0(samples_cln, "_ALT")

DP <- geno(vcf)$DP %>% as.data.frame(check.names = FALSE)
colnames(DP) <- paste0(clean_sample(colnames(DP)), "_DP")

GQ <- geno(vcf)$GQ %>% as.data.frame(check.names = FALSE)
colnames(GQ) <- paste0(clean_sample(colnames(GQ)), "_GQ")

QUAL <- qual(vcf)
MQ   <- info(vcf)$MQ

df <- bind_cols(
  site_df,
  tibble(QUAL = QUAL, MQ = as.numeric(MQ)),
  GT, GQ, DP, AD_REF, AD_ALT
)

# 仅保留转换：G>A 或 C>T
df <- df %>% filter((Ref == "G" & Alt == "A") | (Ref == "C" & Alt == "T"))
if (nrow(df) == 0) stop("没有满足 G>A / C>T 的位点！")

# 计算 Index（自动定位列）
parent_ref_col <- find_col(df, parent_key, "REF")
parent_alt_col <- find_col(df, parent_key, "ALT")
sample_ref_col <- find_col(df, sample_key, "REF")
sample_alt_col <- find_col(df, sample_key, "ALT")

message("[5/6] 计算 Index 与 Delta_Index")
df <- df %>%
  mutate(
    Parent_Index = .data[[parent_alt_col]] /
      pmax(.data[[parent_ref_col]] + .data[[parent_alt_col]], 1e-9),
    Sample_Index = .data[[sample_alt_col]] /
      pmax(.data[[sample_ref_col]] + .data[[sample_alt_col]], 1e-9),
    Delta_Index  = Sample_Index - Parent_Index
  )

# 绘图（仅 chr1..n）
message("[6/6] 绘图与导出 -> ", pdf_file, " / ", xlsx_file)
p <- df %>%
  filter(grepl("^chr\\d+$", Chr)) %>%
  ggplot(aes(x = Pos, y = Delta_Index, color = Chr)) +
  geom_point(size = 0.9, alpha = 0.8) +
  scale_x_continuous(labels = mb_fmt) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1)) +
  labs(x = NULL, y = paste0(sample_key, " - ", parent_key)) +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
  theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.spacing = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 9, family = "Arial"),
    axis.title.y = element_text(margin = margin(r = 7)),
    legend.position = "none"
  )

ggsave(pdf_file, p, width = 23, height = 10, units = "cm", device = cairo_pdf)
write_xlsx(df, xlsx_file)

message("✅ 完成：", tag)

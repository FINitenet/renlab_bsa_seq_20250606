library(VariantAnnotation)
library(tidyverse)
library(ragg)

variant <- readVcf("4_variant/YA10_variants_filtered.vcf.gz",
                   param = ScanVcfParam(info = c("INDEL", "MQ"),
                                        geno = c("GT", "PL", "DP", "AD", "GQ")))

# 1. 保留质量较好的，只有一种变异的SNP
variant <- variant[rowRanges(variant)$FILTER == "PASS" & elementNROWS(alt(variant)) == 1 & info(variant)$INDEL == F]

# 2. 粗mapping的时候，选则高质量的NA5_GT为0/0、1/1
GT <- geno(variant)$GT %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GT")) %>%
  filter(YE_L1_UDI164.MarkDup.bam.PosSorted.MarkDup.bam_GT == "0/0")
variant <- variant[rownames(variant) %in% rownames(GT)]

# 3. AD
# 提取 AD 数据
AD_raw <- geno(variant)$AD

# 清洗样本名（去除冗余路径/后缀）
sample_names <- colnames(AD_raw)
clean_names <- gsub(".*\\/|_dedup.*|_sorted.*|\\.bam", "", sample_names)

# 转为数据框并加上 sample_AD 后缀
AD_df <- as.data.frame(AD_raw)
colnames(AD_df) <- paste0(clean_names, "_AD")

# 转换每列为 “REF,ALT” 字符串
AD_df <- AD_df %>%
  mutate(across(everything(), ~ sapply(., function(x) paste(x, collapse = ","))))

# 拆分每个样本的 REF 和 ALT 到新列
for (sample in clean_names) {
  AD_df <- AD_df %>%
    separate(col = paste0(sample, "_AD"),
             into = c(paste0(sample, "_REF"), paste0(sample, "_ALT")),
             sep = ",", convert = TRUE)
}

# 转为数值
AD <- AD_df %>%
  mutate(across(everything(), as.numeric))
# 4. DP
DP <- geno(variant)$DP %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_DP"))

# 5. GQ
GQ <- geno(variant)$GQ %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GQ"))

# 6.
df <- reduce(list(GT, DP, AD, GQ), cbind) %>%
  mutate(QUAL = qual(variant), MQ = info(variant)$MQ) %>%
  rownames_to_column(var = "id") %>%
  separate(col = "id", into = c("Chr", "Pos", "Ref", "Alt"), sep = "[:_\\/]") %>%
  select(Chr:Alt, QUAL, MQ, contains("GT"), contains("GQ"), contains("YE"), contains("Ya10")) %>%
  filter((Ref == "G" & Alt == "A") | (Ref == "C" & Alt == "T")) %>%
  mutate(Pos = as.integer(Pos),
         YE_Index = `YE_L1_UDI164.MarkDup.PosSorted.MarkDup_ALT`/(YE_L1_UDI164.MarkDup.PosSorted.MarkDup_REF + YE_L1_UDI164.MarkDup.PosSorted.MarkDup_ALT),
         YA10_Index = Ya10_L7_UDI089_ALT/(Ya10_L7_UDI089_REF + Ya10_L7_UDI089_ALT),
         Delta_Index = YA10_Index - YE_Index)

# 7. 
gg_Pic <- df %>% filter(grepl("chr\\d", Chr)) %>% 
  ggplot() +
  geom_point(aes(x = Pos, y = Delta_Index, color = Chr), size = 1) +
  scale_x_continuous(labels = scales::unit_format(accuracy = 1, scale = 1e-6, unit = "")) +
  scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
  labs(x = NULL, y = "YA10 - YE") +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
  theme_bw(base_family = "Arial", base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.spacing = unit(2, "mm"),
        axis.text = element_text(color = "black", size = 9, family = "Arial"),
        axis.title.y = element_text(margin = margin(r = 7)),
        legend.position = "none")
gg_Pic

ggsave("1_Mapping.pdf", gg_Pic, width = 23, height = 8, units = "cm", device = cairo_pdf)
writexl::write_xlsx(df, "1_Mapping.xlsx")


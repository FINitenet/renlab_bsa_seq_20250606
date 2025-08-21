setwd("~/Documents/Project/20241212_TMM_hen1BSA/")

library(VariantAnnotation)
library(tidyverse)
library(ragg)

variant <- readVcf("20241212_TMM_hen1BSA/2_variants/variants_filtered.vcf.gz",
                   param = ScanVcfParam(info = c("INDEL", "MQ"),
                                        geno = c("GT", "PL", "DP", "AD", "GQ")))

# 1. 保留质量较好的，只有一种变异的SNP
variant <- variant[rowRanges(variant)$FILTER == "PASS" & elementNROWS(alt(variant)) == 1 & info(variant)$INDEL == F]

# 2. 粗mapping的时候，选则高质量的hen1_8_GT为0/0、1/1
GT <- geno(variant)$GT %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GT")) %>%
  filter(hen1_8_GT == "0/0")
variant <- variant[rownames(variant) %in% rownames(GT)]

# 3. AD
AD <- geno(variant)$AD %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_AD")) %>%
  mutate(hen1_8_AD = map(hen1_8_AD, paste0, collapse = ",") %>% unlist(),
         `2k_4116_9_AD` = map(`2k_4116_9_AD`, paste0, collapse = ",") %>% unlist()) %>%
  separate(col = hen1_8_AD, into = c("hen1_8_REF", "hen1_8_ALT"), sep = ",") %>%
  separate(col = `2k_4116_9_AD`, into = c("2k_4116_9_REF", "2k_4116_9_ALT"), sep = ",") %>%
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
  select(Chr:Alt, QUAL, MQ, contains("GT"), contains("GQ"), contains("hen1_8"), contains("2k_4116_9")) %>%
  filter((Ref == "G" & Alt == "A") | (Ref == "C" & Alt == "T")) %>%
  mutate(Pos = as.integer(Pos),
         hen1_8_Index = hen1_8_ALT/(hen1_8_REF + hen1_8_ALT),
         `2k_4116_9_Index` = `2k_4116_9_ALT`/(`2k_4116_9_REF` + `2k_4116_9_ALT`),
         Delta_Index = `2k_4116_9_Index` - hen1_8_Index)

# 7. 
gg_Pic <- ggplot(df %>% filter(grepl("Chr\\d", Chr))) +
  geom_point(aes(x = Pos, y = Delta_Index, color = Chr), size = 0.35) +
  scale_x_continuous(labels = scales::unit_format(accuracy = 1, scale = 1e-6, unit = "")) +
  scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
  labs(x = NULL, y = "2k_4116_9 - hen1_8") +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
  theme_bw(base_family = "Arial", base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.spacing = unit(2, "mm"),
        axis.text = element_text(color = "black", size = 9, family = "Arial"),
        axis.title.y = element_text(margin = margin(r = 7)),
        legend.position = "none")
ggsave("1_Mapping.pdf", gg_Pic, width = 23, height = 8, units = "cm", device = cairo_pdf)
writexl::write_xlsx(df, "1_Mapping.xlsx")


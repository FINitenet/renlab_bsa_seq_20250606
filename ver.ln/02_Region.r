# 本脚本的作用是根据粗mapping的结果，提取出Chr1 16M之后的SNP

setwd("/Users/sunlixue/Documents/Project/20241212_TMM_hen1BSA/")

library(VariantAnnotation)
library(tidyverse)
library(ragg)
library(writexl)

# 本函数的作用是输出可能的G-A、C-T突变位点
get_candidate <- function(vr){
  
  tmp_ref <- ref(vr) %>%
    as.character()
  tmp_alt <- alt(vr) %>%
    as.list()
  
  tmp_index <- c()
  for (i in 1:length(tmp_ref)) {
    tmp_alt2 <- tmp_alt[[i]] %>% as.character()
    if(("G" == tmp_ref[i]) & ("A" %in% tmp_alt2)){
      tmp_index <- c(tmp_index, i)
    } else if(("A" == tmp_ref[i]) & ("G" %in% tmp_alt2)){
      tmp_index <- c(tmp_index, i)
    } else if(("C" == tmp_ref[i]) & ("T" %in% tmp_alt2)){
      tmp_index <- c(tmp_index, i)
    } else if(("T" == tmp_ref[i]) & ("C" %in% tmp_alt2)){
      tmp_index <- c(tmp_index, i)
    }
  }
  
  tmp_GT <- geno(vr[tmp_index])$GT %>%
    as.data.frame() %>%
    rename_with(~ sub(".*\\/", "", .x)) %>%
    rename_with(~ sub("_dedup_realign.bam", "", .x))
  
  tmp_gr <- rowRanges(vr[tmp_index])
  for (j in colnames(tmp_GT)) {
    tmp_gr@elementMetadata[[j]] <- tmp_GT[,j]
  }
  
  return(list(vr[tmp_index], tmp_gr))
  
}

variant <- readVcf("variants_filtered.vcf.gz",
                   param = ScanVcfParam(info = c("INDEL", "DP", "MQ"),
                                        geno = c("GT", "PL", "DP", "AD", "GQ")))

# 0. 查看关键region中alt出现2的条目，如果有必要，单独输出到vcf
gr <- rowRanges(variant)
index <- seqnames(gr) == "Chr1" & start(gr) >= 16000000
variant2 <- variant[index]

variant2 <- variant2[elementNROWS(alt(variant2)) == 2]
res_list <- get_candidate(variant2)
res_list # 感觉里面没有可能的

# 1. 前面结果显示2种alt2的情况没有可能的，所以下面只保留1种SNP的
variant3 <- variant[index]
variant3 <- variant3[elementNROWS(alt(variant3)) == 1 & info(variant3)$INDEL == F]

# 2. 尽可能根据GT过滤一部分没必要注释的SNP（其实也可以直接输出所有variant，最后设置cutoff过滤可能的位点）
GT <- geno(variant3)$GT %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GT")) %>%
  filter(hen1_8_GT %in% c("0/0", "1/1")) %>%
  filter(`2k_4116_9_GT` != "./.") %>%
  filter(!(hen1_8_GT == "1/1" & `2k_4116_9_GT` == "1/1"))
variant3 <- variant3[rownames(variant3) %in% rownames(GT)]

seqlevels(variant3) <- c(1:5, "Mt", "Pt")
writeVcf(variant3, "2_Region.vcf")

# 理论上最可能的位点有：
# hen1_8    2k_4116_9
# 0/0     1/1
# 1/1     0/0

##--------------------------------------------------------------------------

# 整合信息输出

GT <- geno(variant3)$GT %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GT")) %>%
  rownames_to_column(var = "id") %>%
  separate(col = id, into = c("Chr", "Pos", "Ref", "Alt"), sep = "[:_/]", remove = F) %>%
  mutate(Pos = as.integer(Pos))

AD <- geno(variant3)$AD %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_AD")) %>%
  mutate(hen1_8_AD = map(hen1_8_AD, paste0, collapse = ",") %>% unlist(),
         `2k_4116_9_AD` = map(`2k_4116_9_AD`, paste0, collapse = ",") %>% unlist()) %>%
  separate(col = hen1_8_AD, into = c("hen1_8_REF", "hen1_8_ALT"), sep = ",") %>%
  separate(col = `2k_4116_9_AD`, into = c("2k_4116_9_REF", "2k_4116_9_ALT"), sep = ",") %>%
  mutate(across(everything(), as.numeric))

DP <- geno(variant3)$DP %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_DP"))

GQ <- geno(variant3)$GQ %>%
  as.data.frame() %>%
  rename_with(~ gsub(".*\\/|_dedup.*|_sorted.*", "", .x)) %>%
  rename_with(~ paste0(.x, "_GQ"))

df <- reduce(list(GT, DP, AD, GQ), cbind) %>%
  mutate(QUAL = qual(variant3), MQ = info(variant3)$MQ) %>%
  select(id, Chr:Alt, QUAL, MQ, contains("GT"), contains("GQ"), contains("hen1_8"), contains("2k_4116_9")) %>%
  mutate(Pos = as.integer(Pos),
         hen1_8_Index = hen1_8_ALT/(hen1_8_REF + hen1_8_ALT),
         `2k_4116_9_Index` = `2k_4116_9_ALT`/(`2k_4116_9_REF` + `2k_4116_9_ALT`),
         Delta_Index = `2k_4116_9_Index` - hen1_8_Index)
write_xlsx(df, "2_Region.xlsx")


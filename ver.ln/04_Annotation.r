setwd("/Users/sunlixue/Documents/Project/20241212_TMM_hen1BSA/")

library(readxl)
library(tidyverse)
library(writexl)

# 本函數的作用是整合gene description
get_gene_des <- function(gene_list){
  
  genes <- str_split(gene_list, pattern = ";") %>%
    unlist() %>%
    str_trim()
  
  if(length(genes) == 0){ # important
    return("")
  }
  
  tmp_anno <- gene_anno %>%
    filter(name %in% genes)
  
  if(nrow(tmp_anno) == 0){ # important
    return("")
  }
  
  tmp_anno <- tmp_anno %>%
    mutate(des = paste0("@", name, " : ", symbol, " : ", description))
  tmp_des <- paste(tmp_anno$des, collapse = " \n ")
  
  return(tmp_des)
}

# 1. 整合SNP信息和Annovar注释结果
snp <- read_excel("2_Region.xlsx")
annovar <- read_tsv("3_Annovar.AT_multianno.txt", show_col_types = F) %>%
  select(id = Otherinfo6, Filter = Otherinfo10, Func.refGene:AAChange.refGene)

df <- left_join(snp, annovar) %>%
  select(-id)

# 2. 再整合基因注释信息
gene_anno <- read_excel("/Users/sunlixue/Documents/database/TAIR_Data_20200630/20200630_Alias_Description.xlsx") %>%
  select(name = ID, symbol, description = Description)
gene_desc <- map(df$Gene.refGene, get_gene_des) %>%
  unlist()
df$Gene_Annotation <- gene_desc

# 3. 输出
write_xlsx(df, "4_Annotation.xlsx")

# hen1_8  2k_4116_9
# 0/0    1/1
# 1/1    0/0

# 4. 最可能的
a <- df %>%
  filter(hen1_8_Index <= 0.2, `2k_4116_9_Index` >= 0.8) %>%
  filter((Ref == "G" & Alt == "A") | (Ref == "C" & Alt == "T"))
b <- df %>%
  filter(hen1_8_Index >= 0.8, `2k_4116_9_Index` <= 0.2) %>%
  filter((Ref == "A" & Alt == "G") | (Ref == "T" & Alt == "C"))
write_xlsx(rbind(a, b), "4_Annotation_Candidate.xlsx")


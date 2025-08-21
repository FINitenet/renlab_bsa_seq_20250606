# devtools::install_github("bmansfeld/QTLseqr")
library("QTLseqr")


#Set sample and file names
HighBulk <- "Ya10_L7_UDI089"
LowBulk <- "Na10_L7_UDI096"
file <- "YA10_raw.filtered.table"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- paste0(rep("chr", 5), 1:5)

#Import SNP data from file
df <-
  importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 10,
    maxTotalDepth = 10000,
    minSampleDepth = 10,
    minGQ = 99, 
    verbose = TRUE
  )


df_filt <- df_filt %>%
  filter((REF == "G" & ALT == "A") | (REF == "C" & ALT == "T"))


#Run G' analysis
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 1,
  popStruc = "F2",
  bulkSize = c(25, 25),
  replications = 10000,
  intervals = c(95, 99)
)

#Plot

ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,1000)

ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH))

plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")

plotQTLStats(SNPset = df_filt, var = "nSNPs")


plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

#export summary CSV
getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")



ggplot(df_filt, aes(x = POS, y = deltaSNP)) +
  geom_point(aes(color = CHROM)) +
  facet_wrap(~CHROM, scales = "free_x") +
  theme_bw() +
  labs(y = "deltaSNP", x = "Genomic Position")

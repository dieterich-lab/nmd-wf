#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DRIMSeq)
  library(DEXSeq)
  library(BiocParallel)
  library(stringr)
  library(tidyr)
})

d <- readRDS("data/drimseq_data.RDS")

BPPARAM <- SnowParam(
  workers = 20,
  type = "SOCK",
  progressbar = TRUE,
  RNGseed = 123,
  log = TRUE)

dxd <- DEXSeqDataSet(
  countData = round(counts(d)[-c(1, 2)]),
  sampleData = DRIMSeq::samples(d),
  design = ~ sample + exon + group:exon,
  featureID = counts(d)$feature_id,
  groupID = counts(d)$gene_id
)
# dxd <- estimateSizeFactors(dxd)
# dxd <- estimateDispersions(dxd, quiet = FALSE)
fullModel <- ~ sample + exon + group:exon
reducedModel <- ~ sample + exon

contrast <- c(
  'HEK_SMG5-KD_Z023'='HEK_NoKO_SMG5KD-KD_Z023',
  'HEK_Luc-KD_Z023'='HEK_NoKO_LucKD-KD_Z023',
  'HEK_SMG6+SMG7-KD_Z023'='HEK_NoKO_SMG6+SMG7KD-KD_Z023',
  'HEK_Luc-KD_Z023'='HEK_NoKO_LucKD-KD_Z023',
  'HEK_SMG7-KO_SMG5-KD_Z245'='HEK_SMG7KO_SMG5KD-KD_Z245',
  'HEK_SMG7-KO_Luc-KD_Z245'='HEK_SMG7KO_LucKD-KD_Z245',
  'HEK_SMG7-KO_SMG6-KD_Z245'='HEK_SMG7KO_SMG6KD-KD_Z245',
  'HEK_SMG7KO_LucKD-KD_Z245'='HEK_SMG7KO_LucKD-KD_Z245',
  'HEK_SMG7-KO_SMG5-KD_Z319'='HEK_SMG7KO_SMG5KD-KD_Z319',
  'HEK_SMG7-KO_Luc-KD_Z319'='HEK_SMG7KO_LucKD-KD_Z319',
  'HEK_SMG7-KO_SMG6-KD_Z319'='HEK_SMG7KO_SMG6KD-KD_Z319',
  'HEK_SMG7-KO_Luc-KD_Z319'='HEK_SMG7KO_LucKD-KD_Z319',
  'HeLa_SMG6+SMG7-KD_Z021'='HeLa_NoKO_SMG6+SMG7KD-KD_Z021',
  'HeLa_Luc-KD_Z021'='HeLa_NoKO_LucKD-KD_Z021',
  'MCF7_SMG6+SMG7-KD'='MCF7_NoKO_SMG6+SMG7KD-KD',
  'MCF7_Luc-KD'='MCF7_NoKO_LucKD-KD',
  'U2OS_SMG6+SMG7-KD'='U2OS_NoKO_SMG6+SMG7KD-KD',
  'U2OS_Luc-KD'='U2OS_NoKO_LucKD-KD'
)
levels(colData(dxd)$group) <- contrast[levels(colData(dxd)$group)]
contrast <- as.data.frame(matrix(contrast, ncol=2, byrow = T))
dexeq_results <- function(dexseq_obj, num, denom) {
  message(str_glue("Computing DEU for {num} vs {denom}"))
  dexseq_obj <- dexseq_obj[, colData(dexseq_obj)$group %in% c(num, denom)]
  dexseq_obj <- testForDEU(dexseq_obj, fullModel = fullModel, reducedModel = reducedModel)
  dexseq_obj <- estimateExonFoldChanges(dexseq_obj, fitExpToVar = "group")
  DEXSeqResults(dexseq_obj) %>%
    data.frame()
}

res <- apply(contrast, 1, function(x) {
  dexeq_results(dxd, x[1], x[2])
})
names(res) <- paste0(contrast[, 1], '-vs-', contrast[, 2])
res <- bind_rows(res, .id = "contrasts") %>% 
  rename(
    "gene_id" = groupID, 
    "transcript_id" = featureID)
rownames(res) <- NULL
saveRDS(res, "data/dte_results.RDS")
saveRDS(dxd, "data/dexseq_data.RDS")


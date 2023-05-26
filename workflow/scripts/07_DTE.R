#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DRIMSeq)
  library(DEXSeq)
  library(here)
  library(openxlsx)
  library(BiocParallel)
})

d <- readRDS("data/drimseq_data.RDS")

BPPARAM <- SnowParam(
  workers = 10,
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
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE, BPPARAM=BPPARAM)
fullModel <- ~ sample + exon + group:exon
reducedModel <- ~ sample + exon

dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges(dxd,fitExpToVar='group',  BPPARAM=BPPARAM)
saveRDS(dxd,  "data/dte_results.RDS")
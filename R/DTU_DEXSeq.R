suppressPackageStartupMessages({
  library(DRIMSeq)
  library(DEXSeq)
  library(here)
  library(openxlsx)
})

d <- readRDS(here("phase2", "data", "tx_counts.RDS"))

dxd <- DEXSeqDataSet(
  countData = round(counts(d)[-c(1, 2)]),
  sampleData = DRIMSeq::samples(d),
  design = ~ sample + exon + group:exon,
  featureID = counts(d)$feature_id,
  groupID = counts(d)$gene_id
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
fullModel <- ~ sample + exon + group:exon
reducedModel <- ~ sample + exon

#save.image('DTU_DEXSeq.Rdata')
#load('phase2/DTU_DEXSeq.Rdata')

ref <- rtracklayer::import(here("../DFG_seq_Nanopore/GRCh38_90_SIRV_Set3.gtf"))
ref <- ref %>% 
  mcols() %>% 
  as.data.frame() %>% 
  filter(type == 'transcript') %>% 
  dplyr::select('gene_id', 'gene_name', 'transcript_id', "transcript_name") 

dirloc.out <- here('phase2/results/')

dexseq_results <- function(dexseq_obj, num, denom) { 
  dexseq_obj <- dexseq_obj[,colData(dexseq_obj)$group %in% c(num, denom)]
  dexseq_obj <- testForDEU(dexseq_obj, fullModel = fullModel, reducedModel = reducedModel)
  dexseq_obj <- estimateExonFoldChanges(dexseq_obj, fitExpToVar = "group")
  res <- as.data.frame(DEXSeqResults(dexseq_obj))
  
  wb <- createWorkbook()
  sheet_name <- str_trunc(str_glue("{num}_vs_{denom}"), 30, "right")
  addWorksheet(wb, sheetName=sheet_name)
  writeDataTable(wb, sheet=1, x=res)
  saveWorkbook(wb, str_glue("{dirloc.out}/dtu_{num}_vs_{denom}.xlsx"), overwrite=TRUE)
}
contrasts <- data.frame(
  num=c('SMG5KD', 'SMG7KO_2', "SMG7KO_2_SMG5KD", "SMG7KO_2_SMG6KD", "SMG7KO_34",
        "SMG7KO_34_SMG5KD", "SMG7KO_34_SMG6KD",
        "HEK_SMG67", "HeLa_SMG67", "U2OS_SMG67", "MCF7_SMG67", "MCF7_7KO_SMG67"),
  denom=c(rep("control", 7), 
          "HEK_control", "HeLa_control", "U2OS_control", "MCF7_control", 
          "MCF7_7KO_control")) 

results <- apply(contrasts, 1, function(x) {
  dexseq_results(dxd, as.character(x[1]), as.character(x[2]))
})




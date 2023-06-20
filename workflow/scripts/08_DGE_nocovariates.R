suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(here)
  library(IHW)
  library(org.Hs.eg.db)
})

output <- 'dge_results/'

dds <- readRDS(here('phase2', 'data', 'dge_dds.RDS'))
message(nrow(dds))
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
message(nrow(dds))
dds <- DESeq(dds)

contrasts <- data.frame(
  treat=c(
    "HEK_SMG5-KD_Z023",
    "HEK_SMG7-KO_SMG5-KD_Z245",
    "HEK_SMG7-KO_SMG6-KD_Z245",
    "HEK_SMG7-KO_SMG5-KD_Z319",
    "HEK_SMG7-KO_SMG6-KD_Z319",
    "HEK_SMG6+SMG7-KD_Z023",
    "HeLa_SMG6+SMG7-KD_Z021",
    "U2OS_SMG6+SMG7-KD",
    "MCF7_SMG6+SMG7-KD"),
  control=c(
    "HEK_Luc-KD_Z023",
    "HEK_Luc-KD_Z023",
    "HEK_Luc-KD_Z023",
    "HEK_SMG7-KO_Luc-KD_Z319",
    "HEK_SMG7-KO_Luc-KD_Z319",
    "HEK_Luc-KD_Z023",
    "HeLa_Luc-KD_Z021",
    "U2OS_Luc-KD",
    "MCF7_Luc-KD"
  )) 


results <- apply(contrasts, 1, function(x) {
  deseq_results(dds, x[1], x[2])
})
names(results) <- paste0(contrasts$treat, '-vs-', contrasts$control)
results <- bind_rows(results, .id='contrast') 

results %<>% 
  dplyr::rename(
    gene_id_tmp=gene_id, 
    gene_id=gene_name) %>% 
  dplyr::rename(
    gene_name=gene_id_tmp)

write_rds(results, here('phase2/results/dge_results.rds'))


# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, sheetName='contrasts')
# writeDataTable(wb, sheet = 'contrasts', as.data.frame(contrasts))
#sheet=as.character(length(names(wb)) + 1)
#addWorksheet(wb, sheetName=sheet)
# session.info <- capture.output(sessionInfo())
# addWorksheet(wb, 'session info')
# writeDataTable(wb, sheet = 'session info', as.data.frame(session.info))
# 
# saveWorkbook(wb, here('phase2/results/dtu_results.xlsx'), overwrite=TRUE)

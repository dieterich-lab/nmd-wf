suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(here)
  library(IHW)
  library(org.Hs.eg.db)
  library(openxlsx)
})


output <- 'dge_results/'
lfcThreshold <- 1
altHypothesis <- "greaterAbs"
alpha <- 0.05

deseq_results <- function (dds, num, denom) { 
  
  res <- results(dds,
                 contrast=c("group", num, denom),
                 lfcThreshold=1, 
                 altHypothesis=altHypothesis,
                 alpha=alpha,
                 filterFun=ihw)
  
  res$SYMBOL <- mapIds(org.Hs.eg.db,
                       keys=rownames(res),
                       column="ENSEMBL",
                       keytype="SYMBOL",
                       multiVals="first")
  
  res.tib <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()

  wb <- createWorkbook()
  sheet_name <- str_trunc(str_glue("{num}_vs_{denom}"), 30, "right")
  addWorksheet(wb, sheetName=sheet_name)
  writeDataTable(wb, sheet=1, x=res.tib)

  saveWorkbook(wb, str_glue("{output}/{num}_vs_{denom}.xlsx"), overwrite=TRUE)
  
  res
}


dds <- readRDS(here('phase2', 'data', 'dge_dds.RDS'))
message(nrow(dds))
group_counts <- t(rowsum(t(counts(dds)), dds$group))
keep <- rowMin(group_counts) > 10
dds <- dds[keep,]
message(nrow(dds))
dds$group <- factor(gsub(x=dds$group, " ", ""))
dds <- DESeq(dds)

contrasts <- data.frame(
  num=c('SMG5KD', 'SMG7KO_2', "SMG7KO_2_SMG5KD", "SMG7KO_2_SMG6KD", "SMG7KO_34",
          "SMG7KO_34_SMG5KD", "SMG7KO_34_SMG6KD",
          "HEK_SMG67", "HeLa_SMG67", "U2OS_SMG67", "MCF7_SMG67", "MCF7_7KO_SMG67"),
  denom=c(rep("control", 7), 
          "HEK_control", "HeLa_control", "U2OS_control", "MCF7_control", 
          "MCF7_7KO_control")) 

results <- apply(contrasts, 1, function(x) {
  deseq_results(dds, as.character(x[1]), as.character(x[2]))
})

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(DESeq2)
  library(IHW)
})

snakemake@source('utils.R')
set_here(path = "/prj/Niels_Gehring/nmd_transcriptome/phaseFinal")
message(getwd())

message('Importing abundances')

txi.genes <- readRDS(here("data", "gene_counts.RDS"))
tx2gene <- readRDS(here("data", "tx2gene.RDS"))
metadata <- readRDS(here("data", "metadata.RDS"))
metadata <- as.data.frame(metadata)
metadata$Knockout <- replace_na(metadata$Knockout, '_NoKO')
metadata$Knockout <- str_replace(metadata$Knockout, '-', '')
metadata$Knockdown <- str_c(metadata$Knockdown, 'KD')
metadata$group <- metadata %>% 
  dplyr::select(Knockdown, Knockout) %>%
  unite(group, sep = '') %>% 
  pull(group)
metadata$group = as.factor(metadata$group)
metadata$cellline = as.factor(metadata$cellline)

dds <- DESeqDataSetFromTximport(
  txi = txi.genes,
  colData = metadata[, c('group', 'cellline')],
  design = ~cellline+group)

message(nrow(dds))
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
message(nrow(dds))
dds <- DESeq(dds)

contrasts <- data.frame(
  treat=c(
    "SMG5KD_NoKO",
    "SMG5KD_SMG7KO",
    "SMG6KD_SMG7KO",
    "SMG6+SMG7KD_NoKO"
  ),
  control=c(
    "LucKD_NoKO",
    "LucKD_NoKO",
    "LucKD_NoKO",
    "LucKD_NoKO"
  )) 

results <- apply(contrasts, 1, function(x) {
  deseq_results(dds, x[1], x[2])
})
names(results) <- paste0(contrasts$treat, '-vs-', contrasts$control)
results <- bind_rows(results, .id='contrast') 

tx2gene <- tx2gene %>% select(gene_id, gene_name) %>% filter(!is.na(gene_name)) %>% distinct()
results <- merge(results, tx2gene, all.x=TRUE)
saveRDS(results,  "data/dge_results.RDS")
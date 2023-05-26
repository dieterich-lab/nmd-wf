suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(DESeq2)
  library(IHW)
  library(PCAtools)
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

vst <- assay(vst(dds))
p <- pca(vst, metadata = colData(dds), removeVar = 0.1)
suppressWarnings({horn <- parallelPCA(vst)}) 
screeplot(
  p,
  components = getComponents(p, 1:30),
  vline = c(horn$n)) +
  geom_label(
    aes(x = horn$n + 1, y = 50, label = 'Horn\'s', vjust = -1, size = 8))

pdf("results/dge_pca.pdf", height = 10, width = 20)
biplot(
  p, 
  showLoadings = FALSE,
  lab = NULL,
  pointSize = 3, 
  colby = 'group', 
  shape = 'cellline',
  legendPosition = 'right', 
  hline = 0, 
  vline = 0) 
dev.off()

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

message("Number of genes:", nrow(dds))
dds <- DESeq(dds)

results <- apply(contrasts, 1, function(x) {
  deseq_results(dds, x[1], x[2])
})
names(results) <- paste0(contrasts$treat, '-vs-', contrasts$control)
results <- bind_rows(results, .id='contrast') 

tx2gene <- tx2gene %>% select(gene_id, gene_name) %>% filter(!is.na(gene_name)) %>% distinct()
results <- merge(results, tx2gene, all.x=TRUE)
saveRDS(results,  "data/dge_results.RDS")
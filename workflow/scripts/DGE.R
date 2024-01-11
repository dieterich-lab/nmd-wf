# see /beegfs/prj/Niels_Gehring/nmd_transcriptome/dge/deseq2.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(here)
  library(IHW)
  library(org.Hs.eg.db)
})


set_here(path = "/prj/Niels_Gehring/nmd_transcriptome/phaseFinal")
source('workflow/scripts/utils.R')
message(getwd())

message('Importing abundances')

txi.genes <- readRDS(here("data", "gene_counts.RDS"))
tx2gene <- readRDS(here("data", "tx2gene.RDS"))
metadata <- readRDS("data/metadata.RDS")
metadata <- as.data.frame(metadata)
metadata$Knockout <- replace_na(metadata$Knockout, '_NoKO')
metadata$Knockout <- str_replace(metadata$Knockout, '-', '')
metadata$Knockdown <- str_c(metadata$Knockdown, 'KD')
metadata$group <- metadata %>% 
  dplyr::select(Knockdown, Knockout) %>%
  unite(group, sep = '') %>% 
  pull(group)
metadata$cellline <- word(metadata$Cell_line, 1)
metadata$group <- str_glue_data(metadata, "{cellline}{Knockout}_{Knockdown}-KD_{clone}", .na='')
metadata$group <- gsub(x=metadata$group, "_$", '')
metadata$group <- as.factor(metadata$group)

dds <- DESeqDataSetFromTximport(
  txi = txi.genes,
  colData = metadata |> select('group'),
  design = ~group)


# dput(x, control = NULL)
x <- data.frame(
  group = c(
    'HEK_NoKO_SMG5KD-KD_Z023',
    'HEK_NoKO_LucKD-KD_Z023',
    'HEK_NoKO_SMG6+SMG7KD-KD_Z023',
    'HEK_NoKO_LucKD-KD_Z023',
    'HEK_SMG7KO_SMG5KD-KD_Z245',
    'HEK_SMG7KO_LucKD-KD_Z245',
    'HEK_SMG7KO_SMG6KD-KD_Z245',
    'HEK_SMG7KO_LucKD-KD_Z245',
    'HEK_SMG7KO_SMG5KD-KD_Z319',
    'HEK_SMG7KO_LucKD-KD_Z319',
    'HEK_SMG7KO_SMG6KD-KD_Z319',
    'HEK_SMG7KO_LucKD-KD_Z319',
    'HeLa_NoKO_SMG6+SMG7KD-KD_Z021',
    'HeLa_NoKO_LucKD-KD_Z021',
    'MCF7_NoKO_SMG6+SMG7KD-KD',
    'MCF7_NoKO_LucKD-KD',
    'U2OS_NoKO_SMG6+SMG7KD-KD',
    'U2OS_NoKO_LucKD-KD')
)

contrasts <- matrix(x$group, ncol =2 ,byrow = T)
x$contrasts <- paste0(contrasts[, 1], '-vs-', contrasts[, 2]) %>% rep(each=2)


message("Number of genes:", nrow(dds))

dds <- DESeq(dds)

results <- apply(contrasts, 1, function(x) {
  deseq_results(dds, x[1], x[2])
})
names(results) <- paste0(contrasts$treat, '-vs-', contrasts$control)
results <- bind_rows(results, .id='contrast') 

tx2gene <- tx2gene %>% select(gene_id, gene_name) %>% distinct()
results <- merge(results, tx2gene, all.x=TRUE)
saveRDS(results,  "data/dge_results.RDS")
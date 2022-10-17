library(tximport)
library(here)
library(tidyverse)
library(GenomicRanges)
message('Import abundances')

ref <- rtracklayer::import(here("phase2/stringtie_merge/merged_each.gtf"))

tx2gene <- ref %>%
  subset(type == 'transcript' ) %>%
  mcols() %>% 
  as.data.frame() %>% 
  filter(!is.na(gene_name)) %>% 
  dplyr::select(transcript_id, gene_name)

files <- Sys.glob(here("phase2/salmon/*/quant.sf"))
names(files) <- str_split(files, "/", simplify = TRUE)[, 8]

metadata <- read.csv(here('phase2/config/metadata_w_files.csv'))
metadata$files <- files[match(metadata$CCG_Sample_ID, names(files))]
metadata <- metadata[!is.na(metadata$files), ]
metadata$group <- gsub(x=metadata$Condition, " ", "")
metadata$cellline <- word(metadata$Cell_line, 1)

txi <- tximport(metadata$files, type = "salmon", tx2gene = tx2gene)

library(DESeq2)
message('Import create deseq object')

dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = metadata['group'],
  design = ~group)

saveRDS(dds, here('phase2', 'data', 'dge_dds.RDS'))
message('done')
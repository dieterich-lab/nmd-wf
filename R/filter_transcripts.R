suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(here)
  library(GenomicRanges)
  library(yaml)
  library(DRIMSeq)
})

set_here(path = "/prj/Niels_Gehring/nmd_transcriptome")

message('Import abundances')

ref <- rtracklayer::import(here("phase2/stringtie_merge/merged_mix.gtf"))

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

txi <- tximport(metadata$files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)

message("Start DRIMSeq filter")

samples <- data.frame(
  sample_id = make.unique(metadata$group, sep = '_'),
  group = metadata$group
)

cts <- as.data.frame(txi$counts)
colnames(cts) <- make.unique(metadata$group, sep = '_')
cts$feature_id <- rownames(cts)
cts$gene_id <- tx2gene$gene_name[match(cts$feature_id, tx2gene$transcript_id)]
cts <- cts[!is.na(cts$gene_id), ]

d <- DRIMSeq::dmDSdata(counts = cts, samples = samples)

n <- 3
n_small <- 3
message("Filtering dataset data")
message("Total number of tx: ", sum(elementNROWS(d@counts)))
message("Total number of genes: ", length(d@counts))
d <- dmFilter(
  d,
  min_samps_feature_prop = n_small,
  min_feature_prop = 0.1,
  min_samps_feature_expr = n_small,
  min_feature_expr = 10,
  min_samps_gene_expr = n,
  min_gene_expr = 10,
  run_gene_twice = TRUE
)
message("Filtered dataset data")
message("Total number of tx: ", sum(elementNROWS(d@counts)))
message("Total number of genes: ", length(d@counts))


dir.create(here("phase2", "data"))
saveRDS(d, here("phase2", "data", "tx_counts.RDS"))
saveRDS(samples, here("phase2", "data", "samples.RDS"))
saveRDS(tx2gene, here("phase2", "data", "tx2gene.RDS"))


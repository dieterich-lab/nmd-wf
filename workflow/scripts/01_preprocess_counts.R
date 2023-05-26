suppressPackageStartupMessages({
  library(tximport)
  library(here)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
})

snakemake@source('utils.R')
set_here(path = "/prj/Niels_Gehring/nmd_transcriptome/phaseFinal")
message(getwd())

ref <- rtracklayer::import(here("stringtie_merge/merged_each.fix.gtf"))

tx2gene <- ref %>%
  subset(type == 'transcript' ) %>%
  mcols() %>% 
  as.data.frame() %>% 
  dplyr::select(transcript_id, gene_id, gene_name, ref_gene_id)

files <- Sys.glob(here("salmon/*/quant.sf"))
names(files) <- str_split(files, "/", simplify = TRUE)[, 8]

metadata <- read.csv(here('config/metadata_w_files.csv'))
metadata <- filter(metadata, !CCG_Sample_ID %in% c(
  "A006200074_120335", 
  "A006200074_120337",
  "A006200074_120339",
  "A006200074_120341",
  "A006200074_120343",
  "A006200074_120345"
  ))
metadata$files <- files[match(metadata$CCG_Sample_ID, names(files))]
metadata <- metadata[!is.na(metadata$files), ]

message('Import transcript abundances')

txi <- tximport(
  metadata$files, 
  type = "salmon", 
  tx2gene = tx2gene, 
  txOut = TRUE, 
  countsFromAbundance="dtuScaledTPM")

colnames(txi$counts) <- make_groups_unique(metadata$group)

cts <- as.data.frame(txi$counts)

cts[cts$gene_id == 'ENSG00000161547']

cts$feature_id <- rownames(cts)
cts$gene_id <- tx2gene$gene_id[match(cts$feature_id, tx2gene$transcript_id)]
cts <- cts[!is.na(cts$gene_id), ]

message("Start DRIMSeq filter")
samples <-  data.frame(
  sample_id = colnames(txi$counts),
  group = str_sub(colnames(txi$counts), 1, -3)
)

d <- DRIMSeq::dmDSdata(counts = cts, samples = samples)

n <- 3
n_small <- 3
message("Filtering dataset data")
message("Total number of tx: ", sum(elementNROWS(d@counts)))
message("Total number of genes: ", length(d@counts))
d <- DRIMSeq::dmFilter(
  d,
  min_samps_feature_prop = n_small,
  min_feature_prop = 0.1,
  min_samps_feature_expr = n_small,
  min_feature_expr = 1,
  min_samps_gene_expr = n,
  min_gene_expr = 10,
  run_gene_twice = TRUE
)
message("Filtered dataset data")
message("Total number of tx: ", sum(elementNROWS(d@counts)))
message("Total number of genes: ", length(d@counts))

dir.create(here("data"))

tx2gene$keep <- tx2gene$transcript_id %in% rownames(d@counts@unlistData)

message('Import gene abundances')
txi.genes <- tximport(
  metadata$files, 
  type = "salmon", 
  tx2gene = tx2gene[tx2gene$keep, ], 
  txOut = FALSE,
  countsFromAbundance="lengthScaledTPM")
colnames(txi.genes$counts) <- make_groups_unique(metadata$group)

saveRDS(d, here("data", "drimseq_data.RDS"))
saveRDS(ref, here("data", "gtf.RDS"))
saveRDS(txi.genes, here("data", "gene_counts.RDS"))
saveRDS(txi, here("data", "tx_counts.RDS"))
saveRDS(samples, here("data", "samples.RDS"))
saveRDS(tx2gene, here("data", "tx2gene.RDS"))
saveRDS(metadata, here("data", "metadata.RDS"))

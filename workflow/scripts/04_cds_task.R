#! /usr/bin/env Rscript
## ---------------------------
## Script name: cds_task.R
## Purpose of script: Combined CDS from multiple sources and 
## annotated the transcript PTC status
## Author: Thiago Britto-Borges
## Date Created: 2023-05-22
## Copyright (c) Thiago Britto-Borges, DieterichLab 2023
## Email: thiago.brittoborges@uni-heidelberg.de
suppressPackageStartupMessages({
  snakemake@source('utils.R')
  
  library(plyranges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(here)
  library(dplyr)
})

# Read reference annotation, CDS sequence info, genes,and combined CDS data
tx2gene <- readRDS("data/tx2gene.RDS")
reference_annotation <- readRDS("phaseFinal/data/gtf.RDS")
cds_sequence_info <- readRDS("phaseFinal/data/cds_seq.RDS")
combined_cds <- readRDS('phaseFinal/data/combined_cds.RDS')

tx2gene <- tx2gene[tx2gene$keep, ]
reference_annotation <- reference_annotation %>% filter(transcript_id %in% tx2gene$transcript_id) 
cds_sequence_info <- cds_sequence_info %>% subset(transcript %in% tx2gene$transcript_id)

# Set levels of 'source' column and sort combined CDS data
combined_cds$source <- factor(combined_cds$source, levels = c(
  "ensembl", "hek293gao", "ribotish", "openprot", "transdecoder"))
combined_cds <- sort(combined_cds, by = ~source, decreasing = FALSE, ignore.strand = FALSE)

# Filter exons from reference annotation and select transcript IDs matching CDS sequence info
exons <- reference_annotation %>%
  filter(type == 'exon') %>% 
  plyranges::select(transcript_id) %>%
  filter(transcript_id %in% cds_sequence_info$transcript)

# Split the filtered exons by transcript ID
exons_by_transcript <- split(exons, ~transcript_id)
# Map exons coordinates to cDNA 
exons_t <- pmapToTranscripts(
  exons,
  exons_by_transcript[exons$transcript_id]
)
exons_t <- sort(exons_t)

cds_granges <- GRanges(cds_sequence_info$transcript, cds_sequence_info$ranges)
cds_granges <- mapFromTranscripts(
  cds_granges, 
  exons_by_transcript)
cds_granges2 <- pintersect(
  exons_by_transcript[cds_granges$transcriptsHits],
  cds_granges,
  drop.nohit.ranges=TRUE)
names(cds_granges2) <- names(cds_granges)
cds_granges <- cds_granges2
rm(cds_granges2)

cds_tranges <- GRanges(
  cds_sequence_info$transcript, 
  cds_sequence_info$ranges)
cds_tranges$cds_id <- names(cds_tranges)

cds_tranges2 <- findOverlapPairs(
  cds_tranges,
  exons_t,
  ignore.strand=TRUE)
cds_tranges2 <- pintersect(cds_tranges2)

exons_t <- split(ranges(exons_t), seqnames(exons_t))
stopifnot(all(all(exons_t %>% heads(1) %>% start() == 1)))
cds_tranges <- split(ranges(cds_tranges2), names(cds_tranges2))
rm(cds_tranges2)

cds_sequence_info$cds_granges <- cds_granges[cds_sequence_info$cds_id]
cds_sequence_info$cds_tranges <- cds_tranges[cds_sequence_info$cds_id]


# Filter out transcripts with unique strand and seqnames
filtered_exons_by_transcript <- exons_by_transcript[
  lengths(runValue(strand(exons_by_transcript))) == 1 &
    lengths(runValue(seqnames(exons_by_transcript))) == 1]

# Ensure the number of filtered exons and exons with unique strand and seqnames match
stopifnot(length(exons_by_transcript) == length(filtered_exons_by_transcript))

# Check if the total number of exons in filtered exons matches the number of exons in the original dataset
stopifnot(sum(lengths(filtered_exons_by_transcript)) == length(exons))

# Check if all transcript IDs in CDS sequence info are present in the filtered exons
stopifnot(isEmpty(setdiff(cds_sequence_info$transcript, names(filtered_exons_by_transcript))))
rm(filtered_exons_by_transcript)

# Anchor the all CD  to their start codon and add info as metadata

# we dont need this column anymore because we have cds_tranges (spliced and ungapped)
cds_sequence_info$ranges <- NULL
query_granges <- cds_sequence_info$cds_granges %>% range() %>% unlist() %>% anchor_5p() %>% mutate(width = 3)
mcols(query_granges) <- cds_sequence_info

# Anchor the CDS from sources to their start codon and annotatea with gene names
subject_granges <- combined_cds %>% anchor_5p() %>% mutate(width = 3)

# Perform join between CDS from sources and all CDS with a minimum overlap of 3
result <- join_overlap_intersect_directed(query_granges, subject_granges, minoverlap = 3)
result <- mcols(result)

stopifnot(all(sum(width(result$cds_granges)) == width(result$cds_seq)))

result$leej <- exons_t[result$transcript] %>% tails(1) %>% start() - 3
result$stop_pos <- result$cds_tranges %>% tails(1) %>% end() %>% unlist2()
# +3 because the CDS does not include the stop codon
result$ptc <- unlist(result$leej) - unlist(result$stop_pos) > 50 + 3

saveRDS(result, snakemake@output[[1]])
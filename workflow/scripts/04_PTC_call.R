#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  
})


if (length(commandArgs(trailingOnly=TRUE)) != 2) {
  cat("Usage: Rscript script.R input_file output_file\n")
  q(status = 1)
}

ungap <- function(.x) {
  IRanges(end=cumsum(width(.x)), width = width(.x))
}

ungap_keep_start <- function(.x) {
  shift(ungap(.x), start(.x)[1] -1)
}

dist_to_stop <- function(x) {
  leej <- x$c_exons |>
    tails(1) |>
    start() - 3 |>
    unlist()
  
  stop_pos <- x$c_thick |>
    tails(1) |>
    end()
  empty_stop <- lapply(stop_pos, isEmpty) |> unlist()
  stop_pos[empty_stop] <- leej[empty_stop]
  stop_pos <- stop_pos |> unlist()
  # +3 because the CDS does not include the stop codon
  unlist(leej) - stop_pos
}


#' Remap thick and blocks to genomic, transcript and cDNA coord
#'
#' This function takes a GRanges object and extracts thick exons by intersecting them
#' with a set of transcript blocks. It also applies a gap removal function to both
#' the original blocks and the thick exons.
#' g = genomic coordinate
#' t = transcript coordinate 
#' c = cDNA coordinate (i.e., transcript without introns)
#'
#' @param x A GRanges object representing bed12 transcriptome with thick column.
#' @return A modified GRanges object with additional columns
#'
#' @import GenomicRanges
#' @export
map_bed12_features <- function(.x) {
  times <- lengths(.x$blocks)
  blocks <- unlist(.x$blocks)
  blocks <- GRanges(
    rep(.x$name, times=times),
    blocks
  )
  transcripts <- .x
  mcols(transcripts) <- NULL
  names(transcripts) <- .x$name
  
  # bed12 blocks are not stranded, but we need the strandness
  g_exons <- GenomicFeatures::mapFromTranscripts(
    blocks,
    transcripts, 
    ignore.strand=TRUE)
  strand(g_exons) <- strand(.x)[g_exons$transcriptsHits]
  names(g_exons) <- .x$name[g_exons$transcriptsHits]

  # Create a GRanges for CDS
  g_thick <- GRanges(
    seqnames(.x),
    .x$thick,
    strand(.x)
  )
  names(g_thick) <- .x$name
  
  # Intersect exons with CDS
  g_thick <- pintersect(
    g_exons,
    g_thick[g_exons$transcriptsHits],
    drop.nohit.ranges=TRUE
  )

  # Map thick exons back to transcripts
  t_thick <- GenomicFeatures::pmapToTranscripts(
    g_thick,
    transcripts[g_thick$transcriptsHits]
  ) |> ranges()
  t_thick <- split(t_thick, names(g_thick)) 
  t_thick <- sort(t_thick)
  
  # Produce the exon block with strandness
  t_exons <- GenomicFeatures::pmapToTranscripts(
    g_exons,
    transcripts[g_exons$transcriptsHits], 
    ignore.strand=FALSE)
  t_exons <- split(ranges(t_exons), seqnames(t_exons))
  t_exons <- sort(t_exons)
  
  # Mover from transcript to cDNA coordinates
  c_exons <- endoapply(t_exons, ungap)
  c_thick <- endoapply(t_thick, ungap_keep_start)
  
  .x$thick <- NULL
  .x$blocks <- NULL
  .x$score <- NULL
  .x$itemRgb <- NULL
  .x$NA. <- NULL
  .x$c_exons <- c_exons[.x$name]
  .x$t_exons <- t_exons[.x$name]
  .x$g_exons <- split(g_exons, .x$name[g_exons$transcriptsHits])[.x$name]
  .x$c_thick <- c_thick[.x$name]
  .x$t_thick <- t_thick[.x$name]
  .x$g_thick <- split(g_thick, .x$name[g_thick$transcriptsHits])[.x$name]
  return(.x)
}

input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]
if (!endsWith(output_file, '.RDS')) {
  output_file <- paste0(output_file, ".RDS")
}

if (!file.exists(input_file)) {
  cat("Input file not found: ", input_file, "\n")
  q(status = 1)
}

data <- rtracklayer::import.bed(input_file)
message('Data loaded.')
data <- map_bed12_features(data)
message('Features mapped.')
data$dist <- dist_to_stop(data)
data$ptc <- data$dist > 50
saveRDS(data, output_file)
message('Finished.')


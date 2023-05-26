library(tidyverse)
library(rtracklayer)
library(plyranges)

make_groups_unique <- function(.group) {
  paste(
    .group,
    ave(rep(NA, length(.group)), .group, FUN=seq_along),
    sep = '_' 
  )
  
}

ungap <- function(.x) {
  IRanges(end=cumsum(width(.x)), width = width(.x))
}

deseq_results <- function (dds, treat, control) { 
  message(str_glue("computing for {treat} vs {control}") )
  res <- results(dds,
                 contrast=c("group", treat, control),
                 lfcThreshold=0.5, 
                 altHypothesis="greaterAbs",
                 alpha=0.05,
                 filterFun=ihw)

  res %>%
    data.frame() %>%
    rownames_to_column(var="gene_id") %>% 
    as_tibble()
  
}

rule_50nt <- function(gr, dist_cutoff = 50) {
  if (length(gr) == 1) {
    return(FALSE)
  }

  if (all(strand(gr) == "-")) {
    gr <- gr[order(-end(gr)), ]

    block <- subset(gr, type == "exon")
    thick <- subset(gr, type == "CDS")
    # shift to the last nt of stop codon
    stop <- min(start(thick)) - 3
    last_ejc <- start(block[length(block) - 1, ])
    stop - last_ejc > dist_cutoff
  } else {
    gr <- gr[order(start(gr)), ]

    block <- subset(gr, type == "exon")
    thick <- subset(gr, type == "CDS")

    stop <- max(end(thick)) + 3
    last_ejc <- end(block[length(block) - 1, ])
    last_ejc - stop > dist_cutoff
  }
}

make_tx_list <- function(gr, feature = "exon") {
  feat <- subset(gr, type == feature)
  stopifnot(length(feat) > 0)
  return(split(feat, feat$transcript_id))
}

filter_multi_exon <- function(gtf) {
  stopifnot(is(gtf, "GRanges"))
  
  ex <- subset(gtf, type == "exon")
  stopifnot(length(ex) > 0)
  
  # discard single exons transcripts
  multi_ex <- table(ex$transcript_id) > 1
  ex <- subset(ex, mcols(ex)$transcript_id %in% names(multi_ex[multi_ex]))
  ex_tx <- split(ex, ex$transcript_id)
  ex_tx
}

make_introns <- function(gr) {
  tmp <- range(gr)
  tmp[elementNROWS(tmp) >= 2] <- heads(tmp[elementNROWS(tmp) >= 2], 1) 
  tmp <- unlist(tmp, use.names = FALSE)
  psetdiff(tmp, gr)
}

process_transdecoder <- function(transdecoder) {
  tran_type <- transdecoder$type
  levels(tran_type)[levels(tran_type) == "mRNA"] <- "transcript"
  transdecoder$type <- tran_type
  transdecoder$transcript_id <- as.character(transdecoder$Parent)
  is_tx <- transdecoder$type == "transcript"
  transdecoder[is_tx]$transcript_id <- transdecoder[is_tx]$ID
  transdecoder
}

get_seq <- function(gr, genome_seq) {
  if (as.logical(strand(gr) == "-")[[1]]) {
    gr <- gr[order(end(gr), decreasing = T)]
  } else {
    gr <- gr[order(start(gr))]
  }
  getSeq(genome_seq, gr)
}


find_stop_inframe <- function(
  .gr,
  .gr_seq,
  stop_pat = DNAStringSet(c("TAA", "TAG", "TGA"))
) {
  stop_c <- matchPDict(stop_pat, .gr_seq)
  stop_hits <- rep(as.character(stop_pat), times = lengths(stop_c))
  stop_c <- unlist(stop_c)
  mcols(stop_c)$hits <- stop_hits
  # filter for matches that are in frame with the start codon
  stop_c <- stop_c[stop_c@start %% 3 == 1]
  stop_c <- stop_c[order(stop_c@start)]
  # the end coord for the first stop codon in phase match
  end(stop_c)[1]
}

check_stop_inframe <- function(i) {
  tx_id <- runValue(seqnames(mapped_start_codon[i]))
  tx_strand <- runValue(strand(mapped_start_codon[i]))
  tx_seq <- cdna[[tx_id]]
  len_seq <- length(tx_seq)
  start_seq <- start(mapped_start_codon)[i]
  stop_seq <- find_stop_inframe(
    Tx[[tx_id]],
    tx_seq[start_seq:len_seq]
  )
  if (is.na(stop_seq)) {
    message(str_glue("{tx_id} iteration {i} no stop codon in-frame"))
    return(0)
  }
  
  trans_seq <- translate(tx_seq[start_seq:(start_seq + stop_seq - 1)])
  
  if (is.na(stop_seq)) {
    message(str_glue("{tx_id} iteration {i} no stop codon in-frame"))
    return(0)
  } else if (!safeExplode(as.character(trans_seq))[1] == "M") {
    message(str_glue("{tx_id} iteration {i} missing M at first position"))
    return(0)
  } else if (!safeExplode(as.character(trans_seq))[length(trans_seq)] == "*") {
    message(str_glue("{tx_id} iteration {i} missing * at last position"))
    return(0)
  }
  return(stop_seq)
}

#' derived from https://support.bioconductor.org/p/101245/#101776
#' cds_by_tx GRangeList of CDS per transcript
addCdsPhase <- function(cds_by_tx) {
  cds_phase <- pc(
    rep(IntegerList(0), length(cds_by_tx)),
    IRanges::heads((3L - (cumsum(width(cds_by_tx)) %% 3L)) %% 3L, n = -1L)
  )
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names = FALSE)
  # there are missing CDS that needs to be fixed.
  mcols(unlisted_cds_by_tx)$phase <- unlist(cds_phase, use.names = FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}

#' Tests whether cdna from transdecoder is the same as the one BSgenome::getSeq
test_cDNA_seq <- function(
  Tx) {
  dna <- readDNAStringSet("/prj/Niels_Gehring/smg_ko/analysis/baltica/GRCh38_90_SIRV_Set3_oneCol.fa")
  novel_seq <- lapply(Tx, BSgenome::getSeq, x=dna)
  stopifnot(
    all(
      lapply(names(novel_seq), function(i) {
        cdna[[i]] == unlist(novel_seq[[i]])
      }
      )
    )
  )
}


#' novel_cds and transcriptome should be aligned 
#' 
single_cds_to_many <- function(novel_cds, transcriptome){
  cds <- pintersect(
    novel_cds, 
    transcriptome, 
    drop.nohit.ranges=TRUE)
  # names(cds) <- seq_along(cds)
  cds <- unlist(cds)
  this_rle <- rle(names(cds))
  this_rle$values <- seq_along(this_rle$values)
  mcols(cds) <- mcols(novel_cds)[inverse.rle(this_rle), ]
  split(unname(cds), names(cds))
}

make_db <- function(file, cols = c('gene_name', 'gene_id', 'ref_gene_id')){
  .gr <- rtracklayer::import(file)
  .gr <- subset(.gr, strand != "*")
  .gr_DB <- makeTxDbFromGRanges(.gr, )
  .gr_by_tx <- exonsBy(.gr_DB, by = "tx")
  names(.gr_by_tx) <- id2name(.gr_DB, "tx")
  .gr_by_tx
}

#' based on plyranges::bind_ranges
bind_ranges <- function(.x, .id='transcript_id'){
  stopifnot(is(.x, "GRangesList"))
  .x <- unlist(.x)
  mcols(.x)[.id] <- names(.x)
  unname(.x)
}

addCdsPhase <- function (cds_by_tx) 
{
  cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)), IRanges::heads((3L - 
                                                                            (cumsum(width(cds_by_tx))%%3L))%%3L, n = -1L))
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names = FALSE)
  mcols(unlisted_cds_by_tx)$phase <- unlist(cds_phase, use.names = FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}

.grList2bed12 <- function(.gr) {
  stopifnot(inherits(.gr, "GRangesList"))
  
  len_ids <- lengths(.gr)
  names(len_ids) <- names(.gr)
  .gr <- unlist(.gr)
  .gr_type <- .gr$type
  mcols(.gr) <- NULL
  mcols(.gr)$Parent <- rep(names(len_ids), times = len_ids)
  names(.gr) <- NULL
  
  .ex <- subset(.gr, .gr_type == "exon") 
  
  .tx <- split(.ex, ~Parent) %>% range()
  
  .ex <- GenomicFeatures::pmapToTranscripts(.ex, .tx[.ex$Parent])
  .ex <- split(ranges(.ex), seqnames(.ex))

  .cds <- subset(.gr, .gr_type == "CDS")
  .cds <- GenomicFeatures::pmapToTranscripts(.cds, .tx[.cds$Parent])
  .cds <- split(ranges(.cds), seqnames(.cds))
  
  .tx <- unlist(.tx)

  mcols(.tx)$thick <- IRanges(0, 0)
  mcols(.tx)$blocks <- .ex[names(.tx)]
  mcols(.tx)$thick <- .cds[names(.tx)]
  score(.tx) <- 0
  
  mcols(.tx)$blocks <- endoapply(mcols(.tx)$blocks, sort)
  mcols(.tx)$thick <- endoapply(mcols(.tx)$thick, sort)
  
  .tx
}

grList2bed12 <- function(.gr) {
  stopifnot(inherits(.gr, "GRangesList"))
  
  len_ids <- lengths(.gr)
  names(len_ids) <- names(.gr)
  .gr <- unlist(.gr)
  mcols(.gr)$Parent <- rep(names(len_ids), times = len_ids)
  names(.gr) <- NULL
  
  .ex <- subset(.gr, type == "exon")
  .ex <- split(.ex, ~Parent)
  .tx <- range(.ex)
  .ex <- rtracklayer::asBED(.ex)
  
  .cds <- subset(.gr, type == "CDS")
  
  .cds <- pmapToTranscripts(
    split(.cds, ~Parent),
    .tx)
  .cds <- unlist(range(.cds))
  #.cds <- split(.cds, ~Parent)
  
  mcols(.ex)$thick <- IRanges(0, 0)
  if (length(.cds) == 0) {
    return(.ex)
  }
  idx <- match(names(.cds), .ex$name, nomatch = 0)
  mcols(.ex[idx, ])$thick <- range(.cds)
  score(.ex) <- 0
  
  return(.ex)
}

library(plyranges)

bed122bigBed <- function(.gr, output, chr_size="trackhub/chr.size") {

  .gr <- .gr[ start(.gr) < start(.gr$thick) & end(.gr) > end(.gr$thick)]
  
  mcols(.gr)$itemRgb <- "black"
  seqlevelsStyle(.gr) <- "UCSC"
  .gr <- rtracklayer:::sortBySeqnameAndStart(.gr)
  .gr <- keepStandardChromosomes(.gr, pruning.mode = "coarse")
  mcols(.gr)$score <- 0
  
  seq_len <- deframe(read.table("trackhub/chr.size"))
  names(seq_len)[names(seq_len) == "chrMT"] <- "chrM"
  seqlengths(.gr) <- seq_len[names(seqlengths(.gr))]
  
  .gr %>%
    select(score, name, thick, itemRgb, blocks) %>%
    export.bb(., output)
}

filter_thick <-  function(.x) { 
  .x[start(.x) < start(.x$thick) & end(.x) > end(.x$thick)]
}

# based on https://stackoverflow.com/a/68172404/1694714
make_unique <- function(.x, sep = '_') {
  vec <- as.numeric(ave(.x, .x, FUN = seq_along))
  .x <- paste(.x, vec, sep = sep)
  .x
}

process_transdecoder_score <- function() {
  cds_scores <- read_table("transdecoder/longest_orfs.cds.scores")
  
  colnames(cds_scores)[1] <- "transcript_id" 
  cds_scores <- cds_scores %>% 
    pivot_longer(
      cols = -c(transcript_id, Markov_order, seq_length), 
      names_to = 'grp',
      values_to = 'score') %>% 
    group_by(transcript_id) %>% 
    summarise(score = max(score))
  cds_scores

}

self_split <- function(.x) split(.x, seq_along(.x))

draw_bed12 <- function(.x) {
  
  blk <- .x$blocks %>% as.data.frame() 
  cds <- .x$thick %>% as.data.frame()
  
  
  ggplot(
    data = blk,
    aes(
      xstart = start,
      xend = end,
      y = group_name),
  ) + 
    geom_range( 
      fill = "white",
      height = 0.25
    ) +
    geom_range(
      data = cds,
      height = 0.50
    ) + 
    geom_intron(
      data = to_intron(blk, "group_name"),
      arrow = NULL
    ) + 
    theme_bw() +
    labs(
      x = "Rel. distance TSS (nt)",
      y = ""
    )
}

fix_seqname <- function(ref, other) {
  other <- keepStandardChromosomes(other, pruning.mode = 'coarse')
  seqlevelsStyle(other) <- seqlevelsStyle(ref)
  other
  
}

#' By the longest
#' 
annotate_gene_names <- function(ref, other) {
  hits <- findOverlaps(other, ref)
  hits <- as(hits, "Hits")
  mcols(hits)$gene_width <- width(ref[to(hits)])
  hits <- sort(hits, by = ~gene_width)
  hits <- breakTies(hits, "last")
  other$gene_name <- NA
  other$gene_name[from(hits)] <- ref[to(hits)]$gene_name
  other

}

make_ranks <- function(x){
  unlist(lapply(runLength(Rle(x)), function(.x) seq.int(1, .x)))
}

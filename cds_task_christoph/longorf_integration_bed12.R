#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicFeatures)
  library(dplyr)
  library(plyranges)
  library(stringr)
})

write_bed12 <- function(.gr, output, chr_size="trackhub/chr.size", write_bed=TRUE) {
  seqlevelsStyle(.gr) <- "UCSC"
  .gr <- rtracklayer:::sortBySeqnameAndStart(.gr)
  .gr <- keepStandardChromosomes(.gr, pruning.mode = "coarse")
  mcols(.gr)$score <- 0
  
  seq_len <- tibble::deframe(read.table("trackhub/chr.size"))
  names(seq_len)[names(seq_len) == "chrMT"] <- "chrM"
  seqlengths(.gr) <- seq_len[names(seqlengths(.gr))]

  .gr %>%
    select(score, name, thick, itemRgb, blocks) %>%
    export.bb(., output)
  
  if (write_bed) {
    system(
      paste(
        "/biosw/ucsc/2021-04-20/bigBedToBed",
        output,
        gsub('.bb', '.bed', output),
        sep = ' '))
  }
}


#' Extract Start seqs
#'
#' This function extracts the start seqs from a given set of seqs based 
#' on the positions of mature RNA.
#'
#' @param seqs A named StringSet object
#' @param mature_rna_position A named IRangesList object with CDS positions.
#' 
#' @return Returns a subset of the input seqs corresponding to the start positions
#'         defined by the `mature_rna_position`.
#' 
#' table(unlist(starts))
#' 
#' @examples
#' # Assuming seqs and mature_rna_position are properly defined:
#' start_seqs <- extract_start(seqs, mature_rna_position)
#' table(unlist(start_seqs))
#'   
#' @export
extract_start <- function(seqs, mature_rna_position) {
  seqs <- seqs[names(mature_rna_position)]
  mature_rna_position <- heads(mature_rna_position, 1) |> unlist() |> resize(3)
  extractAt(seqs, split(mature_rna_position))  
}

 
flag_ptc <- function(.x) {
  leej <- .x$cdna_blocks |>
    tails(1) |>
    start() - 3 |>
    unlist()
  stop_pos <- .x$cdna_thick |>
    tails(1) |>
    end()
  empty_stop <- lapply(stop_pos, isEmpty) |> unlist()
  stop_pos[empty_stop] <- leej[empty_stop]
  stop_pos <- stop_pos |> unlist()
  # +3 because the CDS does not include the stop codon
  unlist(leej) - stop_pos
}
# -----------------------
# 1. DATA LOADING
# -----------------------
ref_match <- read.table(
  "/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/compared/fix_comp_ref.tracking",
  col.names = c("query_id", "locus_id", "ref", "class_code", "sample_info")
)
ref_match$transcript_id <- str_split(ref_match$sample_info, "\\|", simplify = TRUE)[, 2]
ref_match <- ref_match %>%
  mutate(class_code = case_when(
    class_code == '=' ~ 'same_intron_chain',
    class_code == 'n' ~ 'IR',
    class_code %in% c('c', 'j', 'k') ~ 'splicing_variants',
    TRUE ~ 'other'
  ))
ref_match <- tibble::deframe(ref_match %>% select(transcript_id, class_code))

# 2. Import ref, ensembl_102 and ref_bed12
ref <- import("/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/stringtie_merge/merged_each.fix.gtf")
ensembl_102 <- rtracklayer::import("/beegfs/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")
ref_bed12 <- rtracklayer::import.bed("/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/cds_task_christoph/ref.bed12")

# 3. Load and process cdna sequence
fa <- Biostrings::readAAStringSet("/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/fasta/transcriptome.clean.fa")
names(fa) <- gsub(" .*$", "", names(fa))
cds <- setNames(nm = ref_bed12$name, GRanges(seqnames(ref_bed12), ref_bed12$thick, strand = strand(ref_bed12)))

# -----------------------
# 2. DATA WRANGLING
# -----------------------

txdb <- ref |> subset(type == 'exon') |> dplyr::select(transcript_id) |> split(~transcript_id)

add_cdna_blocks <- function(tx) {
  # # Adjust blocks for negative strand
  neg_idx <- which(strand(tx) == "-")
  tx$blocks <- revElements(tx$blocks, neg_idx)
  tx$cdna_blocks <- endoapply(
    tx$blocks,
    function(.x) IRanges(end = cumsum(width(.x)), width = width(.x))
  )
  
  tx
}

process_data_from_ensembl <- function(txdb, ref_match, cds) {
  # Convert to bed12 format
  tx <- asBED(txdb)
  
  tx$gffcompare <- ref_match[tx$name]
  
  # Filter common elements
  common <- intersect(names(cds), names(txdb))
  common <- common[common %in% names(ref_match[ref_match == 'same_intron_chain'])]
  cds <- cds[common]
  
  # Drop cds that don't overlap with transcripts
  cds <- pintersect(cds, txdb[names(cds)], drop.nohit.ranges=TRUE)
  cds <- range(cds) |> unlist()
  
  # Match Ensembl
  cds <- suppressWarnings(c(cds, GRanges(seqnames='EMPTY',ranges=IRanges(1, 0))))
  matches <- match(tx$name, names(cds), nomatch = length(cds))
  tx$thick <- ranges(cds[matches])
  tx$source <- ""
  tx$source[width(tx$thick) > 1] <- 'canonical'
  
  tx
}

compute_ensembl_canonical <- function(
    tx, 
    gffread_path = "/biosw/gffread/0.12.6/gffread", 
    bed_filename = "ensembl_canonical.bed", 
    gtf_filename = "ensembl_canonical.gtf") {
  # Filter tx based on thick width and export to BED format
  filtered_tx <- tx[width(tx$thick) > 1, ]
  rtracklayer::export(filtered_tx |>  select(name, blocks, thick), bed_filename)
  
  # Convert BED to GTF using gffread
  cmd <- sprintf(
    "%s -T --in-bed %s -o %s",
    gffread_path, 
    bed_filename,
    gtf_filename)
  
  try(system(cmd))
  
  # Import the GTF file and process
  imported_tx <- rtracklayer::import(gtf_filename)
  exon_groups <- imported_tx |> 
    subset(type == 'exon') |> 
    split(~transcript_id)
  bed12_format <- asBED(exon_groups)
  
  cds_data <- imported_tx |> 
    subset(type == 'CDS') |> 
    split(~transcript_id)
  
  cds_transcripts2 <- pmapToTranscripts(
    unlist(cds_data), 
    exon_groups[names(unlist(cds_data))])
  
  # Process and shift the ranges
  split_cds_ranges <- split(ranges(cds_transcripts2), names(cds_transcripts2))
  split_cds_ranges <- sort(split_cds_ranges)
  stopifnot(isEmpty(gaps(split_cds_ranges)))
  
  bed12_format$thick <- range(ranges(cds_data[bed12_format$name]))
  bed12_format$cdna_thick <- split_cds_ranges[bed12_format$name]
  
  return(bed12_format)
}

tx <- process_data_from_ensembl(txdb, ref_match, cds)
bed12_format <- compute_ensembl_canonical(tx)
# starts <- extract_start(fa[bed12_format$name], range(bed12_format$cdna_thick))
# sort(-prop.table(table(unlist(starts))))
bed12_format <- add_cdna_blocks(bed12_format) 
bed12_format$is_ptc <- flag_ptc(bed12_format) > 50
bed12_format$itemRgb <- ifelse(bed12_format$is_ptc, 'red', 'black')

bed12_format <- rtracklayer:::sortBySeqnameAndStart(bed12_format)
bed12_format$source <- factor('canonical', levels = c("canonical", "ensembl", "riboseq",  "openprot"))
bed12_format$thick <- unlist(bed12_format$thick)
# removes stop codon from CDS canonical sources
bed12_format$thick <- resize(bed12_format$thick, width = width(bed12_format$thick) - 2)

load_orfs_data <- function(input_file = "longorf2_output.bed") {
  # 1. Import the bed file
  orfs_cdna <- rtracklayer::import.bed(input_file)
  
  # 2. Map from transcripts
  orfs_geno <- mapFromTranscripts(orfs_cdna, txdb)
  
  # 3. Extract and assign source, transcript ID, CDS ID and type
  orfs_geno$source <- str_extract(orfs_cdna$name[orfs_geno$xHits], 'ensembl|openprot|riboseq')
  orfs_geno$transcript_id <- seqnames(orfs_cdna)[orfs_geno$xHits]
  orfs_geno$cds_id <- make.unique(as.character(orfs_geno$transcript_id), '_')
  orfs_geno$type <- 'CDS'
  
  # 4. Convert txdb to BED and update its attributes
  bed12 <- asBED(txdb)
  bed12$thick <- NA
  bed12 <- bed12[orfs_geno$transcriptsHits]
  bed12$thick <- orfs_geno
  mcols(bed12$thick) <- NULL
  bed12$source <- orfs_geno$source
  bed12$cds_id <- orfs_geno$cds_id
  
  # 5. Adjust start positions
  start(bed12) <- start(bed12) + 1
  start(bed12$thick) <- start(bed12$thick) + 1
  
  return(bed12)
}


## Read longorf2 output 
longorf2 <- import('phaseFinal/cds_task_christoph/longorf2_output.bed')
read_and_process_seqs <- function(filepath, source) {
  # Read the DNA or AA StringSet depending on the file extension
  if (grepl("ensembl", filepath)) {
    seqs <- Biostrings::readDNAStringSet(filepath)
    seqs <- Biostrings::translate(seqs)
  } else {
    seqs <- Biostrings::readAAStringSet(filepath)
  }
  
  
  # Clean up names
  input_string <- names(seqs)
  values <- str_extract_all(input_string, "[A-Z0-9.]+|\\d+|\\+|ATG")
  mcols(seqs) <- do.call(rbind, values)
  colnames(mcols(seqs)) <- c("id", "length", "start", "end", "orientation", "ATG", "ATG_position", "tx_length")
  mcols(seqs) <- mcols(seqs)[c("id", "start", "end", "ATG_position")]
  names(seqs) <- mcols(seqs)$ENST_ID
  
  # Set metadata columns
  mcols(seqs)$source <- source
  mcols(seqs)$n <- 1: length(seqs)
  
  return(seqs)
}

a <- read_and_process_seqs("phaseFinal/cds_task_christoph/start_codons_ensembl_longorf2", 'ensembl')
b <- read_and_process_seqs("phaseFinal/cds_task_christoph/start_codons_riboseq_orfs_longorf2", 'riboseq')
c <- read_and_process_seqs("phaseFinal/cds_task_christoph/start_codons_human-openprot_longorf2", 'openprot')

longorf_seqs <- c(a, b, c)

longorf <- mcols(c(a, b, c))
longorf$width <- width(c(a, b, c))
longorf <- longorf %>% 
  as_tibble() %>% 
  mutate(
    source = factor(source, levels = c('ensembl', 'riboseq', 'openprot')),
    seq = as.character(longorf_seqs)) 

longorf_unique <- longorf %>%
  group_by(id, source) %>% 
  arrange(source, desc(width)) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  distinct(seq, .keep_all = TRUE) %>% 
  ungroup()

longorf_final <- asBED(txdb[longorf_unique$id])
longorf_final$source <- longorf_unique$source
longorf_final$cdna_thick <- IRanges(as.integer(longorf_unique$ATG_position), as.integer(longorf_unique$end))
longorf_final <- add_cdna_blocks(longorf_final)
# thick <- pmapFromTranscripts(
#   GRanges(longorf_final$name, longorf_final$cdna_thick), txdb[longorf_final$name])

thick <- pmapFromTranscripts(
  GRanges(longorf_final$name, longorf_final$cdna_thick), txdb[longorf_final$name])
thick_id <- rep(seq(thick), lengths(thick))
thick <- unlist(thick)
thick$id <- thick_id
thick <- subset(thick, hit == TRUE)
thick <- split(thick, ~id)
thick_tx <- heads(thick, 1) %>% unlist() %>% mcols()
thick <- range(thick) 
thick <- unlist(thick)
mcols(thick) <- thick_tx

longorf_final$thick <- thick
longorf_final$is_ptc <- flag_ptc(longorf_final) > 50
longorf_final$itemRgb <- ifelse(longorf_final$is_ptc, 'red', 'black')
longorf_final$thick <- ranges(longorf_final$thick)
longorf_final <- longorf_final[
  start(longorf_final) < start(longorf_final$thick) &
    end(longorf_final) > end(longorf_final$thick)]

write_bed12(
  bed12_format, 'phaseFinal/data/canonical.bb')
write_bed12(
  subset(longorf_final, source == 'ensembl'), 'phaseFinal/data/ensembl.bb')
write_bed12(
  subset(longorf_final, source == 'riboseq'), 'phaseFinal/data/riboseq.bb')
write_bed12(
  subset(longorf_final, source == 'openprot'), 'phaseFinal/data/openprot.bb')

final <- c(bed12_format, longorf_final)
final$blocks <- NULL
final$cdna_blocks <- lapply(unname(final$cdna_blocks), as.data.frame)
final$cdna_thick <- range(final$cdna_thick) %>% as.data.frame() %>% select(-c(group, group_name))
final$thick <- as.character(final$thick) %>% unname()
final <- rtracklayer:::sortBySeqnameAndStart(final)
final <- as.data.frame(final)
saveRDS(final, 'phaseFinal/longorf_bed12.RDS')
load('longorf_integration_bed12.Rdata')

dist <- data.frame(
  name = longorf_final$name,
  dist=flag_ptc(longorf_final)
)

suppressPackageStartupMessages({
  library(dplyr)
  library(plyranges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(here)
})

snakemake@source('utils.R')

ref <- readRDS(snakemake@input[[1]])
cds_seq <- readRDS(snakemake@input[[2]])
cds <- readRDS(snakemake@input[[3]])

cds$source <- factor(cds$source, levels = c(
  "ensembl", "hek293gao", "ribotish", "openprot", "transdecoder"))
cds <- sort(cds, by=~source, decreasing=FALSE, ignore.strand=FALSE)


gr <- ref %>%  
  filter(type == 'exon') %>% 
  plyranges::select(transcript_id) %>% 
  filter(transcript_id %in% cds_seq$transcript)

grl <- split(gr, ~transcript_id) 
grl2 <- grl[lengths(runValue(strand(grl))) == 1 & lengths(runValue(seqnames(grl))) == 1]
stopifnot(length(grl) == length(grl2))

gr <- unlist2(grl)
names(gr) <- NULL
stopifnot(sum(lengths(grl)) == length(gr))
stopifnot(isEmpty(setdiff(cds_seq$transcript, names(grl))))

cds_seq$granges <- pmapFromTranscripts(cds_seq$ranges, grl[cds_seq$transcript] )

tr <- pmapToTranscripts(gr, grl[gr$transcript_id])
trl <- split(tr, seqnames(tr))

all_cds <- readRDS("phase2/data/cds_seq.RDS")
all_cds_t <- GRanges(seqnames = all_cds$transcript, all_cds$ranges, use.names = FALSE)
all_cds_g <- pmapFromTranscripts(all_cds_t, grl[seqnames(all_cds_t)] )

meta <- mcols(all_cds_g)
mcols(all_cds_g) <- NULL

hits <- findOverlaps(
  all_cds_g %>% unlist2() %>% anchor_5p() %>% mutate(width = 3),
  cds %>% anchor_5p() %>% mutate(width = 3),
  select = 'first'
)

res <- DataFrame(
  cds = all_cds_g[!is.na(hits)],
  source = cds[ hits[!is.na(hits)], ]$source,
  transcript = names(grl)[meta$transcriptsHits[!is.na(hits)]],
  cds_id = all_cds_names[ meta$xHits[!is.na(hits)]],
  cds_tx_pos = all_cds[meta$xHits[!is.na(hits)]]
)

grl <- grl[res$transcript]
grl2 <- unlist(grl)
grl2 <- pmapToTranscripts(grl2, grl[names(grl2)])
res$grl <- unname(grl2[res$transcript])
res$stop <- res$cds %>% end() %>% unlist2()
res$last_ejc <- res$grl %>% tails(2) %>% heads(1) %>% end() %>% unlist2()
res$ptc <- res$last_ejc - 1 - res$stop > 50 

table(res$source, res$ptc) %>% 
  kbl() %>% 
  kable_styling()

saveRDS(res, here('data', '02_task_simple.RDS'))



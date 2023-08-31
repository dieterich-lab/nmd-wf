library(rtracklayer) 
library(plyranges)

ref_gtf <- rtracklayer::import("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")
ref <- import.bed('ref.bed12')
biotype <- filter(ref_gtf, type == "transcript") %>% as.data.frame() %>% plyranges::select(transcript_id, transcript_biotype) %>% deframe()
mcols(ref)$itemRgb <- ifelse(biotype[ref$name] == 'nonsense_mediated_decay', 'red', 'black') 
seqlevelsStyle(ref) <- "UCSC"
ref <- rtracklayer:::sortBySeqnameAndStart(ref)
ref <- keepStandardChromosomes(ref, pruning.mode = "coarse")
mcols(ref)$score <- 0

seq_len <- deframe(read.table("trackhub/chr.size"))
names(seq_len)[names(seq_len) == "chrMT"] <- "chrM"
seqlengths(ref) <- seq_len[names(seqlengths(ref))]

ref %>% 
  plyranges::select(name, score, itemRgb, thick, blocks) %>% 
  export.bb(., "/prj/trackhubs/nmd_transcriptome/hg38/ref.bigBed")

new_gtf <- rtracklayer::import("phaseFinal/data/gtf_annotated.gtf")
new <- import.bed('new.bed12')
color <- filter(new_gtf, type == "transcript") %>% as.data.frame() %>% plyranges::select(transcript_id, color) %>% deframe()
mcols(new)$itemRgb <- ifelse(color[new$name] == '#FF0000', 'red', 'black') 

seqlevelsStyle(new) <- "UCSC"
new <- rtracklayer:::sortBySeqnameAndStart(new)
new <- keepStandardChromosomes(new, pruning.mode = "coarse")
mcols(new)$score <- 0
seqlengths(new) <- seq_len[names(seqlengths(new))]

new %>% 
  plyranges::select(name, score, itemRgb, thick, blocks) %>% 
  export.bb(., "/prj/trackhubs/nmd_transcriptome/hg38/new.bigBed")

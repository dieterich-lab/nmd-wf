suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(plyranges)
  library(ggplot2)
  library(here)
})

snakemake@source('utils.R')

ref <- import(snakemake@input[[2]])
genes <- ref %>% filter(type == 'gene')

# from ensembl
ensembl <- ref %>% filter(type == 'CDS') %>% split(~transcript_id) %>% range() %>% unlist()
ensembl$source <- 'ensembl'

# from Christoph
hek293gao <- import("cds_source/hek293-gao-unique.filtered.predicted-orfs.gtf")
hek293gao <- hek293gao %>% filter(type == 'CDS') %>% split(~transcript_id) %>% range() %>% unlist()
hek293gao$source <- 'hek293gao'
names(hek293gao) <- gsub(x=names(hek293gao), '(.*)_.*', '\\1') 
hek293gao <- fix_seqname(ref, hek293gao)
hek293gao <- annotate_gene_names(genes, hek293gao)

# from openprot
openprot <- import("cds_source/human-openprot-r1_6-refprots+altprots+isoforms_min_2_pep-.bed")
openprot <- split(openprot, ~name) %>% range() %>% unlist()
openprot$source <- 'openprot'
openprot <- fix_seqname(ref, openprot)
openprot <- annotate_gene_names(genes, openprot)

transdecoder <- import("transdecoder/longest_orfs.cds.best_candidates.genome.gff3")
transdecoder_score <- gsub(transdecoder$Name, pattern = ".*,score=([-0-9.]+)$", replacement = "\\1") %>% 
  as.numeric() %>% 
  zoo::na.locf()
transdecoder <- transdecoder[ transdecoder_score > 10, ]

transdecoder_complete <- transdecoder$Name %>% 
  zoo::na.locf() %>% 
  grepl(x=., pattern = ':complete')
transdecoder <- transdecoder[transdecoder_complete, ]

transdecoder <- subset(transdecoder, type == 'CDS')

transdecoder$Parent <- transdecoder$Parent %>% unlist()
transdecoder <- split(transdecoder, ~Parent) %>% range() %>% unlist()
transdecoder$source <- 'transdecoder'
transdecoder <- fix_seqname(ref, transdecoder)
transdecoder <- annotate_gene_names(genes, transdecoder)

# from https://github.com/boehmv/NMD-Transcriptome/files/4270509/pred.tar.gz
ribotish <- read.table('cds_source/zhangribotish.txt', header = 1)
ribotish <- GRanges(ribotish$GenomePos, name = ribotish$Tid)
names(ribotish) <- ribotish$name
ribotish$name <- NULL
ribotish$source <- 'ribotish'
ribotish <- fix_seqname(ref, ribotish)
ribotish <- annotate_gene_names(genes, ribotish)

gr <- c(ensembl, hek293gao, openprot, transdecoder, ribotish)
gr$source <- factor(gr$source, levels = c(
  "ensembl", "hek293gao", "ribotish", "openprot", "transdecoder"))
gr <- sort(gr, by=~source, decreasing=FALSE, ignore.strand=FALSE)
message('Table with redudant hits')
before <- table(gr$source)

hits <- findOverlaps(
  gr %>% anchor_5p() %>% mutate(width = 1), 
  drop.redundant=TRUE, 
  drop.self=TRUE)

message('Table after removing redudant hits')
gr <- gr[-to(hits)]  
table(gr$source)
gr <- sortSeqlevels(gr)
gr <- sort(gr)
after <- table(gr$source)

saveRDS(gr, snakemake@output[[1]])

data <- tibble(
  redudant = as.numeric(before),
  unique = as.numeric(after),
  names = names(after)) %>% 
  pivot_longer(cols = -names)


data %>% 
  mutate(names = factor(names, levels= c("ensembl", "hek293gao", "ribotish", "openprot", "transdecoder"))) %>% 
  ggplot(aes(fill=name, y=value, x=names)) + 
  scale_fill_brewer(palette = "Spectral", name="") +
  geom_bar(position="dodge", stat="identity") +
  geom_text(
    aes(label = scales::comma(value)),
    size = 3,
    position=position_dodge(width=0.9), vjust=-0.25
  ) +
  theme_light() +
  labs(
    x= "", 
    y="", 
    title='Number of CDS features', 
    subtitle = 'Redundat means shared 5p') +
  theme(legend.position='top') +
  guides(fill=guide_legend(ncol=2))

ggsave(snakemake@output[[2]], height = 5, width =5*1.4, dpi = 300)

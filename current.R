library(Biostrings)
library(rtracklayer)
library(plyranges)
library(stringr)

setwd('/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal')
fasta <- 'fasta/transcriptome.fa'
cdna_seq <- Biostrings::readDNAStringSet("fasta/transcriptome.fa")
ref <- import("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")
ref_grl <- ref %>% 
  filter(type == 'exon') %>% 
  select(transcript_id) %>%
  split(., ~transcript_id) 
ref12 <- ref_grl%>% asBED()
cds <- ref %>% filter(type == 'CDS') %>% split(., ~transcript_id) %>% range()
cds_c <- GenomicFeatures::mapToTranscripts(unlist(cds), ref_grl[names(cds)])

ref_match <- read.table(
  "compared/fix_comp_ref.tracking",
  col.names = c("query_id", "locus_id", "ref", "class_code", "sample_info"))
ref_match$transcript_id <- str_split(ref_match$sample_info, "\\|", simplify = TRUE)
ref_match$transcript_id  <- ref_match$transcript_id[, 2]
ref_match <- ref_match %>% 
  mutate(class_code = case_when(
    class_code == '=' ~ 'same_intron_chain',
    class_code == 'n' ~ 'IR',
    class_code %in% c('c', 'j', 'k') ~ 'splicing_variants',
    TRUE ~ 'other'
))

gtf_data <- import("stringtie_merge/merged_each.fix.gtf")
# Filter for 'exon' type and select transcript IDs
exon_list <- gtf_data %>%
  filter(type == 'exon') %>%
  select(transcript_id) %>%
  split(., ~transcript_id)
# Convert to BED format
bed_data <- asBED(exon_list)
names(bed_data) <- bed_data$name
# add gene name class code etc

bed_data$thick <- IRanges(start = 1, width = 0)
bed_data[names(cds_c)]$thick <- cds_c

# ensembl
# if
# riboseq
# openprot
# ab_initio_longest_orf
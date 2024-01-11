## DTE kernel ###
##
source('phaseFinal/scripts/plot_utils.R')

# abstract
db$dte %>% 
  filter(padj < 0.05, log2fold > 1) %>%
  group_by(transcript_id) %>%
  tally() %>% 
  filter(n >= 3) %>% 
  nrow()


# De novo transcriptome assembly

db$anno$transcript_id %>% n_distinct()

db$anno$gene_id %>% n_distinct()

db$anno %>% 
  filter(lr_support == TRUE) %>% 
  pull(transcript_id) %>% 
  n_distinct()

# Identification of coding sequences and PTC status
db$bed12 %>% filter(source == 'canonical') %>% select(name) %>% n_distinct()

db$bed12 %>% filter(source != 'canonical') %>% select(name) %>% n_distinct()

unnested_source <- db[["anno"]] %>% unnest(source)

unnested_source %>% 
  filter(
    match == "same_intron_chain",
    source == "canonical") %>%
  nrow()


ptc_per_source = left_join(
  db[["bed12"]], db$anno %>% select(transcript_id, match),
  by = c("name"="transcript_id")) %>%
  group_by(match, source, is_ptc) %>%
  summarise(n = n())

ptc <- ptc_per_source  %>% 
  filter(is_ptc == TRUE)

sum(ptc$n)


db$dge %>%  
   filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
   group_by(contrasts) %>%
   summarise(n = n())

db$dte %>% 
  filter(padj < 0.05, abs(log2fold) > 1) %>%
  group_by(contrasts) %>%
  summarise(n = n())



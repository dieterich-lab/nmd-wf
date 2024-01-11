source('plot_utils.R')

library(circlize)
library(ComplexHeatmap)

lt <- db$anno %>% 
  select(transcript_id, source) %>%
  tidyr::unnest(source) %>% 
  mutate(value = TRUE) %>% 
  tidyr::pivot_wider(
    names_from = source,
    values_from = value, 
    values_fill = list(value = FALSE))
m1 = make_comb_mat(lt)

tx_anno <- db$anno %>% 
  select(transcript_id, match, lr_support) %>%
  group_by(transcript_id) %>% 
  summarise(is_novel = match != 'same_intron_chain', 
            lr_support = any(lr_support)) %>% 
  filter(transcript_id %in% lt$transcript_id) %>% 
  dplyr::arrange(match(transcript_id, lt$transcript_id))

comb_elements = lapply(comb_name(m1), function(nm) extract_comb(m1, nm))

col_fun = colorRamp2(c(0, .8), c("white", "red"))


us <- UpSet(
  m1,
  top_annotation =  upset_top_annotation(m1, add_numbers = T),
  bottom_annotation = HeatmapAnnotation(
    prop_novel = sapply(comb_elements, function(ind) mean(tx_anno$is_novel[ind])),
    prop_lr_support = sapply(comb_elements, function(ind) mean(tx_anno$lr_support[ind])),
    annotation_name_side = "left",
    col = list(prop_novel = col_fun, prop_lr_support = col_fun)
  )
)

png("fig11.png")
draw(us)
dev.off()

pdf(file="fig11.pdf",
    width = 7.87402,
    height = 5.51181,
    useDingbats = FALSE
)
draw(us)
dev.off()

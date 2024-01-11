source('phaseFinal/scripts/plot_utils.R')


unnested_source <- db[["anno"]] %>% unnest(source)
metadata <- load_metadata(db)
dte <- left_join(db$dte, db$bed12, c('transcript_id'='name'))
dte <- dte %>% left_join(metadata %>% select(contrasts, name))
dte <- dte %>% mutate(name = str_replace_all(name, "_", " "))

dte_mean <- dte %>%
  filter(!is.na(is_ptc), padj < 0.05) %>%
  select(is_ptc, log2fold, name) %>%
  group_by(is_ptc, name) %>%
  mutate(mean_l2fc = median(log2fold, na.rm = TRUE))

dte %>%
  filter(!is.na(is_ptc), padj < 0.05) %>%
  ggplot(aes(x = log2fold, fill=is_ptc, alpha=.5)) +
  scale_fill_manual(
    values = c("black", "firebrick", "gray50"),
    name = "is_PTC"
  ) +
  geom_density() +
  geom_vline(data = dte_mean, aes(xintercept=mean_l2fc, color=is_ptc), linetype = 2) +
  scale_color_manual(
    values = c("black", "firebrick", "gray50"),
    name = "is_PTC"
  ) +
  labs(x = "Transcript L2FC", y = "Density") +
  theme_linedraw() +
  theme_Publication() +
  facet_wrap(~name) +
  coord_cartesian(xlim=c(-5, 5)) +
  guides(alpha = FALSE)


dte %>%
  filter(padj < 0.05, abs(log2fold) > 1) %>%
  group_by(contrasts) %>%
  summarise(n = n())

ggsave(
  "fig5.pdf", 
  width = 10, 
  height = 7, 
  dpi = 300)

ggsave(
  "fig5.png", 
  width = 10, 
  height = 7, 
  dpi = 300)

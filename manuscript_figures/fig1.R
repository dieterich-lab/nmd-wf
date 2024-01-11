source('phaseFinal/scripts/plot_utils.R')

x <- left_join(db[["anno"]], db[["bed12"]], by = c("transcript_id"="name")) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  count(match, lr_support) %>% 
  mutate(prop = n / sum(n))

x %>%
  mutate(
    match = forcats::fct_reorder(match, desc(n)),
    label1 = scales::comma_format()(n),
    label2 = scales::label_percent()(prop)
    ) %>%
  ggplot(aes(
    fill = lr_support,
    y = n,
    x = match)
  ) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label=label1), size = 5, position = position_dodge(0.9), color = "black", vjust = -0.1) +
  geom_text(aes(label=label2), size = 4, position = position_dodge(0.9), color = "white", vjust = 1.2) +
  scale_fill_manual(
    values = c("black", "firebrick", "gray50"),
    name = "Has LR support"
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    n.breaks = 5L,
    expand = c(0, 0, 0.15, 0)
  ) +
  theme_linedraw() +
  theme_Publication() +
  labs(
    x = "Reference match",
    y = "Numb. of transcripts"
  )

ggsave(
  "fig1.pdf", 
  width = 10, 
  height = 7, 
  dpi = 300)

ggsave(
  "fig1.png", 
  width = 10, 
  height = 7, 
  dpi = 300)

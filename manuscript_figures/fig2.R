source('plot_utils.R')

x <- left_join(
  db[["bed12"]], db$anno %>% select(transcript_id, match),
  by = c("name"="transcript_id")) %>%
  group_by(match, source, is_ptc) %>%
  summarise(n = n())

x %>%
  mutate(match = forcats::fct_reorder(match, n)) %>%
  ggplot(aes(
    fill = is_ptc,
    y = n,
    x = match,
    label = scales::comma_format()(n))
  ) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(size = 3, position = position_dodge(0.9), color = "black", vjust = -0.2) +
  scale_fill_manual(
    values = c("black", "firebrick", "gray50"),
    labels = c("FALSE", "TRUE", "NO CDS"),
    name = "is_PTC"
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    n.breaks = 5L,
    expand = c(0, 0, 0.20, 0)
  ) +
  theme_linedraw() +
  theme_Publication() +
  facet_wrap( ~ source, scales = "free_x") +
  labs(
    x = "",
    y = "Numb. of transcripts"
  )

ggsave(
  "fig2.pdf", 
  width = 10, 
  height = 7, 
  dpi = 300)

ggsave(
  "fig2.png", 
  width = 10, 
  height = 7, 
  dpi = 300)

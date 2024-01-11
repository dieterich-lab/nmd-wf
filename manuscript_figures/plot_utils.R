suppressPackageStartupMessages({
  library(tidyverse)
  library(grid)
  library(ggthemes)
  library(ggplot2)
})

db <- readRDS('database.RDS')

theme_Publication <- function(base_size=14, base_family = "Helvetica") {

  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.spacing = unit(1, "lines"),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1.5)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

load_metadata <- function(db) {
  db[["metadata"]] %>%
    select(contrasts, Knockdown, Knockout, cellline, clone) %>%
    filter(Knockdown != "LucKD") %>%
    distinct() %>%
    mutate(
      name = str_replace(contrasts, "-vs-.*", ""),
      name = str_replace(name, "KD-KD", "-KD"),
      Knockout = str_replace(Knockout, "_", ""),
      contrast_label = str_glue_data(
        .,
        "<b> {Knockdown} </b> <br> <i>Cell-line</i>: {cellline}; <i>KO</i>: {Knockout} <i>clone</i>: {clone}"
      )
    )
}


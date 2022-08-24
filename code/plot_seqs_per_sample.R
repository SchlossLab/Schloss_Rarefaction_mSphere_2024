#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(magrittr)

set.seed(19760620)

datasets <- c("bioethanol", "human", "lake", "marine", "mice", "peromyscus",
              "rainforest", "rice", "seagrass", "sediment", "soil", "stream")

pretty_datasets <- c("Bioethanol", "Human", "Lake", "Marine", "Mice", "Peromyscus",
              "Rainforest", "Rice", "Seagrass", "Sediment", "Soil", "Stream")

group_count <- glue("data/{datasets}/data.group_count") %>%
    set_names(datasets) %>%
    map_dfr(.x =. , ~read_tsv(.x, col_names = c("group", "n_seqs"),
                              col_types = cols(group = col_character(),
                                               n_seqs = col_double())),
            .id = "dataset")


remove_groups <- glue("data/{datasets}/data.remove_accnos") %>%
    set_names(datasets) %>%
    map_dfr(.x =. , ~read_tsv(.x, col_names = c("group"),
                              col_types = cols(group = col_character())),
            .id = "dataset") %>%
    mutate(state = "remove")

jitter_plot <- left_join(group_count, remove_groups, by = c("dataset", "group")) %>%
    mutate(state = if_else(is.na(state), "keep", "remove"),
           dataset = factor(dataset, levels = rev(datasets),
                            labels = rev(pretty_datasets))) %>%
    ggplot(aes(x = n_seqs, y = dataset, color = state)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.3) +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 1e4, 1e5, 1e6),
                  labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000"),
                  limits = c(0.8, 1e6)) +
    scale_color_manual(name = NULL,
                       breaks = c("keep", "remove"),
                       labels = c("Kept", "Removed"),
                       values = c("black", "red")) +
    labs(y = NULL,
         x = "Number of sequences per sample") +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color= "gray"),
          legend.key = element_blank(),
          axis.line = element_line())

ggsave("figures/seqs_per_sample.tiff", jitter_plot, width = 7, height =4,
       compression = "lzw")

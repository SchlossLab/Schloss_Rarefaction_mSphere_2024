#!/usr/bin/env Rscript

library(tidyverse)
library(tidyverse)
library(glue)
library(patchwork)

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}


datasets <- c("bioethanol", "human", "lake", "marine", "mice", "peromyscus",
              "rainforest", "rice", "seagrass", "sediment", "soil", "stream")

pretty_datasets <- tibble(
  plain = datasets,
  pretty = capitalize(plain)
)

beta_metrics <- c(
  "rare_jaccard",
  "raw_jaccard",
  "relabund_jaccard",
  "srs_jaccard",
  "metagenomeseq_jaccard",

  "rare_bray",
  "raw_bray",
  "relabund_bray",
  "srs_bray",
  "metagenomeseq_bray",

  "rare_euclidean",
  "raw_euclidean",
  "relabund_euclidean",
  "rclr_euclidean",
  "oclr_euclidean",
  "nclr_euclidean",
  "zclr_euclidean",
  "deseq2_euclidean"
)

beta_labels <- c(
  "Rarefied",
  "Raw",
  "Rel. abundance",
  "SRS Normalized",
  "CSS Normalized",

  "Rarefied",
  "Raw",
  "Rel. abundance",
  "SRS Normalized",
  "CSS Normalized",

  "Rarefied",
  "Raw",
  "Rel. abundance",
  "Robust CLR",
  "One CLR",
  "Nudge CLR",
  "Zero CLR",
  "DeSeq2"
)

pretty_beta_calcs <- c(
  jaccard = "Jaccard",
  bray = "Bray-Curtis",
  euclidean = "Euclidean"
)

pretty_simulations <- c(
  null = "Unbiased null model",
  size = "Biased null model",
  effect = "Skewed abundance model"
)

cluster_summary_files <- glue("data/{datasets}/cluster_analysis.tsv")
names(cluster_summary_files) <- datasets

beta_composite <- map_dfr(cluster_summary_files, read_tsv, .id = "dataset") %>%
  mutate(process_calc = glue("{corrections}_{distances}")) %>%
  inner_join(., tibble(process_calc = beta_metrics),
             by = "process_calc") %>%
  mutate(process_calc = factor(process_calc, levels = beta_metrics),
         distances = factor(distances, levels = names(pretty_beta_calcs)),
         simulation = case_when(model == "r" & design == "r" ~ "null",
                                model == "r" & design == "s" ~ "size",
                                model == "e" & design == "e" ~ "effect"),
         simulation = factor(simulation, levels = c("null", "size", "effect"))
         ) %>%
  drop_na(accuracies) %>%
  group_by(simulation, dataset, corrections, distances, process_calc) %>%
  summarize(mean = mean(accuracies), sd = sd(accuracies), .groups = "drop")

null <- beta_composite %>%
  ggplot(aes(x = process_calc, y = mean,
             color = dataset, shape = dataset)) +
    geom_hline(yintercept = 0.48, size = 0.75, color = "black") +
    # geom_linerange(aes(ymin = mean - sd, ymax = mean + sd),
    #                position = position_dodge(width = 0.3)) +
    geom_point(position = position_dodge(width = 0.3)) +
    facet_grid(simulation ~ distances,
               scales = "free_x", space = "free_x",
              labeller = labeller(distances = pretty_beta_calcs,
                                  simulation = pretty_simulations)
                ) +
    scale_color_manual(name = NULL,
                       values = rep(c('#1b9e77','#d95f02','#7570b3'), 4),
                       breaks = pretty_datasets$plain,
                       labels = pretty_datasets$pretty
                       ) +
    scale_shape_manual(name = NULL,
                       values = rep(c(15, 17, 19, 25), 3),
                       breaks = pretty_datasets$plain,
                       labels = pretty_datasets$pretty
                       ) +
    scale_x_discrete(breaks = beta_metrics, labels = beta_labels) +
    scale_y_continuous(limits = c(0.48, 1.0), breaks = seq(0.5, 0.9, 0.2),
                       expand = c(0,0)) +
    labs(x = NULL, y = "Accuracy") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(12, "pt"),
      legend.position = "bottom",
     # axis.line.x = element_line()
      )

ggsave("figures/cluster_analysis.tiff",
       null,
       width = 6, height = 7,
       compression = "lzw+p")

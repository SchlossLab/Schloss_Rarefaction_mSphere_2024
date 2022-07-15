#!/usr/bin/env Rscript

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

## Alpha diversity panel...

alpha_metrics <- c(
  "sobs_raw",
  "sobs_rarefy",
  "sobs_srs",
  "ace_raw",
  "chao_raw",
  "default_breakaway",
  "poisson_breakaway",

  "shannon_raw",
  "shannon_rarefy",
  "shannon_srs",
  "npshannon_raw",

  "simpson_raw",
  "simpson_rarefy",
  "simpson_srs"

#"ace_rarefy", "chao_rarefy", "coverage_rarefy", "npshannon_rarefy",
#"nseqs_rarefy", "npshannon_srs", "coverage_raw", "coverage_srs", "nseqs_srs",
#"nseqs_raw", "ace_srs", "chao_srs",
)

alpha_classes <- c(
  "richness",
  "richness",
  "richness",
  "richness",
  "richness",
  "richness",
  "richness",

  "shannon",
  "shannon",
  "shannon",
  "shannon",

  "inv_simpson",
  "inv_simpson",
  "inv_simpson"
)

alpha_labels <- c(
  "Raw",
  "Rarefied",
  "SRS Normalized",
  "ACE Estimate",
  "Chao1 Estimate",
  "BA Default",
  "BA Poisson",

  "Raw",
  "Rarefied",
  "SRS Normalized",
  "Estimate",

  "Raw",
  "Rarefied",
  "SRS Normalized"
)

pretty_alpha_classes <- c(
  richness = "Richness",
  shannon = "Shannon",
  inv_simpson = "Inverse\nSimpson"
)

alpha_summary_files <- glue("data/{datasets}/data.c_calpha_kw")
names(alpha_summary_files) <- datasets

alpha_composite <- map_dfr(alpha_summary_files, read_tsv, .id = "dataset") %>%
  mutate(metric_method = glue("{metric}_{method}")) %>%
  inner_join(., tibble(metric_method = alpha_metrics, class = alpha_classes),
             by = "metric_method") %>%
  mutate(metric_method = factor(metric_method, levels = alpha_metrics),
         class = factor(class, levels = names(pretty_alpha_classes))
         )

alpha <- alpha_composite %>%
  ggplot(aes(x = metric_method, y = frac, color = dataset, shape = dataset)) +
    # geom_hline(yintercept = 0.05, size = 0.25, color = "gray") +
    geom_line(aes(group = dataset), position = position_dodge(width = 0.3),
              color = "gray", size = 0.1) +
    geom_point(position = position_dodge(width = 0.3)) +
    facet_grid(. ~ class,
               scales = "free_x", space = "free_x",
               labeller = labeller(class = pretty_alpha_classes)
                ) +
    scale_color_manual(name = NULL,
                       values = rep(c("#1b9e77", "#d95f02", "#7570b3"), 4),
                       breaks = pretty_datasets$plain,
                       labels = pretty_datasets$pretty
                       ) +
    scale_shape_manual(name = NULL,
                       values = rep(c(15, 17, 19, 25), 3),
                       breaks = pretty_datasets$plain,
                       labels = pretty_datasets$pretty
                       ) +
    scale_x_discrete(breaks = alpha_metrics, labels = alpha_labels) +
    # scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.15, 0.05)) +
    labs(x = NULL, y = "Power") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(0,0,0,0),
      legend.key.size = unit(12, "pt"),
      axis.line.x = element_line()
    )

ggsave("figures/power_cffect_model.pdf", alpha, width = 6, height = 3)


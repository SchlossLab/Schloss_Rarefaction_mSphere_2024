#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
design <- args[1]

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
  "sobs_rarefy",
  "sobs_raw",
  "sobs_srs",
  "ace_raw",
  "chao_sobs_inext",
  "size_sobs_inext",
  "coverage_sobs_inext",
  "default_breakaway",
  "poisson_breakaway",

  "shannon_rarefy",
  "shannon_raw",
  "shannon_srs",
  "chao_shannon_inext",
  "size_shannon_inext",
  "coverage_shannon_inext",

  "simpson_rarefy",
  "simpson_raw",
  "simpson_srs",
  "chao_invsimpson_inext",
  "size_invsimpson_inext",
  "coverage_invsimpson_inext"
)

alpha_classes <- c(
  "richness",
  "richness",
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
  "shannon",
  "shannon",

  "inv_simpson",
  "inv_simpson",
  "inv_simpson",
  "inv_simpson",
  "inv_simpson",
  "inv_simpson"
)

alpha_labels <- c(
  "Rarefaction",
  "Raw",
  "SRS Normalized",
  "ACE Estimate",
  "Chao1 Estimate",
  "iNEXT (Size)",
  "iNEXT (Coverage)",
  "BA Default",
  "BA Poisson",

  "Rarefaction",
  "Raw",
  "SRS Normalized",
  "Estimate",
  "iNEXT (Size)",
  "iNEXT (Coverage)",

  "Rarefaction",
  "Raw",
  "SRS Normalized",
  "Estimate",
  "iNEXT (Size)",
  "iNEXT (Coverage)"
)


pretty_alpha_classes <- c(
  richness = "Richness",
  shannon = "Shannon",
  inv_simpson = "Inverse\nSimpson"
)

alpha_summary_files <- glue("data/{datasets}/data.r_{design}alpha_kw")
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
    geom_hline(yintercept = 0.05, size = 0.25, color = "gray") +
    geom_point(position = position_dodge(width = 0.3), fill = "white") +
    facet_grid(. ~ class,
               scales = "free_x", space = "free_x",
               labeller = labeller(class = pretty_alpha_classes)
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
    scale_x_discrete(breaks = alpha_metrics, labels = alpha_labels) +
    # scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.15, 0.05)) +
    labs(x = NULL, y = "False positive rate") +
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


## Beta diversity panel...
beta_metrics <- c(
  "rare_bray",
  "raw_bray",
  "relabund_bray",
  "srs_bray",
  "metagenomeseq_bray",

  "rare_jaccard",
  "raw_jaccard",
  "relabund_jaccard",
  "srs_jaccard",
  "metagenomeseq_jaccard",

  "rare_euclidean",
  "raw_euclidean",
  "relabund_euclidean",
  "rclr_euclidean",
  "oclr_euclidean",
  "nclr_euclidean",
  "zclr_euclidean",
  "deseq2_euclidean"
)

# beta_classes <- c(
#
# )

beta_labels <- c(
  "Rarefaction",
  "Raw",
  "Rel. abundance",
  "SRS Normalized",
  "CSS Normalized",

  "Rarefaction",
  "Raw",
  "Rel. abundance",
  "SRS Normalized",
  "CSS Normalized",

  "Rarefaction",
  "Raw",
  "Rel. abundance",
  "Robust CLR",
  "One CLR",
  "Nudge CLR",
  "Zero CLR",
  "DeSeq2"

)

pretty_beta_calcs <- c(
  bray = "Bray-Curtis",
  jaccard = "Jaccard",
  euclidean = "Euclidean"
)

beta_summary_files <- glue("data/{datasets}/data.r_{design}amova")
names(beta_summary_files) <- datasets

beta_composite <- map_dfr(beta_summary_files, read_tsv, .id = "dataset") %>%
  mutate(process_calc = glue("{beta_process}_{beta_calculator}")) %>%
  inner_join(., tibble(process_calc = beta_metrics),
             by = "process_calc") %>%
  mutate(process_calc = factor(process_calc, levels = beta_metrics),
         beta_calculator = factor(beta_calculator, levels = names(pretty_beta_calcs))
         )

beta <- beta_composite %>%
  drop_na(significant) %>%
  ggplot(aes(x = process_calc, y = significant,
             color = dataset, shape = dataset)) +
    geom_hline(yintercept = 0.05, size = 0.25, color = "gray") +
    geom_point(position = position_dodge(width = 0.3), show.legend = FALSE) +
    facet_grid(. ~ beta_calculator,
               scales = "free_x", space = "free_x",
              labeller = labeller(beta_calculator = pretty_beta_calcs)
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
    # scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
    labs(x = NULL, y = "False positive rate") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(12, "pt"),
      axis.line.x = element_line()
      )

combo <- alpha + beta +
  plot_layout(guides = "collect", nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold",
                            margin = margin(t = 0, r = -15, b = -15, l = 0)
                            )
    )

ggsave(glue("figures/false_positive_null_model_{design}size.tiff"),
       combo,
       width = 6, height = 7,
       compression = "lzw+p")

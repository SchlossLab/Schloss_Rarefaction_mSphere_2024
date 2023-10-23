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
  "rare_sobs",
  "raw_sobs",
  "srs_sobs",
  "raw_ace",
  "chao_sobs",
  "size_sobs",
  "coverage_sobs",
  "ba_default",
  "ba_poisson",

  "rare_shannon",
  "raw_shannon",
  "srs_shannon",
  "chao_shannon",
  "size_shannon",
  "coverage_shannon",

  "rare_simpson",
  "raw_simpson",
  "srs_simpson",
  "chao_invsimpson",
  "size_invsimpson",
  "coverage_invsimpson"
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

alpha_summary_files <- glue("data/{datasets}/random_alpha_correlation.tsv")

alpha_composite <- map_dfr(alpha_summary_files, read_tsv) %>%
  inner_join(., tibble(metric = alpha_metrics, class = alpha_classes),
             by = "metric") %>%
  mutate(metric = factor(metric, levels = alpha_metrics),
         class = factor(class, levels = names(pretty_alpha_classes)),
         median = if_else(str_detect(metric, "simpson"), -1 * median, median))

alpha <- alpha_composite %>%
  ggplot(aes(x = metric, y = median, color = dataset, shape = dataset)) +
    geom_hline(yintercept = 0, size = 0.25, color = "gray") +
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
    labs(x = NULL, y = "Correlation with\nsequencing depth") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(0,0,0,0),
      legend.key.size = unit(12, "pt"),
      axis.line.x  = element_line()
    )


## Beta diversity panel...
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
  jaccard = "Jaccard",
  bray = "Bray-Curtis",
  euclidean = "Euclidean"
)

beta_summary_files <- glue("data/{datasets}/random_beta_correlation.tsv")

beta_composite <- map_dfr(beta_summary_files, read_tsv) %>%
  mutate(process_calc = glue("{process}_{calculator}")) %>%
  inner_join(., tibble(process_calc = beta_metrics),
             by = "process_calc") %>%
  mutate(process_calc = factor(process_calc, levels = beta_metrics),
         calculator = factor(calculator, levels = names(pretty_beta_calcs))
         )

beta <- beta_composite %>%
  ggplot(aes(x = process_calc, y = median, color = dataset, shape = dataset)) +
      geom_hline(yintercept = 0, size = 0.25, color = "gray") +
    geom_point(position = position_dodge(width = 0.3), show.legend = FALSE) +
    facet_grid(. ~ calculator,
               scales = "free_x", space = "free_x",
              labeller = labeller(calculator = pretty_beta_calcs)
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
    labs(x = NULL, y = "Correlation with difference\nin sequencing depth") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(12, "pt"),
      axis.line.x  = element_line()
      )


combo <- alpha + beta +
  plot_layout(guides = "collect", nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold",
                            margin = margin(t = 0, r = -15, b = -15, l = 0)
                            )
    )

ggsave("figures/alpha_beta_depth_correlation.tiff",
       combo,
       width = 6, height = 7,
       compression = "lzw+p")

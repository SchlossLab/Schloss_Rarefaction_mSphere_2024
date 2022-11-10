#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
# library(patchwork)
library(ggtext)

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

obs_coverage_files <- glue("data/{datasets}/data.otu.obs_coverage")
names(obs_coverage_files) <- datasets

obs_cor <- map_dfr(obs_coverage_files, read_tsv, .id = "dataset",
        col_types = cols(group = col_character(),
                         .default = col_double()))


# datasets <- datasets[datasets != "rice" ]
rare_coverage_files <- glue("data/{datasets}/data.otu.rarefy_coverage")
names(rare_coverage_files) <- datasets

cor_line <- map_dfr(rare_coverage_files, read_tsv, .id = "dataset",
        col_types = cols(.default = col_double())) %>%
  filter(nseqs != 1 & nseqs <= 1e6)


cor_plot <- obs_cor %>%
  ggplot(aes(x = norare_nseqs, y = norare_coverage,
             color = dataset, shape = dataset)) +
  geom_point(fill = "white") +
  geom_boxplot(aes(x=325, y=rare_coverage),
               outlier.shape = NULL, outlier.size = 0.75,
               width = 0.2) +
  geom_line(data = cor_line, aes(x = nseqs, y = mean), color = "black") +
  geom_vline(xintercept = 800, color = "darkgray") +
  scale_x_log10(breaks = c(325, 1e3, 1e4, 1e5, 1e6),
                labels = c("R", "10^3^", "10^4^", "10^5^", "10^6^"),
                limits = c(2e2, 1e6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  scale_color_manual(name = NULL,
                    values = rep(c('#1b9e77','#d95f02','#7570b3'), 4),
                    breaks = pretty_datasets$plain,
                    labels = pretty_datasets$pretty,
                    guide = "none"
                    ) +
  scale_shape_manual(name = NULL,
                  values = rep(c(15, 17, 19, 25), 3),
                  breaks = pretty_datasets$plain,
                  labels = pretty_datasets$pretty,
                    guide = "none"
                  ) +

  facet_wrap(~dataset, labeller = labeller(.default = capitalize)) +
  labs(y = "Good's coverage (%)",
       x = "Number of sequences per sample") +
  theme(
    axis.text.x = element_markdown(),
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.25, linetype = "dotted", color="gray"),
    axis.line = element_line(),
  )
  
ggsave("figures/coverage_plot.tiff", cor_plot,
       width = 7, height = 6,
       compression = "lzw+p")

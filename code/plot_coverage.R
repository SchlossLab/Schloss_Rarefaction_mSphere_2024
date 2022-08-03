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

cor <- map_dfr(obs_coverage_files, read_tsv, .id = "dataset",
        col_types = cols(group = col_character(),
                         .default = col_double())) %>%
  ggplot(aes(x = norare_nseqs, y = norare_coverage, color = dataset, shape = dataset)) +
  geom_point(fill = "white", show.legend = FALSE) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6),
                labels = c("10^3", "10^4", "10^5", "10^6"),
                limits = c(1e3, 1e6)) +
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

  facet_wrap(~dataset, labeller = labeller(.default = capitalize)) +
  labs(y = "Good's coverage",
       x = "Number of sequences per sample") +
  theme(
    axis.text.x = element_markdown(),
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.25, linetype = "dotted", color="gray"),
    axis.line = element_line(),
    
  )
  
ggsave("figures/obs_coverage_plot.pdf", cor, width = 7, height = 6)



# map_dfr(obs_coverage_files, read_tsv, .id = "dataset",
#         col_types = cols(group = col_character(),
#                          .default = col_double())) %>%
#   mutate(rare_norare = rare_coverage / norare_coverage) %>%
#   group_by(dataset) %>%
#   summarize(mean_ratio = mean(rare_norare),
#             sd_ratio = sd(rare_norare))
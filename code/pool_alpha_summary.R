#!/usr/bin/env Rscript

library(tidyverse)

input <- commandArgs(trailingOnly = TRUE)

summary_files <- input
# summary_files <- list.files(path = "data/soil",
#                             pattern = ".*r_rare_alpha.summary",
#                             full.names = TRUE)

output_file <- str_replace(summary_files[1],
                           "data.otu.*",
                           "data.otu.alpha_depth.summary")

reader <- function(x) {

  seed <- str_replace(x, ".*\\.(\\d*)\\.r_rare_.*.summary", "\\1")
  dataset <- str_replace(x, "data/(.*)/data.otu.*", "\\1")
  # dmethod <- str_replace(x, ".*\\.\\d*\\.r_rare_(.*).summary", "\\1")

  read_tsv(x,
           col_types = cols(metric = col_character(),
                            .default = col_double())) %>%
  mutate(seed = seed, dataset = dataset) %>%
  select(seed, dataset, method = metric, n_seqs, mean, sd, n_samples)

}

map_dfr(summary_files, reader) %>%
  write_tsv(output_file)
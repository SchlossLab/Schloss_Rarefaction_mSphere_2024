#!/usr/bin/env Rscript

library(tidyverse)

input <- commandArgs(trailingOnly = TRUE)

summary_files <- input
# summary_files <- list.files(path = "data/mice",
#                             pattern = ".*r_rare_.*.summary",
#                             full.names = TRUE)

output_file <- str_replace(summary_files[1],
                           "data.otu.*",
                           "data.otu.beta_depth.summary")

reader <- function(x) {

  seed <- str_replace(x, ".*\\.(\\d*)\\.r_rare_.*.summary", "\\1")
  dataset <- str_replace(x, "data/(.*)/data.otu.*", "\\1")
  dmethod <- str_replace(x, ".*\\.\\d*\\.r_rare_(.*).summary", "\\1")

  read_tsv(x,
           col_types = cols(.default = col_double())) %>%
<<<<<<< HEAD
    mutate(seed = as.numeric(seed), dataset = dataset, method = dmethod,
           .before = n_seqs)
=======
  mutate(seed = seed, dataset = dataset, method = dmethod, .before = n_seqs)
>>>>>>> paper

}

map_dfr(summary_files, reader) %>%
<<<<<<< HEAD
  arrange(seed, method, n_seqs) %>%
=======
>>>>>>> paper
  write_tsv(output_file)
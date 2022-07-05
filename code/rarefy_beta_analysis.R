#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)
# library(furrr)

# plan(multicore, workers = 8)

input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/mice/data.otu.1.rshared"
nseqs_file <- input[2] # nseqs_file <- "data/mice/data.group_count"
output_file <- input[3] # output_file <- "data/mice/data.otu.1.r_rare_bray.summary"

shared_df <- fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames("Group")

groups <- rownames(shared_df)

max_nseqs <- read_tsv(nseqs_file,
                  col_names = c("Group", "nseqs")) %>%
  filter(Group %in% groups) %>%
  arrange(nseqs) %>% 
  slice(n = nrow(.) - 4) %>% # insures there are at least 10 distances in 
  pull(nseqs)                     # final class


nseqs <- c(1000, seq(2000, max_nseqs, 2000))

distance <- str_replace(output_file, ".*r_rare_(.*)\\.summary", "\\1")

run_rarefy_beta <- function(df, x, distance) {

  avgdist(df, sample = x, dmethod = distance)  %>%
      as.matrix() %>%
      as_tibble(rownames = "samples") %>%
      pivot_longer(cols = -samples) %>%
      filter(samples < name) %>%
      summarize(mean = mean(value),
                sd = sd(value),
                n_samples = (1 + sqrt(1 + n() * 8)) / 2) %>%
      mutate(n_seqs = x, .before = mean)

}

#future_map_dfr(nseqs,
map_dfr(nseqs,
               ~run_rarefy_beta(df = shared_df,
                                x = .x,
                                distance = distance)) %>%
  write_tsv(output_file)
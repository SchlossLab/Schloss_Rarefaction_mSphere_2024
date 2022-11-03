#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)
library(glue)
library(magrittr)


input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/mice/data.otu.1.rshared"
nseqs_file <- input[2] # nseqs_file <- "data/mice/data.group_count"
remove_file <- input[3] # remove_file <- "data/mice/data.remove_accnos"
output_file <- input[4] # output_file <- "data/mice/data.otu.1.1.r_rare_bray.summary"

index <- str_replace(output_file,
            "data/.*/data\\.otu\\.(\\d*)\\.\\d*\\.r_.*summary",
            "\\1") %>%
  as.numeric()

to_remove <- read_tsv(remove_file, col_names="Group")

max_nseqs <- read_tsv(nseqs_file,
                      col_types = cols(Group = col_character()),
                      col_names = c("Group", "nseqs")) %>%
  anti_join(., to_remove) %>%
  arrange(nseqs) %>% 
  slice(n = nrow(.) - 4) %>% # insures there are at least 10 distances in 
  pull(nseqs)                     # final class



nseqs <- 10^seq(log10(1000), log10(max_nseqs), length.out = 100) %>%
  as.integer() %>%
  extract(index)

distance <- str_replace(output_file, ".*r_rare_(.*)\\.summary", "\\1")
distance <- case_when(distance == "bray" ~ "braycurtis",
            distance == "jaccard" ~ "jclass")

temp_shared <- str_replace(output_file, "summary", "shared")
file.copy(shared_file, temp_shared, overwrite = TRUE)

mothur_cmd <- glue("mothur '#dist.shared(shared={temp_shared}, subsample={nseqs},
                                         calc={distance}, iters=100, processors=1, output=square)'")
system(mothur_cmd)

temp_dist <- str_replace(temp_shared, "shared", glue("{distance}.0.03.square.ave.dist"))

read_tsv(temp_dist, skip = 1, col_names = FALSE) %>%
  set_colnames(c("samples", all_of(.[["X1"]]))) %>%
      pivot_longer(cols = -samples) %>%
      filter(samples < name) %>%
      summarize(mean = mean(value),
                sd = sd(value),
                n_samples = (1 + sqrt(1 + n() * 8)) / 2) %>%
      mutate(n_seqs = nseqs, .before = mean) %>%
      write_tsv(output_file)
     
file.remove(temp_shared)
file.remove(temp_dist)
file.remove(str_replace(temp_dist, "square.ave.dist", "square.std.dist"))

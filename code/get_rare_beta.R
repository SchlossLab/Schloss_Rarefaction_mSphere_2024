#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)

set.seed(19760620)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] #e.g. dist_file <- "data/soil/data.otu.1.r_rare_jaccard.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

tidy_shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu"))
  
min_n_seqs <- tidy_shared %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(n = sum(value)) %>%
  summarize(n = min(n)) %>%
  pull(n)
  
tidy_shared %>%
  column_to_rownames("Group") %>%
  avgdist(sample = min_n_seqs, dmethod = method) %>%
  as.matrix() %>%
  as_tibble(rownames = "Group") %>%
  write_tsv(dist_file)

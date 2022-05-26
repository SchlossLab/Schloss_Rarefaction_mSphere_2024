#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)

group_count_file <- input[1] # e.g. group_count_file <- "data/mice/data.group_count"
remove_accnos_file <- input[2] #e.g. remove_accnos_file <- "data/mice/data.remove_accnos"
seed <- as.numeric(input[3]) # e.g. seed <- 87

samplesize_design_file <- str_replace(group_count_file,
                             "group_count",
                             glue("{seed}.sdesign"))

remove_accnos <- scan(remove_accnos_file, sep = "\n", what = character())

group_count <- read_tsv(group_count_file,
         col_names = c("Group", "n_seqs"),
         col_types = cols(Group = col_character(), n_seqs = col_double())
         ) %>%
  filter(! Group %in% remove_accnos)

n_samples <- nrow(group_count)

n_a <- floor(n_samples / 2)
n_b <- n_samples - n_a

group_count %>%
  arrange(n_seqs) %>%
  mutate(treatment = (c(rep("A", n_a), rep("B", n_b)))) %>%
  select(Group, treatment) %>%
  write_tsv(., samplesize_design_file)

#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)

group_count_file <- input[1] # e.g. group_count_file <- "data/human/data.group_count"
remove_accnos_file <- input[2] #e.g. remove_accnos_file <- "data/human/data.remove_accnos"
seed <- as.numeric(input[3]) # e.g. seed <- 87

set.seed(seed)

samplesize_design_file <- str_replace(group_count_file,
                             "group_count",
                             glue("{seed}.idesign"))

remove_accnos <- scan(remove_accnos_file, sep = "\n",
                      what = character(), quiet = TRUE)

group_count <- read_tsv(group_count_file,
         col_names = c("Group", "n_seqs"),
         col_types = cols(Group = col_character(), n_seqs = col_double())
         ) %>%
  filter(! Group %in% remove_accnos)

n_samples <- nrow(group_count)

l_quantile <- quantile(group_count$n_seqs, prob = 0.10)
u_quantile <- quantile(group_count$n_seqs, prob = 0.90)

group_count %>%
  mutate(treatment = case_when(n_seqs < l_quantile ~ "A",
                              n_seqs > u_quantile ~ "B")) %>%
  nest(data = -treatment) %>%
  mutate(rand_treat = map2(data, treatment,
                          function(.x, .y) {

                            if (is.na(.y)) {
                              return(sample(rep(c("A", "B"),
                                                length.out = nrow(.x))))
                            }

                            return(.y)

                          })) %>%
  unnest(c(data, rand_treat)) %>%
  select(Group, treatment = rand_treat) %>%
  write_tsv(., samplesize_design_file)

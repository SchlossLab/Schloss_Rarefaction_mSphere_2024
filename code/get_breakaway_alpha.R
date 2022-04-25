#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(breakaway)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/soil/data.otu.1.rshared"

breakaway_file <- str_replace(shared_file,
                             "shared",
                             "_breakaway_alpha")


get_breakaway <- function(x) {

  ba <- breakaway::breakaway(x)

  tibble(ba_est = ba$estimate,
         ba_lci = ba$interval[1], ba_uci = ba$interval[2],
         ba_model = ba$model)

}

get_poisson <- function(x) {

  p <- breakaway::poisson_model(x)

  tibble(p_est = p$estimate,
         p_lci = p$interval[1], p_uci = p$interval[2],
         p_model = p$model)

}

fread(shared_file,
      colClasses = c(Group = "character")) %>%
  pivot_longer(cols = -c(label, Group, numOtus),
               names_to = "otu",
               values_to = "n_seqs") %>%
  filter(n_seqs != 0) %>%
  count(Group, n_seqs) %>%
  nest(data = -Group) %>%
  mutate(n_seqs = map(data, ~sum(.x$n * .x$n_seqs)),
         sobs = map(data, ~sum(.x$n)),
         ba = map(data, ~get_breakaway(.x)),
         p = map(data, ~get_poisson(.x))) %>%
  unnest(c(n_seqs, sobs, ba, p)) %>%
  select(-data) %>%
  write_tsv(breakaway_file)

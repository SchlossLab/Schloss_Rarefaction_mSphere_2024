#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)

read_raw <- function(x) {

# x <- "data/human/data.otu.1.r_raw_alpha" 

  read_tsv(x, col_types = cols(group = col_character())) %>%
    select(group, nseqs,
           raw_sobs = sobs, raw_shannon = shannon,
           raw_simpson = simpson, raw_chao = chao,
           raw_ace = ace, raw_npshannon = npshannon,
           raw_coverage = npshannon)
}


read_rarefy <- function(x) {

# x <- "data/human/data.otu.1.r_rarefy_alpha" 

  read_tsv(x, col_types = cols(group = col_character())) %>%
    filter(method == "ave") %>%
    select(group, rare_nseqs = nseqs,
           rare_sobs = sobs, rare_shannon = shannon,
           rare_simpson = simpson, rare_chao = chao,
           rare_ace = ace, rare_npshannon = npshannon,
           rare_coverage = npshannon)

}



read_srs <- function(x) {

# x <- "data/human/data.otu.1.r_srs_alpha"

  read_tsv(x, col_types = cols(group = col_character())) %>%
    select(group, srs_nseqs = nseqs,
           srs_sobs = sobs, srs_shannon = shannon,
           srs_simpson = simpson, srs_chao = chao,
           srs_ace = ace, srs_npshannon = npshannon,
           srs_coverage = npshannon)

}


read_breakaway <- function(x) {

# x <- "data/human/data.otu.1.r_breakaway_alpha"

  read_tsv(x, col_types = cols(Group = col_character())) %>%
    select(group = Group, ba_nseqs = n_seqs,
           ba_sobs = sobs, ba_default = ba_est, ba_poisson = p_est)

}


input <- commandArgs(trailingOnly = TRUE)
dataset <- input[1] # e.g. dataset <- "human"

seeds <- 1:100

# methods <- c("raw", "rarefy", "srs", "breakaway")

output_file <- glue("data/{dataset}/random_alpha_correlation.tsv")


raw <- map_dfr(glue("data/{dataset}/data.otu.{seeds}.r_raw_alpha"),
               read_raw, .id = "seed")

rarefy <- map_dfr(glue("data/{dataset}/data.otu.{seeds}.r_rarefy_alpha"),
                  read_rarefy, .id = "seed")

srs <- map_dfr(glue("data/{dataset}/data.otu.{seeds}.r_srs_alpha"),
               read_srs, .id = "seed")

breakaway <- map_dfr(glue("data/{dataset}/data.otu.{seeds}.r_breakaway_alpha"),
               read_breakaway, .id = "seed")

stopifnot(sum(raw$nseqs != breakaway$ba_nseqs) == 0)
stopifnot(sum(rarefy$rare_nseqs != srs$srs_nseqs) == 0)

inner_join(raw, rarefy, by = c("seed", "group")) %>%
  inner_join(., srs, by = c("seed", "group")) %>%
  inner_join(., breakaway, by = c("seed", "group")) %>%
  select(-rare_nseqs, -srs_nseqs, -ba_nseqs) %>%
  pivot_longer(cols = -c(seed, group, nseqs),
               names_to = "metric", values_to = "value") %>%
  nest(data = -c(metric, seed)) %>%
  mutate(test = map(data, ~cor.test(.x$nseqs, .x$value, method = "spearman",
                                    exact = FALSE) %>%
                    tidy())) %>%
  unnest(test) %>%
  select(seed, metric, estimate, p.value) %>%
  group_by(metric) %>%
  summarize(median = median(estimate),
            iqr = IQR(estimate),
            frac_sig = sum(p.value < 0.05)/n()) %>%
  mutate(dataset = dataset) %>%
  write_tsv(output_file)

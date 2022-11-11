#!/usr/bin/env Rscript

library(tidyverse)
library(cluster)
library(glue)

get_accuracy <- function(seed, model, design, dataset, correct, dist) {

  data <- read_tsv(glue("data/{dataset}/data.otu.{seed}.{model}_{correct}_{dist}.dist"),
          col_types = cols(Group = col_character())) %>%
    column_to_rownames("Group")
  
  if(is.na(data[1,2])) {
    return(NA_real_)
  }
  
  data %>%
    pam(k = 2, diss = T) %>%
    pluck("clustering") %>%
    as_tibble(rownames = "Group") %>%
    rename(cluster = value) %>%
    inner_join(., read_tsv(glue("data/{dataset}/data.{seed}.{design}design"),
                          col_types = cols(Group = col_character())),
              by = "Group") %>%
    count(cluster, treatment) %>%
    pivot_wider(names_from = "treatment", values_from = "n", values_fill = 0) %>%
    summarize(a = .[["A"]][1] + .[["B"]][2], b = .[["A"]][2] + .[["B"]][1]) %>%
    pivot_longer(cols = c("a", "b")) %>%
    mutate(accuracy = value / sum(value)) %>%
    slice_max(accuracy, with_ties = FALSE) %>%
    pull(accuracy)

}

input <- commandArgs(trailingOnly = TRUE)
dataset <- input[1]

output_file <- glue("data/{dataset}/cluster_analysis.tsv")

corrections <- c("raw", "relabund", "rare",
                 "deseq2", "metagenomeseq", "srs",
                 "nclr", "oclr", "rclr", "zclr")

distances <- c("bray",
               "euclidean",
               "jaccard")

model <- c("r", "e")

design <- c("r", "s", "e")

seeds <- 1:100

expand_grid(model, design, corrections, distances, seeds) %>% 
  filter((model == "r" & design != "e") | (model == "e" & design == "e")) %>%
  mutate(
      dataset = dataset
  ) %>%
  rowwise() %>%
  mutate(
    accuracies = get_accuracy(seeds, model, design, dataset, corrections, distances)
  ) %>%
  write_tsv(output_file)


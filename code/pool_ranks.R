#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

datasets <- c("bioethanol", "human", "lake", "marine", "mice", "peromyscus",
            "rainforest", "sediment", "soil", "stream", "rice", "seagrass")

read_tsv(glue("data/{datasets}/data.rank")) %>%
  group_by(dataset, simulation) %>%
  summarize(frac_full_rank = mean(rank == n_samples), .groups = "drop") %>%
  write_tsv("data/process/rank_fractions.tsv")

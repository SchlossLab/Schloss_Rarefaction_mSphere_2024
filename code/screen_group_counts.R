#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

data_threshold <- tribble(
  ~data_set,   ~threshold,
  "bioethanol", 3600,
  "human",      10000,
  "lake",       10000,
  "marine",     3600,
  "mice",       1800,
  "peromyscus", 4000,
  "rainforest", 3600,
  "rice",       2500,
  "seagrass",   1800,
  "sediment",   7000,
  "soil",       3600,
  "stream",     3600
)

get_removes <- function(d, t) {

  glue("data/{d}/data.group_count") %>%
    read_tsv(col_names = c("group", "n_seqs"),
             col_types = cols(col_character(), col_double())) %>%
    filter(n_seqs < t) %>%
    pull(group) %>%
    write(glue("data/{d}/data.remove_accnos"))

}

map2(data_threshold$data_set, data_threshold$threshold, get_removes)
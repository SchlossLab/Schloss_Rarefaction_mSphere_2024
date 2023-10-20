#!/usr/bin/env Rscript

library(data.table)
library(Matrix)
library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
dataset <- input[1]

get_matrix_stats <- function(file_name) {

  fread(file_name, colClasses = c(Group = "character")) %>%
      select(Group, starts_with("Otu")) %>%
      column_to_rownames("Group") %>% 
      summarize(rank = rankMatrix(as.matrix(.)),
                n_samples = nrow(.),
                n_otus = ncol(.))

}

shared_file <- glue("data/{dataset}/data.otu.1.eshared")

r_seeds <- rep(1:100, 3)
r_simulations <- rep(c("rshared", "eshared", "cshared"), each = 100)

shared_files <- c(glue("data/{dataset}/data.otu.{r_seeds}.{r_simulations}"),
                  glue("data/{dataset}/data.otu.shared"))

map_dfr(shared_files, get_matrix_stats) %>%
  mutate(dataset = dataset,
        seed = c(r_seeds, 0),
        simulation = c(r_simulations, "observed"),
        .before = rank) %>%
  write_tsv(glue("data/{dataset}/data.rank"))

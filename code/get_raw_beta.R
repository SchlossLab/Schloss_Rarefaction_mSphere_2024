#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] #e.g. dist_file <- "data/soil/data.otu.1.r_raw_bray.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames("Group") %>%
  vegdist(method = method) %>%
  as.matrix() %>%
  as_tibble(rownames = "Group") %>%
  write_tsv(dist_file)

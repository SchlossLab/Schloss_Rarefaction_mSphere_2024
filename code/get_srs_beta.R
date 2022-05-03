#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(SRS)
library(vegan)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] # dist_file <- "data/soil/data.otu.1.r_srs_bray.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu"))

smallest_group <- shared %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(n = sum(value)) %>%
  summarize(min = min(n)) %>%
  pull(min)

shared %>%
  column_to_rownames("Group") %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin = smallest_group) %>%
  t() %>%
  vegdist(method = method) %>%
  as.matrix() %>%
  as_tibble(rownames = "Group") %>%
  write_tsv(dist_file)

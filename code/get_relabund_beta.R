#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)

set.seed(19760620)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] # dist_file <- "data/soil/data.otu.1.r_relabund_bray.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(n = sum(value),
         rel_abund = value / n) %>%
  ungroup() %>%
  select(-value, -n) %>%
  pivot_wider(names_from = name, values_from = rel_abund, values_fill = 0) %>%  
  column_to_rownames("Group") %>%
  vegdist(method = method) %>%
  as.matrix() %>%
  as_tibble(rownames = "Group") %>%
  write_tsv(dist_file)

#!/usr/bin/env Rscript

library(metagenomeSeq)
library(tidyverse)
library(data.table)
library(vegan)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] # dist_file <- "data/soil/data.otu.1.r_metagenomeseq_bray.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")


fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  filter(value > 0) %>%
  pivot_wider(names_from = "name", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Group") %>%
  t() %>% 
  newMRexperiment(counts = .) %>% # samples need to be in columns and otus need
  MRcounts(norm = TRUE) %>%       # to be in rows
  t() %>%
  vegdist(method = method) %>%
  as.matrix() %>%
  as_tibble(rownames = "Group") %>%
  write_tsv(dist_file)

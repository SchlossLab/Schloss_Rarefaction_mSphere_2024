#!/usr/bin/env Rscript

library(DESeq2)
library(tidyverse)
library(data.table)
library(vegan)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] # dist_file <- "data/soil/data.otu.1.r_deseq2_bray.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

if (method == "euclidean") {

  fread(shared_file, colClasses = c(Group = "character")) %>%
    select(Group, starts_with("Otu")) %>%
    pivot_longer(-Group) %>%
    filter(value > 0) %>%
    mutate(value = value + 1) %>%  # can't have zero counts
    pivot_wider(names_from = "name", values_from = "value", values_fill = 0) %>%
    column_to_rownames("Group") %>%
    t() %>%
    varianceStabilizingTransformation(fitType = "mean") %>%
    t() %>%
    vegdist(method = method) %>%
    as.matrix() %>%
    as_tibble(rownames = "Group") %>%
    write_tsv(dist_file)

} else {

  fread(shared_file, colClasses = c(Group = "character")) %>%
    select(Group) %>%
    expand_grid(B = Group) %>%
    mutate(value = if_else(Group == B, 0, NA_real_)) %>%
    pivot_wider(values_from = value, names_from = B) %>%
    write_tsv(dist_file)

}
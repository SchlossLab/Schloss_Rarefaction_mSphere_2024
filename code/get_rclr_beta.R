#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(vegan)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"
dist_file <- input[2] # dist_file <- "data/soil/data.otu.1.r_rclr_euclidean.dist"

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")

geometric_mean <- function(x) {

  exp(mean(log(x[x > 0])))

}

if (method == "euclidean") {

  fread(shared_file, colClasses = c(Group = "character")) %>%
    select(Group, starts_with("Otu")) %>%
    pivot_longer(-Group) %>%
    filter(value > 0) %>%
    group_by(Group) %>%
    mutate(rclr = log(value / geometric_mean(value))) %>%
    ungroup() %>%
    select(-value) %>%
    pivot_wider(names_from = name, values_from = rclr, values_fill = 0) %>%
    column_to_rownames("Group") %>%
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
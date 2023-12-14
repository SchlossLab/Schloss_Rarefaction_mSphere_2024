#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1]
dist_file <- input[2]
# shared_file <- "data/soil/data.otu.1.rshared"
# dist_file <- "data/soil/data.otu.1.r_rclr_euclidean.dist"

temp_otu_file <- str_replace(shared_file, "shared", "shared_gemelli")
temp_dist_file <- str_replace(shared_file, "shared", "dist_gemelli")

method <- str_replace(dist_file, ".*_(\\w*).dist", "\\1")
dataset <- str_replace(dist_file, "data/(.*)/data.*", "\\1")

# couldn't get stream, rice, seagrass to reliably get through run_gemelli even
# with 100 GB of RAM

if (method == "euclidean" && !(dataset %in% c("stream", "rice", "seagrass"))) {

  fread(shared_file, colClasses = c(Group = "character")) %>%
    select(Group, starts_with("Otu")) %>%
    pivot_longer(-Group) %>%
    pivot_wider(names_from = Group, values_from = value, values_fill = 0) %>%
    column_to_rownames("name") %>%
    write.table(temp_otu_file, sep = "\t", quote = FALSE)

  system(glue("./code/run_gemelli.py {temp_otu_file} {temp_dist_file}"))

  fread(temp_dist_file, header = TRUE) %>%
    as_tibble() %>%
    rename(Group = V1) %>%
    write_tsv(dist_file)

  unlink(temp_otu_file, force = TRUE)
  unlink(temp_dist_file, force = TRUE)

} else {

  fread(shared_file, colClasses = c(Group = "character")) %>%
    select(Group) %>%
    expand_grid(B = Group) %>%
    mutate(value = if_else(Group == B, 0, NA_real_)) %>%
    pivot_wider(values_from = value, names_from = B) %>%
    write_tsv(dist_file)

}
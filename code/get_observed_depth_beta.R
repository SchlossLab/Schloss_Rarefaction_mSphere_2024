#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(vegan)
library(data.table)


input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/human/data.otu.shared"
depth <- as.numeric(input[2]) # depth <- 1000
distance <- input[3] # distance <- "bray"

dist_file <- str_replace(shared_file,
                         "shared",
                         glue("obs_{depth}_{distance}.dist"))


fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames("Group") %>%
  avgdist(sample = depth, dmethod = distance)  %>%
      as.matrix() %>%
      as_tibble(rownames = "Group") %>%
      write_tsv(dist_file)

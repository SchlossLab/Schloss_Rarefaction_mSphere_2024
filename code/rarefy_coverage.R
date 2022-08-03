#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/mice/data.otu.shared"
remove_file <- input[2] # e.g. remove_file <- "data/mice/data.remove_accnos"

dataset <- str_replace(shared_file,
                        "data/(.*)/data.otu.shared",
                        "\\1")

composite_shared_file <- str_replace(shared_file,
                                     "shared",
                                     "composite.shared")

removal <- scan(remove_file, quiet = TRUE, what = character())

fread(shared_file, colClasses = c(Group = "character")) %>%
  pivot_longer(cols = -c(label, Group, numOtus),
               names_to = "otu",
               values_to = "n_seqs") %>%
  filter(! Group %in% removal) %>%
  select(Group, otu, n_seqs) %>%
  group_by(otu) %>%
  summarize(n_seqs = sum(n_seqs)) %>%
  pivot_wider(names_from = "otu", values_from = "n_seqs") %>%
  mutate(label = 0.03,
         Group = "composite",
         numOtus = ncol(.),
         .before = everything()) %>%
  write_tsv(composite_shared_file)
  
mothur_cmd <- glue("mothur '#rarefaction.single(shared={composite_shared_file},
                                                calc=coverage, freq=1000,
                                                processors=1)'")

system(mothur_cmd)

mothur_file <- glue("data/{dataset}/data.otu.composite.groups.r_coverage")
output_file <- glue("data/{dataset}/data.otu.rarefy_coverage")

read_tsv(mothur_file,
         skip  = 1,
         col_names = c("nseqs", "mean", "lci", "uci")) %>%
  write_tsv(output_file)

file.remove(composite_shared_file)
file.remove(mothur_file)
file.remove(glue("data/{dataset}/data.otu.composite.composite.rabund"))
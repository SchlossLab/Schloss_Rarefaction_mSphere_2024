#!/usr/bin/env Rscript

# this script takes in a shared file and random number generator seed and
# outputs a null model version of the shared file. the samples in the outputted
# shared file have the same number of sequences that they had initially, but
# are a statistical sampling of the pooled overal abundances of each OTU

library(tidyverse)
library(data.table)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/soil/data.otu.shared"
seed <- as.numeric(input[2]) # e.g. seed <- 1

# e.g. r_shared_file <- "data/soil/data.otu.1.rshared"
r_shared_file <- str_replace(shared_file,
                             "shared",
                             glue("{seed}.rshared"))


#	set the random number generator seed so multiple runs generate the same 
#	randomization
set.seed(seed)

rand_shared <- fread(shared_file,
                     colClasses = c(Group = "character")) %>%
  pivot_longer(cols = -c(label, Group, numOtus),
               names_to = "otu",
               values_to = "n_seqs") %>%
  uncount(n_seqs) %>%
  mutate(otu = sample(otu)) %>%
  count(label, Group, numOtus, otu) %>%
  pivot_wider(names_from = otu, values_from = n, values_fill = 0)

stopifnot(ncol(rand_shared) - 3 == unique(rand_shared$numOtus))

write_tsv(rand_shared, file = r_shared_file)
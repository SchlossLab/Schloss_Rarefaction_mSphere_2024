#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/lake/data.otu.shared"
remove_file <- input[2] # remove_file <- "data/lake/data.remove_accnos"

output_file <- str_replace(shared_file, "shared", "obs_coverage")
temp_shared_file <- str_replace(shared_file, "shared", "temp_coverage.shared")

remove_groups <- scan(remove_file, what = character(), quiet = T)

fread(shared_file) %>%
  filter(!(Group %in% remove_groups)) %>%
  write_tsv(temp_shared_file)

mothur_command <- glue("mothur '#summary.single(shared={temp_shared_file},
                                                calc=nseqs-coverage,
                                                subsample = TRUE);
                       summary.single(shared={temp_shared_file},
                                      calc=nseqs-coverage)'")

system(mothur_command)

rare_file <- str_replace(temp_shared_file, "shared", "groups.ave-std.summary") %>%
  read_tsv() %>%
  filter(method == "ave") %>%
  select(group, rare_nseqs = nseqs, rare_coverage = coverage)
  
norare_file <- str_replace(temp_shared_file, "shared", "groups.summary") %>%
  read_tsv() %>%
  select(group, norare_nseqs = nseqs, norare_coverage = coverage)
  
inner_join(rare_file, norare_file, by = "group") %>%
  write_tsv(output_file)
  

temp_files <- list.files(path = dirname(temp_shared_file),
                         pattern = ".*temp_coverage.*",
                         full.names = TRUE)
unlink(temp_files)

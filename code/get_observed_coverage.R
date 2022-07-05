#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/marine/data.otu.shared"
output_file <- str_replace(shared_file, "shared", "obs_coverage")

temp_shared_file <- str_replace(shared_file, "shared", "temp_coverage.shared")

if (file.exists(temp_shared_file)) {
  unlink(temp_shared_file)
}

file.copy(shared_file, temp_shared_file, overwrite = TRUE)

mothur_command <- glue("mothur '#summary.single(shared={temp_shared_file},
                                                calc=nseqs-coverage)'")

system(mothur_command)


str_replace(temp_shared_file, "shared", "groups.summary") %>%
  read_tsv() %>%
  select(-label) %>%
  write_tsv(output_file)
  
  
temp_files <- list.files(path = dirname(temp_shared_file),
                         pattern = ".*temp_coverage.*",
                         full.names = TRUE)
unlink(temp_files)

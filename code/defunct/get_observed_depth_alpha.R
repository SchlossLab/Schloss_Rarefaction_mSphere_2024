#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/human/data.otu.shared"
depth <- as.numeric(input[2]) # depth <- 1000

mothur_output_file <-
      str_replace(shared_file,
                  "shared",
                  glue("temp_obs_{depth}_alpha.groups.ave-std.summary"))

script_output_file <- str_replace(shared_file,
                                  "shared",
                                  glue("obs_{depth}_alpha.summary"))

temp_shared_file <- str_replace(shared_file, "shared",
                                glue("temp_obs_{depth}_alpha.shared"))

if (file.exists(temp_shared_file)) {
  unlink(temp_shared_file)
}

file.copy(shared_file, temp_shared_file, overwrite = TRUE)




mothur <- glue("mothur '#summary.single(shared = {temp_shared_file},
                        calc=sobs-shannon-simpson-coverage,
                        subsample={depth})'")

system(mothur)

file.rename(mothur_output_file, script_output_file)

temp_files <- list.files(path = dirname(temp_shared_file),
                         pattern = str_replace(basename(temp_shared_file),
                                               ".shared", ""),
                         full.names = TRUE)

unlink(temp_files)

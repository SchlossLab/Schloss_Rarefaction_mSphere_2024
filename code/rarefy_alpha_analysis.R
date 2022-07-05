#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

input <- commandArgs(trailingOnly = TRUE)

shared_file <- input[1] # shared_file <- "data/soil/data.otu.24.rshared"
output_file <- input[2] # output_file <- "data/soil/data.otu.24.r_rare_alpha.summary"


temp_shared_file <- str_replace(shared_file, "rshared", "temp_alpha.rshared")

if (file.exists(temp_shared_file)) {
  unlink(temp_shared_file)
}

file.copy(shared_file, temp_shared_file, overwrite = TRUE)

mothur <- glue("mothur '#rarefaction.single(shared = {temp_shared_file},
                        calc=sobs-shannon-simpson-coverage,
                        freq=1000,
                        processors = 1)'")

system(mothur)

read_rarefaction_files <- function(x) {

  m <- if_else(str_detect(x, "rarefaction"),
               "sobs",
               str_replace(x, ".*temp_alpha.groups.r_", "")
               )
  
  read_tsv(x, col_types = cols(.default = col_double())) %>%
    filter(numsampled == 1000 | (numsampled %% 2000 == 0)) %>%
    select(n_seqs = numsampled, starts_with("0.03")) %>%
    pivot_longer(-n_seqs, names_to = "group", values_to = "value") %>%
    mutate(group = str_replace(group, "0\\.03-", ""),
           metric = m)

}

extensions <- c("r_coverage", "r_shannon", "r_simpson", "rarefaction")
file_stub <- str_replace(temp_shared_file, "rshared", "groups")
files <- glue("{file_stub}.{extensions}")


map_dfr(files, read_rarefaction_files) %>% 
  drop_na() %>%
  group_by(n_seqs, metric) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            n_samples = n(), .groups = "drop") %>%
  write_tsv(output_file)

temp_files <- list.files(path = dirname(temp_shared_file),
                         pattern = str_replace(basename(temp_shared_file),
                                               ".rshared", ""),
                         full.names = TRUE)

unlink(temp_files)

#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(SRS)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # shared_file <- "data/soil/data.otu.1.rshared"


srs_shared_file <- str_replace(shared_file, "(.)shared", "\\1_srs.shared")

shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  select(Group, starts_with("Otu"))

smallest_group <- shared %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(n = sum(value)) %>%
  summarize(min = min(n)) %>%
  pull(min)

srs_df <- shared %>%
  column_to_rownames("Group") %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin = smallest_group) %>%
  t()
  
colnames(srs_df) <- glue("Otu{ str_pad(1:ncol(srs_df),
                                       pad = 0,
                                       width = nchar(ncol(srs_df)))}"
                         )
  
srs_df %>%
  as_tibble(rownames = "Group") %>%
  mutate(label = 0.03, numOtus = ncol(.) - 1) %>%
  select(label, Group, numOtus, everything()) %>%
  write_tsv(srs_shared_file)

system(

  glue('mothur "#summary.single(shared={srs_shared_file},
                    calc=nseqs-sobs-shannon-simpson-chao-ace-npshannon-coverage,
                    subsample=F)"')

)

rabund_file_stub <- str_replace(srs_shared_file, "shared", "")
rabund_files <- glue("{rabund_file_stub}{rownames(srs_df)}.rabund")

file.remove(rabund_files)
file.remove(srs_shared_file)

mothur_output <- str_replace(srs_shared_file,
                             ".shared",
                             ".groups.summary")
  
alpha_output <- str_replace(srs_shared_file,
                            ".shared",
                            "_alpha")

file.rename(mothur_output, alpha_output)
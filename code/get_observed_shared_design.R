#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

shared_file <- "data/human/data.otu.shared"
remove_file <- "data/human/data.remove_accnos"

obs_design_file <- "data/human/data.odesign"
obs_shared_file <- "data/human/data.otu.oshared"

removal <- scan(remove_file, quiet = TRUE, what = character()) %>%
  str_split(., "-") %>%
  unlist()

# some sample names appear twice, but their diagnoses are the same. this was
# because some samples were sequenced twice
read_tsv("data/human/sra_info.tsv",
         col_type=cols(Sample_Name = col_character())) %>%
  mutate(treatment = if_else(
    diagnosis == "Cancer" | diagnosis == "adv Adenoma",
    "srn", "healthy"))          %>%
  count(Sample_Name, treatment) %>%
  select(Group = Sample_Name, treatment) %>%
  filter(! Group %in% removal) %>%
  write_tsv(obs_design_file)

rand_shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  filter(! Group %in% removal) %>%
  write_tsv(obs_shared_file)

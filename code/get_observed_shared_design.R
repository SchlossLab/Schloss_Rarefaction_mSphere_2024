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
         col_type = cols(Sample_Name = col_character())) %>% 
  mutate(srn = diagnosis == "Cancer" | diagnosis == "adv Adenoma",
         obese = BMI >= 30,
         smoker = as.logical(Smoke),
         nsaid = as.logical(NSAID),
         prev_crc = as.logical(Hx_Prev),
         prev_polyps = as.logical(Hx_of_Polyps),
         female = Gender == 'female'
         ) %>%
  select(Group = Sample_Name, srn, obese, female, smoker, nsaid, 
         prev_crc, prev_polyps) %>%
  filter(! Group %in% removal) %>%
  write_tsv(obs_design_file)

rand_shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  filter(! Group %in% removal) %>%
  write_tsv(obs_shared_file)

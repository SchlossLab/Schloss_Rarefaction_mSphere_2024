#!/usr/bin/env Rscript

# This script reads in the name of the directory holding the data (e.g.
# data/marine) and merges the information in the sra_info.tsv file (described
# in code/datasets_download.sh) with the actual sequence files that are
# available to generate the files file that is used by mothur

library(tidyverse)

input <- commandArgs(trailingOnly = TRUE)

path <- input[1] # e.g. path <- 'data/soil'


# Generate a tibble based on the R1 files found in path that contains the name
# of the R1 and R2 files as well as the SRR name (i.e. stub column). This step
# assumes that every sample has two reads. The cases with single read  files
# should have been removed in the code/datasts_download.sh step

fastq_gz <- tibble(read_1 = list.files(path = path, pattern = "*1.fastq.gz"),
                  read_2 = str_replace(read_1, "_1", "_2"),
                  stub = str_replace(read_1, "_1.*", ""))


# If the data were originally deposited into the SRA, the run files will start
# with SRR and the sample names will be stored in the "Sample_Name" column. If
# they were originally deposited into the ENA, then the run files will start
# with ERR and the sample names will be stored in the "Alias" column

sample_id <- ifelse(str_detect(fastq_gz[1, "stub"], "SRR"),
                    "Sample_Name",
                    "Alias")


# We want to read in the sra_info file, extract the relevant column that
# contains the sample name, and then merge that information with the fastq_gz
# tibble. Finally, we select out the columns with the sample names and the
# fastq.gz file names and output this as path/data.files. If sample_id has a
# hyphen, turn it into an underscore

read_tsv(file = paste0(path, "/sra_info.tsv"),
         col_types = cols(.default = col_character())) %>%
  select(all_of(sample_id), Run) %>%
  left_join(fastq_gz, ., by = c("stub" = "Run")) %>%
  rename(sample_id = all_of(sample_id)) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_")) %>%
  select(sample_id, read_1, read_2) %>%
  write_tsv(paste0(path, "/data.files"), col_names = F, quote = "none")

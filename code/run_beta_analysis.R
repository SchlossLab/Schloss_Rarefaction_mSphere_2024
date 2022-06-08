#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(vegan)

run_amova <- function(distance_file, design_file) {

  design <- read_tsv(design_file,
                     col_types = cols(.default = col_character()))

  distance_design <- read_tsv(distance_file,
                              col_types = cols(Group = col_character(),
                                               .default = col_double()
                                               )
                              ) %>%
    inner_join(design, ., by = "Group")

  w_na <- nrow(distance_design)
  wo_na <- distance_design %>%
    drop_na() %>%
    nrow()

  p_value <- NA_real_

  if (w_na == wo_na) {
    distance_matrix <- distance_design %>%
      select(Group, distance_design$Group) %>%
      column_to_rownames("Group") %>%
      as.dist()

    test <- vegan::adonis2(distance_matrix~distance_design$treatment)

    p_value <- test$`Pr(>F)`[1]
  }

  return(p_value)

}

input <- commandArgs(trailingOnly = TRUE)

dataset <- input[1] # dataset <- "data/human"
output_file <- input[2] # output_file <- "data/human/data.e_eamova"

model <- str_replace(output_file, ".*\\.(.)_.amova", "\\1")
design <- str_replace(output_file, ".*\\.._(.)amova", "\\1")

distance_files <- tibble(
    distance_file = list.files(path = dataset,
                               pattern = glue(".*\\.{model}_.*.dist"),
                               full.names = TRUE)) %>%
  mutate(seed = str_replace(distance_file, ".*\\.(\\d+)\\..*", "\\1"))

design_files <- tibble(
    design_file = list.files(path = dataset,
                             pattern = glue(".*.{design}design"),
                             full.names = TRUE)) %>%
  mutate(seed = str_replace(design_file, ".*\\.(\\d+)\\..*", "\\1"))

inner_join(distance_files, design_files, by = "seed") %>%
  rowwise() %>%
  mutate(p_value = run_amova(distance_file, design_file)) %>%
  ungroup() %>%
  separate(distance_file,
           into = c("model", "beta_process", "beta_calculator"),
           sep = "_") %>%
  mutate(model = str_replace(model, ".*\\.", ""),
         beta_calculator = str_replace(beta_calculator, ".dist", "")) %>%
  select(-design_file) %>%
  group_by(model, beta_process, beta_calculator) %>%
  summarize(significant = sum(p_value < 0.05) / n(), .groups = "drop") %>%
  write_tsv(output_file)
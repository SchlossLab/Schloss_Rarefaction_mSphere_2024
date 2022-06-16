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

  result <- tibble(p.value = NA_real_, r_squared = NA_real_)

  if (w_na == wo_na) {
    distance_matrix <- distance_design %>%
      select(Group, distance_design$Group) %>%
      column_to_rownames("Group") %>%
      as.dist()

    test <- vegan::adonis2(distance_matrix~distance_design$treatment)

    result <- tibble(p.value = test$`Pr(>F)`[1],
                     r_squared = test$R2[1])
  }

  return(result)

}

dataset <- "data/human"
output_file <- "data/human/data.o_oamova"

model <- str_replace(output_file, ".*\\.(.)_.amova", "\\1")
design <- str_replace(output_file, ".*\\.._(.)amova", "\\1")

distance_design_files <- tibble(
    distance_file = list.files(path = dataset,
                               pattern = glue(".*\\.{model}_.*.dist"),
                               full.names = TRUE),
    design_file = "data/human/data.odesign")


distance_design_files %>%
  rowwise() %>%
  mutate(result = run_amova(distance_file, design_file)) %>%
  unnest(result) %>%
  ungroup() %>%
  separate(distance_file,
           into = c("model", "beta_process", "beta_calculator"),
           sep = "_") %>%
  mutate(beta_calculator = str_replace(beta_calculator, ".dist", "")) %>%
  select(-design_file, -model) %>%
  write_tsv(output_file)
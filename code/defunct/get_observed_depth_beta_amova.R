#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(vegan)

run_amova <- function(distance_file, treatment) {

  design <- read_tsv("data/human/data.odesign",
                     col_types = cols(.default = col_character())) %>%
    select(Group, "treatment" = all_of(treatment)) %>%
    distinct()

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

  result <- tibble(p.value = NA_real_,
                   r_squared = NA_real_,
                   n_FALSE = NA_real_,
                   n_TRUE = NA_real_
                   )

  if (w_na == wo_na) {
    distance_matrix <- distance_design %>%
      select(Group, distance_design$Group) %>%
      column_to_rownames("Group") %>%
      as.dist()

    test <- vegan::adonis2(distance_matrix~distance_design$treatment)

    tf_count <- distance_design %>% count(treatment) %>%
      pivot_wider(names_from = treatment, values_from = "n")

    result <- tibble(p.value = test$`Pr(>F)`[1],
                     r_squared = test$R2[1],
                     n_FALSE = tf_count$`FALSE`,
                     n_TRUE = tf_count$`TRUE`)
  }

  return(result)

}


design_file <- read_tsv("data/human/data.odesign")
treatments <- colnames(design_file)[-1]

output_file <- "data/human/data.otu.obs_depth_beta_amova.summary"

distance_design <- expand_grid(
    distance_file = list.files(path = "data/human",
                               pattern = "data.otu.obs_\\d*_\\w*.dist",
                               full.names = TRUE),
    treatment = treatments
    )


distance_design %>%
  rowwise() %>%
  mutate(result = run_amova(distance_file, treatment)) %>%
  unnest(result) %>%
  ungroup() %>%
  separate(distance_file,
           into = c("dummy", "nseqs", "beta_calculator"),
           sep = "_") %>%
  mutate(beta_calculator = str_replace(beta_calculator, ".dist", "")) %>%
  select(-dummy) %>%
  write_tsv(output_file)
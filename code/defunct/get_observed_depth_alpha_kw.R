#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(effectsize)
library(broom)

get_rarefy_results <- function(rfile, design) {

  alpha_data <- read_tsv(rfile, col_types = cols(group = col_character(),
                                                 method = col_character(),
                                                 .default = col_double()))

  inner_join(alpha_data,
             design,
             by = "group") %>% 
    filter(method == "ave") %>%
    select(Group = group, treatment, sobs, shannon, simpson, coverage) %>%
    pivot_longer(cols = -c(Group, treatment), names_to = "metric") %>%
    nest(data = -metric) %>%
    mutate(kw = map(data, ~tidy(kruskal.test(.x$value, g = .x$treatment))),
           biserial = map(data, ~rank_biserial(.x$value~as.factor(.x$treatment),
                                                    ci = NULL)[[1]][[1]]),
           median = map(data,
                        ~group_by(.x, treatment) %>%
                          summarize(median = median(value),
                                    n = n())
                          )) %>%
    unnest(c(kw, median, biserial)) %>%
    select(metric, p.value, biserial, treatment, median, n) %>%
    pivot_wider(names_from = "treatment", values_from = c("median", "n"))

}


run_kw_tests <- function(alpha_file, treatment) {

  nseqs <- as.numeric(str_replace(alpha_file, ".*_(\\d*)_.*", "\\1"))

  design_file <- "data/human/data.odesign"

  design <- read_tsv(design_file,
                     col_types = cols(.default = col_character())) %>%
    select(group = Group, treatment = all_of(treatment))

  get_rarefy_results(alpha_file, design) %>%
    mutate("treatment" = treatment) %>%
    mutate(nseqs = nseqs, .after = "metric")

}

#########

design_file <- "data/human/data.odesign"
output_file <- "data/human/data.otu.obs_depth_alpha_kw.sumary"

summary_files <- list.files(path = "data/human",
                            pattern = "data.otu.obs_\\d*_alpha.summary",
                            full.names = TRUE)

metadata <- c("srn", "obese", "female", "smoker", "nsaid", "prev_crc",
              "prev_polyps")


combinations <- expand_grid(metadata, summary_files)

map2_dfr(combinations$summary_files,
         combinations$metadata,
         ~run_kw_tests(.x, .y)) %>%
  mutate(model = "observed") %>%
  write_tsv(output_file)

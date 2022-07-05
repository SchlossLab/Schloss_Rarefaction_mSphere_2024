#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)

get_rarefy_results <- function(rfile, design) {

  alpha_data <- read_tsv(rfile, col_types = cols(group = col_character(),
                                                 method = col_character(),
                                                 .default = col_double()))

  inner_join(alpha_data,
             design,
             by = "group") %>%
    filter(method == "ave") %>%
    select(Group = group, treatment, nseqs, sobs, shannon,
          simpson, chao, ace, npshannon, coverage) %>%
    pivot_longer(cols = -c(Group, treatment), names_to = "metric") %>%
    nest(data = -metric) %>%
    mutate(kw = map(data, ~tidy(kruskal.test(.x$value, g = .x$treatment)))) %>%
    unnest(kw) %>%
    mutate(median = map(data,
                        ~group_by(.x, treatment) %>%
                          summarize(median = median(value))
                          )) %>%
    unnest(median) %>%
    select(metric, p.value, statistic, treatment, median) %>%
    pivot_wider(names_from = "treatment", values_from = "median") %>%
    mutate(method = "rarefy")

}

get_breakaway_results <- function(bafile, design) {

  alpha_data <- read_tsv(bafile, col_types = cols(Group = col_character(),
                                                 ba_model = col_character(),
                                                 p_model = col_character(),
                                                 .default = col_double()))

  inner_join(alpha_data,
             design,
             by = c("Group" = "group")) %>%
    select(Group, treatment, ba_est, ba_model, p_est, p_model) %>%
    pivot_longer(cols = -c(Group, treatment), names_to = c("set", ".value"),
                 names_pattern = "(.*)_(.*)") %>%
    nest(data = -set) %>%
    mutate(kw = map(data, ~tidy(kruskal.test(.x$est, g = .x$treatment)))) %>%
    unnest(kw) %>%
    mutate(median = map(data,
                        ~group_by(.x, treatment) %>%
                          summarize(median = median(est))
                          )) %>%
    unnest(median) %>%
    select(metric = set, p.value, statistic, treatment, median) %>%
    pivot_wider(names_from = "treatment", values_from = "median") %>%
    mutate(method = "breakaway",
          metric = recode(metric, "ba" = "default", "p" = "poisson"))

}

get_raw_results <- function(rawfile, design) {

  alpha_data <- read_tsv(rawfile, col_types = cols(group = col_character(),
                                                 .default = col_double()))

  inner_join(alpha_data,
             design,
             by = "group") %>%
    select(group, treatment, nseqs, sobs, shannon,
          simpson, chao, ace, npshannon, coverage) %>%
    pivot_longer(cols = -c(group, treatment), names_to = "metric") %>%
    nest(data = -metric) %>%
    mutate(kw = map(data, ~tidy(kruskal.test(.x$value, g = .x$treatment)))) %>%
    unnest(kw) %>%
    mutate(median = map(data,
                        ~group_by(.x, treatment) %>%
                          summarize(median = median(value))
                          )) %>%
    unnest(median) %>%
    select(metric, p.value, statistic, treatment, median) %>%
    pivot_wider(names_from = "treatment", values_from = "median") %>%
    mutate(method = "raw")

}

get_srs_results <- function(srsfile, design) {

  alpha_data <- read_tsv(srsfile, col_types = cols(group = col_character(),
                                                 .default = col_double()))

  inner_join(alpha_data,
             design,
             by = "group") %>%
    select(Group = group, treatment, nseqs, sobs, shannon,
          simpson, chao, ace, npshannon, coverage) %>%
    pivot_longer(cols = -c(Group, treatment), names_to = "metric") %>%
    nest(data = -metric) %>%
    mutate(kw = map(data, ~tidy(kruskal.test(.x$value, g = .x$treatment)))) %>%
    unnest(kw) %>%
    mutate(median = map(data,
                        ~group_by(.x, treatment) %>%
                          summarize(median = median(value))
                          )) %>%
    unnest(median) %>%
    select(metric, p.value, statistic, treatment, median) %>%
    pivot_wider(names_from = "treatment", values_from = "median") %>%
    mutate(method = "srs")

}

run_kw_tests <- function(method, treatment) {

  f <- glue("get_{method}_results")

  alpha_file <- glue("data/human/data.otu.o_{method}_alpha")
  design_file <- "data/human/data.odesign"

  design <- read_tsv(design_file,
                     col_types = cols(.default = col_character())) %>%
    select(group = Group, treatment = treatment)

  do.call(f, list(alpha_file, design)) %>%
    mutate("treatment" = treatment)

}

#########

design_file <- "data/human/data.odesign"
output_file <- "data/human/data.o_oalpha_kw"


metadata <- c("srn", "obese", "female", "smoker", "nsaid", "prev_crc",
              "prev_polyps")
methods <- c("rarefy", "breakaway", "raw", "srs")


combinations <- expand_grid(metadata, methods)

map2_dfr(combinations$methods, combinations$metadata, ~run_kw_tests(.x, .y)) %>%
  mutate(model = "observed") %>%
  write_tsv(output_file)


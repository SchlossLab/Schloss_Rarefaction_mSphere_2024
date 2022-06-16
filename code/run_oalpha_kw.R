#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)

get_rarefy_results <- function(rfile, dfile) {

  inner_join(read_tsv(rfile, show_col_types = FALSE),
             read_tsv(dfile, show_col_types = FALSE),
             by = c("group" = "Group")) %>%
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

get_breakaway_results <- function(bafile, dfile) {

  inner_join(read_tsv(bafile, show_col_types = FALSE),
             read_tsv(dfile, show_col_types = FALSE),
             by = c("Group" = "Group")) %>%
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

get_raw_results <- function(rawfile, dfile) {

  inner_join(read_tsv(rawfile, show_col_types = FALSE),
             read_tsv(dfile, show_col_types = FALSE),
             by = c("group" = "Group")) %>%
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
    mutate(method = "raw")

}

get_srs_results <- function(srsfile, dfile) {

  inner_join(read_tsv(srsfile, show_col_types = FALSE),
             read_tsv(dfile, show_col_types = FALSE),
             by = c("group" = "Group")) %>%
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

run_kw_tests <- function(method) {

  f <- glue("get_{method}_results")

  alpha_file <- glue("data/human/data.otu.o_{method}_alpha")
  design_file <- "data/human/data.odesign"

  do.call(f, list(alpha_file, design_file))

}

#########

design_file <- "data/human/data.odesign"

output_file <- "data/human/data.o_oalpha_kw"

methods <- c("rarefy", "breakaway", "raw", "srs")

map_dfr(methods, ~run_kw_tests(.x)) %>%
  mutate(model = "observed") %>%
  write_tsv(output_file)


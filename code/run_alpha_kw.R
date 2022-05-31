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
    select(metric, p.value) %>%
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
    select(metric = set, p.value) %>%
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
    select(metric, p.value) %>%
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
    select(metric, p.value) %>%
    mutate(method = "srs")

}

run_kw_tests <- function(dataset, method, model) {

  f <- glue("get_{method}_results")

  tibble(
      alpha_file = list.files(path = dataset,
                pattern = glue(".*{model}_{method}_alpha"),
                include.dirs = TRUE, full.names = TRUE)) %>%
      mutate(seed = str_replace(alpha_file, ".*\\.(\\d+)\\..*", "\\1")) %>%
      inner_join(., design_files, by = "seed") %>%
      nest(data = -seed) %>%
      mutate(results = map(data,
                          ~do.call(f, list(.x$alpha_file, .x$design_file)))) %>%
      select(results) %>%
      unnest(results) %>%
      group_by(metric, method) %>%
      summarize(frac = sum(p.value < 0.05) / n(), .groups = "drop")
}

#########

input <- commandArgs(trailingOnly = TRUE)

dataset <- input[1] # dataset <- "data/mice"
output_file <- input[2] # output_file <- "data/mice/data.r_ralpha_kw"

model <- str_replace(output_file, ".*\\.(.)_.alpha_kw", "\\1")
design <- str_replace(output_file, ".*\\.._(.)alpha_kw", "\\1")

design_files <- tibble(
    design_file = list.files(path = dataset,
                             pattern = glue(".*.{design}design"),
                             full.names = TRUE)) %>%
  mutate(seed = str_replace(design_file, ".*\\.(\\d+)\\..*", "\\1"))


methods <- c("rarefy", "breakaway", "raw", "srs")

map_dfr(methods, ~run_kw_tests(dataset, .x, model)) %>%
  mutate(model = model) %>%
  write_tsv(output_file)


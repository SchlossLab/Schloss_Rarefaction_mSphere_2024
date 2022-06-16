#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)

# data/human/data.otu.1.r_deseq2_bray.dist
# data/human/data.otu.1.r_deseq2_euclidean.dist
# data/human/data.otu.1.r_deseq2_jaccard.dist
# data/human/data.otu.1.r_metagenomeseq_bray.dist
# data/human/data.otu.1.r_metagenomeseq_euclidean.dist
# data/human/data.otu.1.r_metagenomeseq_jaccard.dist
# data/human/data.otu.1.r_nclr_bray.dist
# data/human/data.otu.1.r_nclr_euclidean.dist
# data/human/data.otu.1.r_nclr_jaccard.dist
# data/human/data.otu.1.r_oclr_bray.dist
# data/human/data.otu.1.r_oclr_euclidean.dist
# data/human/data.otu.1.r_oclr_jaccard.dist
# data/human/data.otu.1.r_rare_bray.dist
# data/human/data.otu.1.r_rare_euclidean.dist
# data/human/data.otu.1.r_rare_jaccard.dist
# data/human/data.otu.1.r_raw_bray.dist
# data/human/data.otu.1.r_raw_euclidean.dist
# data/human/data.otu.1.r_raw_jaccard.dist
# data/human/data.otu.1.r_rclr_bray.dist
# data/human/data.otu.1.r_rclr_euclidean.dist
# data/human/data.otu.1.r_rclr_jaccard.dist
# data/human/data.otu.1.r_relabund_bray.dist
# data/human/data.otu.1.r_relabund_euclidean.dist
# data/human/data.otu.1.r_relabund_jaccard.dist
# data/human/data.otu.1.r_srs_bray.dist
# data/human/data.otu.1.r_srs_euclidean.dist
# data/human/data.otu.1.r_srs_jaccard.dist
# data/human/data.otu.1.r_zclr_bray.dist
# data/human/data.otu.1.r_zclr_euclidean.dist
# data/human/data.otu.1.r_zclr_jaccard.dist

input <- commandArgs(trailingOnly = TRUE)
dataset <- input[1] # e.g. dataset <- "human"

output_file <- glue("data/{dataset}/random_beta_correlation.tsv")

n_seqs <- read_tsv(glue("data/{dataset}/data.group_count"),
                   col_names = c("group", "nseqs"),
                   col_type = cols(group = col_character()))


dist_files <- list.files(path = glue("data/{dataset}"),
                         pattern = "data\\.otu\\.\\d*\\.r_.*dist",
                         include.dirs = TRUE,
                         full.names = TRUE)

get_correlation <- function(x) {

  dataset <- str_replace(x, "data/(.*)/data.otu.*", "\\1")
  seed <- str_replace(x, ".*data.otu.(\\d*)\\..*", "\\1")
  process <- str_replace(x, ".*data.otu.\\d*\\.r_(.*)_.*", "\\1")
  calculator <- str_replace(x, ".*data.otu.\\d*\\.r_.*_(.*)\\.dist", "\\1")

  input <- read_tsv(x, col_types = cols(Group = col_character(),
                               .default = col_double())) %>%
    pivot_longer(-Group, names_to = "cols", values_to = "dist") %>%
    select(rows = Group, everything()) %>%
    filter(rows < cols) %>%
    inner_join(., n_seqs, by = c("rows" = "group")) %>%
    inner_join(., n_seqs, by = c("cols" = "group")) %>%
    mutate(diff = abs(nseqs.x - nseqs.y),
          dataset = dataset,
          seed = seed,
          process = process,
          calculator = calculator)

  if (sum(is.na(input$dist)) == 0) {

    input %>%
      nest(data = -c(dataset, seed, process, calculator)) %>%
      mutate(test = map(data, ~cor.test(.x$diff, .x$dist,
                                        method = "spearman",
                                        exact = FALSE) %>%
                        tidy)) %>%
      unnest(test) %>%
      select(dataset, seed, process, calculator, estimate, p.value)

  } else {

    tibble(dataset = dataset, seed = seed,
           process = process, calculator = calculator,
          estimate = NA_real_, p.value = NA_real_)
  }
}

map_dfr(dist_files, get_correlation) %>%
  drop_na() %>%
  group_by(dataset, process, calculator) %>%
  summarize(median = median(estimate),
            iqr = IQR(estimate),
            frac_sig = sum(p.value < 0.05) / n(),
            .groups = "drop") %>%
  write_tsv(output_file)
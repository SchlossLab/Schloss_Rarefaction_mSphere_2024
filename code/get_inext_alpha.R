#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(iNEXT)

set.seed(19760620)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/mice/data.otu.1.rshared"

inext_output_file <- str_replace(shared_file,
                                "shared",
                                "_inext_alpha")

shared_df <- fread(shared_file,
      colClasses = c(Group = "character")) %>%
  select(-label, -numOtus) %>%
  melt(measure = patterns("^Otu"),
      variable.name = "otu", value.name = "count") %>%
  dcast(otu ~ Group, value.var = "count") %>%
  column_to_rownames("otu")

nseqs <- apply(shared_df, 2, sum) %>%
  as_tibble(rownames = "Group")

asy_richness <- ChaoRichness(shared_df) %>%
  as_tibble(rownames = "Group") %>%
  select(Group, chao_sobs = Estimator)

asy_shannon <- ChaoShannon(shared_df, B = 1) %>%
  as_tibble(rownames = "Group") %>%
  select(Group, chao_shannon = Estimator)

asy_simpson <- ChaoSimpson(shared_df, B = 1, transform = TRUE) %>%
  as_tibble(rownames = "Group") %>%
  select(Group, chao_invsimpson = Estimator)

asy <- inner_join(asy_richness, asy_shannon, by = "Group") %>%
  inner_join(., asy_simpson, by = "Group")

coverage_qd <- estimateD(shared_df,
                         q = c(0, 1, 2),
                         datatype = "abundance",
                         base = "coverage",
                         nboot = 0) %>%
  mutate(Diversity = recode(Order.q,
                            "0" = "Species richness",
                            "1" = "Shannon diversity",
                            "2" = "Simpson diversity")) %>%
  select(Assemblage, coverage_method = Method, Diversity, Coverage = qD) %>%
  pivot_wider(names_from = Diversity, values_from = Coverage) %>%
  select(Group = Assemblage,
        coverage_method,
        coverage_sobs = `Species richness`,
        coverage_shannon = `Shannon diversity`,
        coverage_invsimpson = `Simpson diversity`)

size_qd <- estimateD(shared_df,
                     q = c(0, 1, 2),
                     datatype = "abundance",
                     base = "size",
                     nboot = 0) %>%
  mutate(Diversity = recode(Order.q,
                            "0" = "Species richness",
                            "1" = "Shannon diversity",
                            "2" = "Simpson diversity")) %>%
  select(Assemblage, size_method = Method, Diversity, Size = qD) %>%
  pivot_wider(names_from = Diversity, values_from = Size) %>%
  select(Group = Assemblage,
        size_method,
        size_sobs = `Species richness`,
        size_shannon = `Shannon diversity`,
        size_invsimpson = `Simpson diversity`)

inner_join(nseqs, asy, by = "Group") %>%
  inner_join(., coverage_qd, by = "Group") %>%
  inner_join(., size_qd, by = "Group") %>%
  rename(group = Group) %>%
  write_tsv(inext_output_file)

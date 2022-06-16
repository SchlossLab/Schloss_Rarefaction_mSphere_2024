#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

input <- commandArgs(trailingOnly = TRUE)
shared_file <- input[1] # e.g. shared_file <- "data/human/data.otu.shared"
remove_file <- input[2] # e.g. remove_file <- "data/human/data.remove_accnos"
seed <- as.numeric(input[3]) # e.g. seed <- 87

# could also do 0.10 and 0.05 (i.e. flipping the values below)
remove_factor <- 0.03 # this is the fraction of OTUs to remove

# e.g. s_shared_file <- "data/human/data.otu.87.cshared"
c_shared_file <- str_replace(shared_file,
                             "shared",
                             glue("{seed}.cshared"))

c_design_file <- glue("{dirname(shared_file)}/data.{seed}.cdesign")


#	set the random number generator seed so multiple runs generate the same 
#	randomization
set.seed(seed)

removal <- scan(remove_file, quiet = TRUE, what = character()) %>%
  str_split(., "-") %>%
  unlist()


tidy_shared <- fread(shared_file, colClasses = c(Group = "character")) %>%
  pivot_longer(cols = -c(label, Group, numOtus),
               names_to = "otu",
               values_to = "n_seqs") %>%
  filter(! Group %in% removal) %>%
  select(Group, otu, n_seqs)
  
group_counts <- tidy_shared %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(n_seqs)) %>%
  mutate(treatment = sample(rep(c("A", "B"), length.out = nrow(.))))
  
seq_dist_original <- tidy_shared %>%
  group_by(otu) %>%
  summarize(n_seqs = sum(n_seqs)) %>%
  mutate(prob = n_seqs / sum(n_seqs)) %>%
  select(otu, prob)
  
n_otus <- nrow(seq_dist_original)

rand_perturbed <- c(rep(T, round(remove_factor * n_otus)),
                   rep(F, round(n_otus - round(remove_factor * n_otus)))) %>%
  sample(., n_otus)

seq_dist_perturbed <- seq_dist_original %>%
  mutate(change = rand_perturbed,
         prob = if_else(change,
                        0, #remove OTU if rand_perturbed is TRUE
                        prob)
  ) %>%
  mutate(pert_total = sum(.data$prob[change]),
         no_pert_total = sum(.data$prob[!change]),
         deduct = (1 - pert_total) / no_pert_total,
         prob = if_else(change, prob, prob * deduct)) %>%
  select(otu, prob)
  
seq_counts_original <- group_counts %>%
  filter(treatment == "A") %>%
  mutate(otus = map(n_seqs,
                    ~tibble(otu = sample(seq_dist_original$otu,
                                         size = .x,
                                         prob = seq_dist_original$prob,
                                         replace = TRUE)) %>%
                    count(otu))) %>%
  unnest(otus) %>%
  select(Group, treatment, otu, n)
  
seq_counts_perturbed <- group_counts %>%
  filter(treatment == "B") %>%
  mutate(otus = map(n_seqs,
                    ~tibble(otu = sample(seq_dist_perturbed$otu,
                                         size = .x,
                                         prob = seq_dist_perturbed$prob,
                                         replace = TRUE)) %>%
                    count(otu))) %>%
  unnest(otus) %>%
  select(Group, treatment, otu, n)
  
seq_counts_composite <- bind_rows(seq_counts_original, seq_counts_perturbed) %>%
  pivot_wider(names_from = otu, values_from = n, values_fill = 0)


seq_counts_composite %>%
  mutate(label = "0.03", numOtus = ncol(.) - 2) %>%
  select(label, Group, numOtus, starts_with("Otu")) %>%
  write_tsv(c_shared_file)

seq_counts_composite %>%
  select(Group, treatment) %>%
  write_tsv(c_design_file)
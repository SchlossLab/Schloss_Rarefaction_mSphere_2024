#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)
library(ggtext)

raw_sobs <- "data/human/data.otu.1.r_raw_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
    select(group, nseqs,
           raw = sobs)

rare_sobs <- "data/human/data.otu.1.r_rarefy_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
    filter(method == "ave") %>%
    select(group,
           rare = sobs)

srs_sobs <- "data/human/data.otu.1.r_srs_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
    select(group,
           srs = sobs)

ba_sobs <- "data/human/data.otu.1.r_breakaway_alpha" %>%
  read_tsv(col_types = cols(Group = col_character())) %>%
    select(group = Group,
           ba_default = ba_est,
           ba_poisson = p_est)



tidy_composite <- inner_join(raw_sobs, rare_sobs, by = "group") %>%
  inner_join(., srs_sobs, by = "group") %>%
  inner_join(., ba_sobs, by = "group") %>%
  select(-ba_poisson) %>%
  pivot_longer(cols = -c(group, nseqs),
               names_to = "metric", values_to = "value")
  
correlations <- tidy_composite %>%
  nest(data = -c(metric)) %>%
  mutate(test = map(data, ~cor.test(.x$nseqs, .x$value, method = "spearman",
                                    exact = FALSE) %>%
                    tidy())) %>%
  unnest(test) %>%
  select(metric, estimate, p.value)

get_rho <- function(m, df = correlations) {

  df %>%
    filter(metric == m) %>%
    mutate(estimate = format(round(estimate, digits = 2), nsmall = 2),
           p = case_when(p.value < 0.001 ~ glue("P < 0.001"),
                         p.value < 0.010 ~ glue("P = {format(round(p.value, digits = 3), nsamll = 3)}"),
                         TRUE ~ glue("P = {format(round(p.value, digits = 2), nsamll = 2)}")),
           pretty = glue("{estimate}; {p}")
           ) %>%
    pull(pretty)

}

pretty_labels <- c(
  raw = glue("No rarefaction\n(\u03C1 = {get_rho('raw')})"),
  rare = glue("Rarefaction\n(\u03C1 = {get_rho('rare')})"),
  srs = glue("SRS normalization\n(\u03C1 = {get_rho('srs')})"),
  ba_default = glue("Breakaway estimate\n(\u03C1 = {get_rho('ba_default')})") 
)

aplot <- tidy_composite %>%
  mutate(metric = factor(metric,
                         levels = c("rare", "raw", "srs", "ba_default"))) %>%
  filter(!(metric == "ba_default" & value > 4000)) %>%
  ggplot(aes(x = nseqs, y = value)) +
  geom_point() +
  # geom_smooth(se = FALSE, method = "lm")  +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller(metric = pretty_labels)) +
  scale_x_log10(breaks = c(1e4, 3e4, 1e5, 3e5),
                labels = c("1x10^4", "3x10^4", "1x10^5", "3x10^5"),
                           limits = c(1e4, NA)) +
  labs(x = "Number of sequences per sample",
       y = "Richness") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_markdown())
    
ggsave("figures/example_alpha_cor.tiff", aplot,
       width = 6, height = 6,
       compression = "lzw+p")

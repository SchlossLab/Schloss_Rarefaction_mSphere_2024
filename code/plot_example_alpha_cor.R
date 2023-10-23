#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)
library(ggtext)

raw <- "data/human/data.otu.1.r_raw_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
  select(group, nseqs,
        raw = sobs,
        ace)

rare <- "data/human/data.otu.1.r_rarefy_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
  filter(method == "ave") %>%
  select(group,
          rare = sobs)

srs <- "data/human/data.otu.1.r_srs_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
    select(group,
           srs = sobs)

ba <- "data/human/data.otu.1.r_breakaway_alpha" %>%
  read_tsv(col_types = cols(Group = col_character())) %>%
  select(group = Group,
        ba_default = ba_est,
        ba_poisson = p_est)

inext <- "data/human/data.otu.1.r_inext_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
  select(group,
        chao = chao_sobs,
        size = size_sobs,
        coverage = coverage_sobs
        )

inext_methods <- "data/human/data.otu.1.r_inext_alpha" %>%
  read_tsv(col_types = cols(group = col_character())) %>%
  select(group,
        coverage = coverage_method,
        size = size_method) %>%
  pivot_longer(-group, names_to = "metric", values_to = "method") %>%
  mutate(method = tolower(method))




tidy_composite <- inner_join(raw, rare, by = "group") %>%
  inner_join(., srs, by = "group") %>%
  inner_join(., ba, by = "group") %>%
  inner_join(., inext, by = "group") %>%
  pivot_longer(cols = -c(group, nseqs),
               names_to = "metric", values_to = "value") %>%
  left_join(., inext_methods, by = c("group", "metric")) %>%
  mutate(method = replace_na(method, "rarefaction"),
        metric = factor(metric, levels = c("rare", "size", "coverage",
                                            "raw", "ace", "chao",
                                            "srs", "ba_default",
                                            "ba_poisson")))


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
  ace = glue("ACE estimate\n(\u03C1 = {get_rho('ace')})"),
  chao = glue("Chao1 estimate\n(\u03C1 = {get_rho('chao')})"),
  size = glue("iNEXT estimate (Size)\n(\u03C1 = {get_rho('size')})"),
  coverage = glue("iNEXT estimate (Coverage)\n(\u03C1 = {get_rho('coverage')})"),
  ba_default = glue("Breakaway est. (Default)\n(\u03C1 = {get_rho('ba_default')})"),
  ba_poisson = glue("Breakaway est. (Poisson)\n(\u03C1 = {get_rho('ba_poisson')})") 
)

aplot <- tidy_composite %>%
  filter(!((metric == "ba_default" | metric == "ba_poisson") &
                value > 4000)) %>%
  ggplot(aes(x = nseqs, y = value, color = method)) +
  geom_point(show.legend = FALSE, size = 0.5) +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller(metric = pretty_labels)) +
  scale_x_log10(breaks = c(1e4, 3e4, 1e5, 3e5),
                labels = c("1x10<sup>4</sup>", "3x10<sup>4</sup>",
                           "1x10<sup>5</sup>", "3x10<sup>5</sup>"),
                limits = c(1e4, NA)) +
  scale_color_manual(breaks = c("extrapolation", "rarefaction"),
                    values = c("gray", "black")) +
  labs(x = "Number of sequences per sample",
       y = "Richness") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_markdown())

ggsave("figures/example_alpha_cor.tiff", aplot,
       width = 7, height = 6,
       compression = "lzw+p")

#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(broom)
library(ggtext)
library(patchwork)

dist_files <- c(
    "data/human/data.otu.1.r_deseq2_euclidean.dist",
    "data/human/data.otu.1.r_nclr_euclidean.dist",
    "data/human/data.otu.1.r_oclr_euclidean.dist",
    "data/human/data.otu.1.r_rclr_euclidean.dist",
    "data/human/data.otu.1.r_zclr_euclidean.dist",

    "data/human/data.otu.1.r_rare_bray.dist",
    "data/human/data.otu.1.r_rare_jaccard.dist",

    "data/human/data.otu.1.r_raw_bray.dist",
    "data/human/data.otu.1.r_raw_jaccard.dist",

    "data/human/data.otu.1.r_relabund_bray.dist",
    "data/human/data.otu.1.r_relabund_jaccard.dist",

    "data/human/data.otu.1.r_metagenomeseq_bray.dist",
    "data/human/data.otu.1.r_metagenomeseq_jaccard.dist",

    "data/human/data.otu.1.r_srs_bray.dist",
    "data/human/data.otu.1.r_srs_jaccard.dist"
  )




n_seqs <- read_tsv("data/human/data.group_count",
                   col_names = c("group", "nseqs"),
                   col_type = cols(group = col_character()))


get_dists <- function(x) {

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
          process = process,
          calculator = calculator)

}

distances <- map_dfr(dist_files, get_dists) %>%
  mutate(process = factor(process,
                          levels = c("raw", "rare", "relabund", "srs",
                                     "metagenomeseq", "nclr", "oclr", "rclr",
                                     "zclr", "deseq2")))

correlations <- distances %>%
  nest(data = -c(process, calculator)) %>%
  mutate(test = map(data, ~cor.test(.x$diff, .x$dist,
                                        method = "spearman",
                                        exact = FALSE) %>%
                        tidy)) %>%
  unnest(test) %>%
  select(process, calculator, estimate, p.value)
  
j_distances <- distances %>% filter(calculator == "jaccard")

b_distances <- distances %>% filter(calculator == "bray")

e_distances <- distances %>% filter(calculator == "euclidean")

get_rho <- function(p, m, df = correlations) {

  df %>%
    filter(calculator == m & process == p) %>%
    mutate(estimate = format(round(estimate, digits = 2), nsmall = 2),
           p = case_when(p.value < 0.001 ~ glue("P < 0.001"),
                         p.value < 0.010 ~ glue("P = {format(round(p.value, digits = 3))}"),
                         TRUE ~ glue("P = {format(round(p.value, digits = 2))}")),
           pretty = glue("{estimate}; {p}")
           ) %>%
    pull(pretty)

}


j_beta_labels <- c(
  raw = glue("Raw\n(\u03C1 = {get_rho('raw', 'jaccard')})"),
  rare = glue("Rarefied\n(\u03C1 = {get_rho('rare', 'jaccard')})"),
  relabund = glue("Rel. abundance\n(\u03C1 = {get_rho('relabund', 'jaccard')})"),
  srs = glue("SRS Normalized\n(\u03C1 = {get_rho('srs', 'jaccard')})"),
  metagenomeseq = glue("CSS Normalized\n(\u03C1 = {get_rho('metagenomeseq', 'jaccard')})")
)

b_beta_labels <- c(
  raw = glue("Raw\n(\u03C1 = {get_rho('raw', 'bray')})"),
  rare = glue("Rarefied\n(\u03C1 = {get_rho('rare', 'bray')})"),
  relabund = glue("Rel. abundance\n(\u03C1 = {get_rho('relabund', 'bray')})"),
  srs = glue("SRS Normalized\n(\u03C1 = {get_rho('srs', 'bray')})"),
  metagenomeseq = glue("CSS Normalized\n(\u03C1 = {get_rho('metagenomeseq', 'bray')})")
)

e_beta_labels <- c(
  nclr = glue("NCLR\n(\u03C1 = {get_rho('nclr', 'euclidean')})"),
  oclr = glue("OCLR\n(\u03C1 = {get_rho('oclr', 'euclidean')})"),
  rclr = glue("RCLR\n(\u03C1 = {get_rho('rclr', 'euclidean')})"),
  zclr = glue("ZCLR\n(\u03C1 = {get_rho('zclr', 'euclidean')})"),
  deseq2 = glue("DeSeq2\n(\u03C1 = {get_rho('deseq2', 'euclidean')})")
)

j_plot <- j_distances %>%
# slice_sample(n = 1000) %>% #remove this
  ggplot(aes(x = diff, y = dist)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.25) +
  facet_grid(calculator~process, labeller = labeller(process = j_beta_labels)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous() +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Jaccard distances") +
  theme(
    axis.text.x = element_blank(),
    axis.line.y = element_line()
  )
  
b_plot <- b_distances %>%
#  slice_sample(n = 1000) %>% #remove this
  ggplot(aes(x = diff, y = dist)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.25) +
  facet_grid(calculator~process, labeller = labeller(process = b_beta_labels)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous() +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Bray-Curtis distances") +
  theme(
    axis.text.x = element_blank(),
    axis.line.y = element_line()
  )


e_plot <- e_distances %>%
 # slice_sample(n = 1000) %>% #remove this
  ggplot(aes(x = diff, y = dist)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.25) +
  facet_grid(calculator~process, labeller = labeller(process = e_beta_labels)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 1e5, 2e5, 3e5, 4e5),
                     labels = c(0, 1, 2, 3, 4)) +
  coord_cartesian(clip = "off") +
  labs(x = "Difference in the number of sequences between pairs of samples (x10^5)",
       y = "Euclidean distances") +
  theme(
    axis.line.y = element_line()
  )


beta_plot <- j_plot / b_plot / e_plot &
  theme(
    panel.background = element_blank(),
    axis.title.x = element_markdown(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(face = "bold", hjust = 0, size = 8)
  )
  
#ggsave("figures/example_beta_cor.pdf", beta_plot,  width = 7, height = 7, device = cairo_pdf)
ggsave("figures/example_beta_cor.tiff", beta_plot,  width = 7, height = 7, compression = "lzw+p")

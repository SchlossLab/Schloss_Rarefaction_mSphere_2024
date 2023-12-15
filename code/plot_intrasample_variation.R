#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)
library(ggtext)

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

datasets <- c("bioethanol", "human", "lake", "marine", "mice", "peromyscus",
              "rainforest", "rice", "seagrass", "sediment", "soil", "stream")

pretty_datasets <- tibble(
  plain = datasets,
  pretty = capitalize(plain)
)


get_samples_to_remove <- function(file) {

  scan(file, what = character(), quiet = TRUE) %>%
    as_tibble_col(column_name = "sample")

}

remove_files <- glue("data/{datasets}/data.remove_accnos")
names(remove_files) <- datasets

remove_samples <- map_dfr(remove_files, get_samples_to_remove, .id = "dataset")


nseqs_files <- glue("data/{datasets}/data.group_count")
names(nseqs_files) <- datasets

smallest_sample_size <-
  map_dfr(nseqs_files, read_tsv, .id = "dataset",
                 col_names = c("sample", "n_seqs"),
                 col_types = cols(sample = col_character(),
                                  n_seqs = col_double())) %>%
  anti_join(., remove_samples) %>%
  group_by(dataset) %>%
  summarize(
    n_seqs = 2000 * ceiling(min(n_seqs) / 2000))



alpha_summary_files <- glue("data/{datasets}/data.otu.alpha_depth.summary")
names(alpha_summary_files) <- datasets

alpha_composite <- map_dfr(alpha_summary_files, read_tsv, .id = "dataset") %>%
  filter(method == "shannon" | method == "sobs") %>%
  group_by(dataset, method, n_seqs) %>%
  summarize(mean = mean(mean),
            cov = 100 * mean(sd / mean),
            sd = mean(sd),
            n_samples = mean(n_samples),
            .groups = "drop") %>%
  filter(n_samples >= 5) %>%
  pivot_longer(cols = c(mean, cov, sd),
               names_to = "statistic",
               values_to = "value") %>%
  mutate(dataset = factor(dataset, level = datasets),
         method = factor(method, levels = c("sobs", "shannon")),
         statistic = factor(statistic, levels = c("mean", "sd", "cov"))
         )




beta_summary_files <- glue("data/{datasets}/data.otu.beta_depth.summary")
names(beta_summary_files) <- datasets

beta_composite <- map_dfr(beta_summary_files, read_tsv, .id = "dataset") %>%
  filter(method == "bray") %>%
  group_by(dataset, method, n_seqs) %>%
  summarize(mean = mean(mean),
            cov = 100 * mean(sd / mean),
            sd = mean(sd),
            n_samples = mean(n_samples),
            .groups = "drop") %>%
  filter(n_samples >= 5) %>%
  pivot_longer(cols = c(mean, sd, cov),
               names_to = "statistic",
               values_to = "value") %>%
  mutate(dataset = factor(dataset, level = datasets),
         statistic = factor(statistic, levels = c("mean", "sd", "cov"))
         )


plot_method_statistic <- function(m, s, data = composite,
                                  smallest = smallest_sample_size) {

  filtered_line <- data %>%
    filter(method == m & statistic == s)

  filtered_point <- filtered_line %>%
    inner_join(., smallest, by = c("dataset")) %>%
    mutate(diff = abs(n_seqs.x - n_seqs.y)) %>%
    group_by(dataset) %>%
    slice_min(diff) %>%
    select(dataset, n_seqs = n_seqs.y, value) %>%
    ungroup()


  filtered_line %>%
    ggplot(aes(x = n_seqs, y = value,
               shape = dataset, group = dataset, color = dataset)) +
    geom_line() +
    geom_point(data = filtered_point, fill = "white", size = 2) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5),
                  labels = c("10<sup>3</sup>", "10<sup>4</sup>",
                              "10<sup>5</sup>")) +
    scale_color_manual(name = NULL,
                    values = rep(c("#1b9e77", "#d95f02", "#7570b3"), 4),
                    breaks = pretty_datasets$plain,
                    labels = pretty_datasets$pretty
                    ) +
    scale_shape_manual(name = NULL,
                    values = rep(c(15, 17, 19, 25), 3),
                    breaks = pretty_datasets$plain,
                    labels = pretty_datasets$pretty
                    ) +
    labs(x = NULL,
         y = NULL) +
    theme(
      strip.placement = "outside",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(lineend = "round"),
      axis.text.x = element_markdown()
    )
}

composite <- bind_rows(alpha_composite, beta_composite)

sobs_mean <- plot_method_statistic("sobs", "mean") +
  labs(y = "Mean value\nacross samples",
       title = "Richness") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     labels = c("0", "10,000", "20,000", "30,000"))

sobs_sd <- plot_method_statistic("sobs", "sd") +
  labs(y = "Standard\ndeviation") +
  scale_y_continuous(limits = c(0, NA))

sobs_cov <- plot_method_statistic("sobs", "cov") +
  labs(y = "Coefficient of\nvariation (%)") +
  scale_y_continuous(limits = c(0, 4))


shannon_mean <- plot_method_statistic("shannon", "mean") +
  labs(title = "Shannon")

shannon_sd <- plot_method_statistic("shannon", "sd") +
  scale_y_continuous(limits = c(0, NA))

shannon_cov <- plot_method_statistic("shannon", "cov") +
  scale_y_continuous(limits = c(0, NA))


bray_mean <- plot_method_statistic("bray", "mean") +
  labs(title = "Bray-Curtis")

bray_sd <- plot_method_statistic("bray", "sd") +
  scale_y_continuous(limits = c(0, NA))

bray_cov <- plot_method_statistic("bray", "cov") +
  scale_y_continuous(limits = c(0, NA))


patch <- sobs_mean + shannon_mean + bray_mean +
         sobs_sd + shannon_sd + bray_sd +
         sobs_cov + shannon_cov + bray_cov +
        plot_layout(ncol = 3, guides = "collect") +
        plot_annotation(caption = "Number of sequences sampled",
                        theme = theme(plot.caption =
                                      element_text(size = 12, hjust = 0.45))) &
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.key = element_blank())


ggsave("figures/intrasample_variation.tiff", patch,
       width = 8, height = 6,
       compression = "lzw+p")

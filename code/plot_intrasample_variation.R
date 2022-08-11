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
  filter(n_samples > 5) %>%
  pivot_longer(cols = c(mean, cov),
               names_to = "statistic",
               values_to = "value") %>%
  mutate(dataset = factor(dataset, level = datasets),
         method = factor(method, levels = c("sobs", "shannon")),
         statistic = factor(statistic, levels = c("mean", "cov"))
         )



datasets <- c("bioethanol", "human", "lake", "marine", "mice", "peromyscus",
              "rainforest", #"rice",
              "seagrass", "sediment", "soil"#,
              #"stream"
              )

beta_summary_files <- glue("data/{datasets}/data.otu.beta_depth.summary")
names(beta_summary_files) <- datasets

beta_composite <- map_dfr(beta_summary_files, read_tsv, .id = "dataset") %>%
  filter(method == "bray" | method == "jaccard") %>%
  group_by(dataset, method, n_seqs) %>%
  summarize(mean = mean(mean),
            cov = 100 * mean(sd / mean),
            sd = mean(sd),
            n_samples = mean(n_samples),
            .groups = "drop") %>% 
  filter(n_samples > 5) %>% 
  pivot_longer(cols = c(mean, cov),
               names_to = "statistic",
               values_to = "value") %>%
  mutate(dataset = factor(dataset, level = datasets),
         method = factor(method, levels = c("bray", "jaccard")),
         statistic = factor(statistic, levels = c("mean", "cov"))
         )


plot_method_statistic <- function(m, s, data = composite,
                                  smallest=smallest_sample_size) {

  filtered_line <- data %>%
    filter(method == m & statistic == s)

  filtered_point <- filtered_line %>%
    inner_join(., smallest, by = c("dataset", "n_seqs"))

  filtered_line %>%
    ggplot(aes(x = n_seqs, y = value,
               shape = dataset, group = dataset, color = dataset)) +
    geom_line() +
    geom_point(data = filtered_point, fill = "white") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5),
                  labels = c("10^3", "10^4", "10^5")) +
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

composite <- bind_rows(alpha_composite, beta_composite) %>%
  bind_rows(tibble(dataset = c(rep("rice", 8), rep("stream", 8)),
                   method = rep(rep(c("bray", "jaccard"), each = 4), 2),
                   n_seqs = 1000,
                   sd = 0,
                   n_samples = 100,
                   statistic = rep(c("mean", "cov"), 8),
                   value = 0)
  )

sobs_mean <- plot_method_statistic("sobs", "mean") +
  labs(y = "Mean value\nacross samples",
       title = "Richness") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     labels = c("0", "10,000", "20,000", "30,000"))
  
sobs_cov <- plot_method_statistic("sobs", "cov") +
  labs(y = "Coefficient of\nvariation (%)") +
  scale_y_continuous(limits = c(0, 4))
  
shannon_mean <- plot_method_statistic("shannon", "mean") +
  labs(title = "Shannon")
shannon_cov <- plot_method_statistic("shannon", "cov") +
  scale_y_continuous(limits = c(0, 4))

bray_mean <- plot_method_statistic("bray", "mean") +
  labs(title = "Bray-Curtis")
  
bray_cov <- plot_method_statistic("bray", "cov")

# jaccard_mean <- plot_method_statistic("jaccard", "mean") +
#   labs(title = "Jaccard")

# jaccard_cov <- plot_method_statistic("jaccard", "cov")

patch <- sobs_mean + shannon_mean + bray_mean +
         sobs_cov + shannon_cov + bray_cov +
        plot_layout(ncol = 3, guides = "collect") +
        plot_annotation(caption = "Number of sequences sampled",
                        theme = theme(plot.caption =
                                      element_text(size = 12, hjust = 0.45))) &
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.key = element_blank())


ggsave("figures/intrasample_variation.tiff", patch,
       width = 7, height = 4,
       compression = "lzw+p")

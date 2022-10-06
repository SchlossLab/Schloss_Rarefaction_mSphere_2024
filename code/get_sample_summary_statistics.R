#!/bin/env Rscript

library(tidyverse)

summarize_data <- function(environment) {

	removed_samples <- paste0("data/", environment, "/data.remove_accnos") %>%
		read_lines(.) %>%
		str_split("-") %>%
		unlist() %>%
		enframe(., name=NULL, value="sample") %>%
		mutate(sample = as.character(sample))

	sra_study <- paste0("data/", environment, "/sra_info.tsv") %>%
		read_tsv() %>%
		pull("SRA_Study") %>%
		unique()


	paste0("data/", environment, "/data.group_count") %>%
		read_tsv(col_names=c("sample", "n_seqs"),
						 col_types=cols(sample=col_character(), n_seqs=col_integer())) %>%
		anti_join(., removed_samples, by="sample") %>%
		summarize(
			n_samples=n(), total_seqs=sum(n_seqs), mean=mean(n_seqs),
			min=min(n_seqs), l_quartile = quantile(n_seqs, 0.25),
			median=median(n_seqs), u_quartile = quantile(n_seqs, 0.75),
			max=max(n_seqs),
			n_removed = removed_samples %>% filter(sample != "") %>% nrow(),
			fold_difference = max/min) %>%
		mutate(sra_study = sra_study)
}


tibble(
		directory=c(
			"bioethanol", "human", "lake", "marine",
			"mice", "peromyscus", "rainforest",	"rice",
			"seagrass", "sediment", "soil", "stream"
		),
		nice_name=c(
			"Bioethanol", "Human", "Lake", "Marine",
			"Mice", "Peromyscus", "Rainforest", "Rice",
			"Seagrass", "Sediment", "Soil", "Stream"
		),
		reference=c(
			"[@Li2015]", "[@Baxter2016]", "[@Beall2015]", "[@Henson2016]",
			"[@Kozich2013]", "[@Baxter2014]", "[@LevyBooth2018]", "[@Edwards2015]",
			"[@Ettinger2017]", "[@Graw2018]", "[@Johnston2016]", "[@Hassell2018]"
		)
	) %>%
	group_by(directory) %>%
	mutate(summary_data = map(directory, summarize_data)) %>%
	unnest(summary_data) %>%
	ungroup() %>%
	write_tsv("data/process/study_summary_statistics.tsv")

#!/usr/bin/env Rscript

library("rmarkdown")

if (!tinytex::is_tinytex()) {
  tinytex::install_tinytex(force = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
rmd <- args[1]

render(rmd, output_format = "all", clean = FALSE)

file.remove(gsub("Rmd", "log", rmd))

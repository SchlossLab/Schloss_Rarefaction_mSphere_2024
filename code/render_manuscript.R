#!/usr/bin/env Rscript

library('rmarkdown')

if(!tinytex::is_tinytex()) {
  tinytex::install_tinytex()
} 

render('submission/manuscript.Rmd', clean=FALSE)

file.rename("submission/manuscript.knit.md", "submission/manuscript.md")

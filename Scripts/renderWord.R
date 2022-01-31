#!/usr/bin/env Rscript

bookdown::render_book("index.Rmd","bookdown::word_document2", clean=FALSE)

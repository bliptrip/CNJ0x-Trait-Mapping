#!/usr/bin/env Rscript

bookdown::render_book("index.Rmd","redoc::redoc", clean=FALSE)

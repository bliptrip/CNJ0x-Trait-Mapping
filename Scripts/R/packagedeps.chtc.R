#!/usr/bin/env Rscript
#Only includes R package dependencies needed for running the permutation and other remote code
#on the UW-Madison CHTC condor clusters.

install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

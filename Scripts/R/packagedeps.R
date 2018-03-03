#!/usr/bin/env Rscript

##### For analyzing phenotypes #####
install.packages(c("dplyr"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("pryr"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("openxlsx","RColorBrewer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lm"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

##### For qtl pipeline #####
#Only install the following packages when running for the first time
install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("sommer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lattice"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#IRanges is a part of bioconductor package: ## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")
biocLite("GenomicRanges")

install.packages(c("intervals"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#The following may need to be run as a superuser for the first time on a macosx.  Also, one may manually need to install openmpi.
install.packages(c("Rmpi"), repos = "http://mirror.las.iastate.edu/CRAN/", configure.args="--with-Rmpi-include=/opt/openmpi/include --with-Rmpi-libpath=/opt/openmpi/lib --with-Rmpi-type=OPENMPI", dependencies=TRUE, verbose=TRUE)
install.packages(c("snow"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("doSNOW"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

install.packages(c("jsonlite"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

#For performing clustering of qtls to group them appropriately.
install.packages(c("hclust"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("factoextra"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

#For generating interesting, interactive plots
install.packages(c("plotly"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#For generating correlation-type heatmaps
install.packages(c("GGally"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)


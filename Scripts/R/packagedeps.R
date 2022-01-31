#!/usr/bin/env Rscript

##### For analyzing phenotypes #####
install.packages(c("Hmisc"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("tidyverse"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("pryr"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("car","openxlsx","RColorBrewer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lm"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("seriation"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#
###### For qtl pipeline #####
##Only install the following packages when running for the first time
install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("devtools"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#require(devtools)
#install_version("sommer", version = "4.1.3", repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE) #Install specific version of sommer since every release breaks shit
install.packages(c("sommer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lattice"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#IRanges is a part of bioconductor package: ## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
BiocManager::install(version = "3.12")
BiocManager::install(c("IRanges","GenomicRanges"))
install.packages(c("intervals"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE) #The following may need to be run as a superuser for the first time on a macosx.  Also, one may manually need to install openmpi.
install.packages(c("Rmpi"), repos = "http://mirror.las.iastate.edu/CRAN/", configure.args="--with-Rmpi-include=/opt/openmpi/include --with-Rmpi-libpath=/opt/openmpi/lib --with-Rmpi-type=OPENMPI", dependencies=TRUE, verbose=TRUE)
install.packages(c("snow"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("doSNOW"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

install.packages(c("jsonlite"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#For performing clustering of qtls to group them appropriately.
#For generating interesting, interactive plots
install.packages(c("plotly"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#For generating correlation-type heatmaps
install.packages(c("ggcorrplot"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE) 
install.packages(c("corrplot"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE) 
install.packages(c("GGally"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE) #Themese package for ggplot2
install.packages(c("ggthemes"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("factoextra"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
remotes::install_github("bliptrip/kableExtra")
install.packages(c("formattable"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("ggfittext"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("cowplot"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("rlist"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

remotes::install_github("noamross/redoc")

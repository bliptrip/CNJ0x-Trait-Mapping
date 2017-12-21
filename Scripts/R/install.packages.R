#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script simply installs dependencies needed for QTL analysis
install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("sommer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lattice"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#IRanges is a part of bioconductor package: 
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")
biocLite("GenomicRanges")

install.packages(c("intervals"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#The following may need to be run as a superuser for the first time on a macosx.  Also, one may manually need to install openmpi.
install.packages(c("Rmpi"), repos = "http://mirror.las.iastate.edu/CRAN/", configure.args="--with-Rmpi-include=/opt/openmpi/include --with-Rmpi-libpath=/opt/openmpi/lib --with-Rmpi-type=OPENMPI", dependencies=TRUE, verbose=TRUE)
install.packages(c("snow"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("doSNOW"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

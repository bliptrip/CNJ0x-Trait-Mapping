#!/bin/sh

if [ "$#" != "4" ]; then
    echo "FAILED: Incorrect number of input parameters on commandline! $#" >&2
fi

#Untar the R install and set the path
tar -xzf $5

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH

# Command arguments per process launched are the unique seed per permutation and the number of permutations per process
mtraits=\"$(cat $1)\"
qtlmethod=\"$(cat $2)\"
maxqtl=\"$(cat $3)\"
# run R, with the name of your  R script
R CMD BATCH '--args query_mtraits='$mtraits' query_method='$qtlmethod' query_max_qtl='$maxqtl'' qtl_pipeline_07_stepwiseqtl.R

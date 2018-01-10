#!/bin/sh

if [ "$#" != "7" ]; then
    echo "FAILED: Incorrect number of input parameters on commandline! $#" >&2
fi

#Untar the R install and set the path
tar -xzf $7

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH

# Command arguments per process launched are the unique seed per permutation and the number of permutations per process
seed=$(cat $2)
nperms=$(cat $3)
model=\"$(cat $4)\"
mtraits=\"$(cat $5)\"
qtlmethod=\"$(cat $6)\"
# run R, with the name of your  R script
R CMD BATCH '--args process='$1' seed='$seed' n.perms='$nperms'  query_model='$model'  query_mtraits='$mtraits' query_method='$qtlmethod'' qtl_pipeline_05_runperms2D.R

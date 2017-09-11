#!/bin/sh

if [ $# -neq 7 ]; then
    echo "FAILED: Incorrect number of input parameters on commandline! $#" >&2
fi

#Untar the R install and set the path
tar -xzf $7

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH

#Untar the datafiles, which includes the csv files containing the BLUP-converted phenotypic data for each of the trait types, etc.
tar -xzf qtl_pipeline_05_runperms2D.datafiles.tar.gz

# Command arguments per process launched are the unique seed per permutation and the number of permutations per process
seed=$(cat $2)
nperms=$(cat $3)
mmers=\"$(cat $4)\"
traits=\"$(cat $5)\"
subtraits=\"$(cat $6)\"
# run R, with the name of your  R script
R CMD BATCH '--args process='$1' seed='$seed' n.perms='$nperms'  query_mmers='$mmers'  query_traits='$traits'  query_subtraits='$subtraits'' qtl_pipeline_05_runperms2D.R

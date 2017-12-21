#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

args = commandArgs(trailingOnly=TRUE)

# loading libraries
library(qtl)

if(length(args)==0) {
    print("No arguments supplied.")
    workflow="../../Workflows/1"
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source(paste0(workflow,"/configs/model.cfg"))

traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T, stringsAsFactors=F)
for( i in 1:length(traits.df[,1]) ) {
    trait.cfg       <- traits.df[i,]
    if ( is.na(trait.cfg$mask) || (trait.cfg$mask != "TRUE") ) {
        traits           <- unlist(strsplit(trait.cfg$mtraits,","))
        trait.names      <- paste0(traits,collapse="__")
        trait_subfolder  <- paste0(c(trait.cfg$model,trait.names),collapse="--")
        trait_subfolder_fpath <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        cross <- read.cross(format='csv', file=paste0(trait_subfolder_fpath,"/cross.csv"), genotypes=NULL)
        cross <- calc.genoprob(cross, step=0, map.function="kosambi")
        scan.one  <- scanone(cross, pheno.col=traits, method=qtl_method, verbose=FALSE)
        scan.one.perms <- readRDS(paste0(trait_subfolder_fpath, "/operms.rds"))
        scan.one.summary <- summary(scan.one, perms=scan.one.perms, alpha=0.2, pvalues=TRUE)
        saveRDS(scan.one, file=paste0(trait_subfolder_fpath, "/scanone.rds"), compress=TRUE)
        saveRDS(scan.one.summary, file=paste0(trait_subfolder_fpath, "/scanone.summary.rds"), compress=TRUE)
    }
}

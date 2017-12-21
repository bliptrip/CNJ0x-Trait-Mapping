#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
library(qtl)

args = commandArgs(trailingOnly=TRUE)

workflow="../../Workflows/1"
if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source(paste0(workflow,"/configs/model.cfg"))

#Loop over all mmers, trait groups, subgroups, and perform stepwiseqtl()
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
for( i in 1:length(traits.df[,1]) ) {
    trait.cfg       <- traits.df[i,]
    if ( is.na(trait.cfg$mask) || (trait.cfg$mask != "TRUE") ) {
        traits           <- unlist(strsplit(trait.cfg$mtraits,","))
        trait.names      <- paste0(traits,collapse="__")
        trait_subfolder  <- paste0(c(trait.cfg$model,trait.names),collapse="--")
        trait_subfolder_fpath <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        cross <- read.cross(format='csv', file=paste0(trait_subfolder_fpath,"/cross.csv"), genotypes=NULL)
        cross <- calc.genoprob(cross, step=0, map.function="kosambi")
        scan.two.perms <- readRDS(paste0(trait_subfolder_fpath, "/operms.2D.rds"))
        pens  <- calc.penalties(scan.two.perms)
        for( j in 1:length(traits) ) {
            trait <- traits[j]
            print(paste0("Running stepwiseqtl() with model: ",trait.cfg$model," | traits: ",trait.names, " | trait: ", trait))
            trait_subsubfolder_fpath = paste0(trait_subfolder_fpath, '/', trait) 
            dir.create(trait_subsubfolder_fpath, showWarnings = FALSE)
            scan.sw <- stepwiseqtl(cross,pheno.col=trait,max.qtl=qtl_max,method=qtl_method,penalties=pens[j,],additive.only=FALSE)
            saveRDS(scan.sw, file=paste0(trait_subsubfolder_fpath, "/scansw.rds"), compress=TRUE)
        }
    }
}

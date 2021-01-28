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
source('./usefulFunctions.R')

#Loop over all mmers, trait groups, subgroups, and perform stepwiseqtl()
stepwiseQtlCB <- function(trait.cfg, trait.path, funArgs) {
    cross <- readRDS(file=paste0(trait.path,"/cross.rds"))
    trait <- trait.cfg$trait
    if( any(cross$pheno[trait] != 0) ) {
		scan.two.perms <- readRDS(paste0(trait.path, "/operms.2D.rds"))
		pens  <- calc.penalties(scan.two.perms)
		print(paste0("Running stepwiseqtl() with model: ",trait.cfg$model, " | trait: ", trait))
		trait_subsubfolder_fpath = paste0(trait.path, '/', trait) 
		dir.create(trait_subsubfolder_fpath, showWarnings = FALSE)
		scan.sw <- stepwiseqtl(cross,pheno.col=trait,max.qtl=qtl_max,method=qtl_method,penalties=pens,additive.only=FALSE, keeplodprofile=TRUE, keeptrace=TRUE)
		saveRDS(scan.sw, file=paste0(trait_subsubfolder_fpath, "/scansw.rds"), compress=TRUE)
	}
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, stepwiseQtlCB)

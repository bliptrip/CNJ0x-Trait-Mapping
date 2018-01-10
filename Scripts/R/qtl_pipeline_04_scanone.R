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
source('./usefulFunctions.R')

traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T, stringsAsFactors=F)

merge1D <- function(trait.cfg, trait.names, traits, trait.path, funArgs) {
    cross            <- readRDS(file=paste0(trait.path,"/cross.rds"))
    scan.one         <- scanone(cross, pheno.col=traits, method=qtl_method, verbose=FALSE)
    scan.one.perms   <- readRDS(paste0(trait.path, "/operms.rds"))
    scan.one.summary <- summary(scan.one, perms=scan.one.perms, alpha=0.2, pvalues=TRUE)
    saveRDS(scan.one, file=paste0(trait.path, "/scanone.rds"), compress=TRUE)
    saveRDS(scan.one.summary, file=paste0(trait.path, "/scanone.summary.rds"), compress=TRUE)
}

loopThruTraits(workflow, merge1D)

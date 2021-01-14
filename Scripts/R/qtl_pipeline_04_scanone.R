#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

args = commandArgs(trailingOnly=TRUE)

# loading libraries
library(qtl)

workflow="../../Workflows/8"

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source(paste0(workflow,"/configs/model.cfg"))
source('./usefulFunctions.R')

scanoneCB <- function(trait.cfg, trait.path, funArgs) {
    cross            <- readRDS(file=paste0(trait.path,"/cross.rds"))
    scan.one.perms   <- readRDS(paste0(trait.path, "/operms.rds"))
    trait <- trait.cfg$trait
    print(paste0("Running scanone() with model: ",trait.cfg$model, " | trait: ", trait))
    trait_subsubfolder_fpath = paste0(trait.path, '/', trait) 
    dir.create(trait_subsubfolder_fpath, showWarnings = FALSE)
    scan.one         <- scanone(cross, pheno.col=trait, method=qtl_method, verbose=FALSE)
    scan.one.sum <- summary(scan.one, perms=scan.one.perms[,trait], alpha=scanone_alpha, pvalues=TRUE)
    #Now perform a refineqtl() to make the scanone object compatible with stepwiseqtl() results
    scan.one.qtl <- refineqtl(cross, chr=scan.one.sum$chr, pos=scan.one.sum$pos, keeplodprofile=TRUE, method=qtl_method)
    attr(scan.one.qtl, "formula") <-  paste0("y ~ ",paste0(scan.one.qtl$altname,collapse=" + "))
    saveRDS(scan.one, file=paste0(trait_subsubfolder_fpath, "/scanone.rds"), compress=TRUE)
    saveRDS(scan.one.sum, file=paste0(trait_subsubfolder_fpath, "/scanone.summary.rds"), compress=TRUE)
    saveRDS(scan.one.qtl, file=paste0(trait_subsubfolder_fpath, "/scanone.qtl.rds"), compress=TRUE)
}

loopThruTraits(workflow, scanoneCB)

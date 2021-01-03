#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on HTCondor clusters and is meant to run the stepwiseqtl() function.

args = commandArgs(trailingOnly=TRUE)

a_process  <- 0
a_mmer	   <- "2011"
a_traits   <- "total_berry_weight"
a_subtrait <- "total_berry_weight"
a_max_qtl  <- 10
a_method   <- "hk"
if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

#Read in the scantwo file 
cross  <- readRDS(file="cross.rds")

#Extract the penalties and run stepwiseqtl()
scan.two.perms <- readRDS(paste0(trait.path, "operms.2D.rds"))
pens           <- calc.penalties(scan.two.perms)
outsw		   <- stepwiseqtl(cross,pheno.col=a_subtrait,max.qtl=a_max_qtl,method=a_method,penalties=pens,additive.only=FALSE, keeplodprofile=TRUE, keeptrace=TRUE)

saveRDS(outsw, filename, compress=TRUE)

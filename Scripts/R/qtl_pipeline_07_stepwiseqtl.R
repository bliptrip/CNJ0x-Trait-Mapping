#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on HTCondor clusters and is meant to run the stepwiseqtl() function.
library(qtl)

args = commandArgs(trailingOnly=TRUE)

query_mtraits   <- "total_berry_weight"
query_max_qtl  <- 10
query_method   <- "hk"
if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

#Read in the scantwo file 
traits         <- unlist(strsplit(query_mtraits,','))
traits.string  <- paste(traits,collapse="__")
cnjpop.cross   <- (read.cross(format = "csv", file='cross.csv', genotypes = NULL))
cnjpop.cross   <- calc.genoprob(cnjpop.cross,step=0,map.function="kosambi") 
scan.two.perms <- readRDS("operms.2D.rds")

#May need to convert this first
pens              <- calc.penalties(scan.two.perms)
scan.sw           <- stepwiseqtl(cnjpop.cross,pheno.col=traits,max.qtl=query_max_qtl,method=query_method,penalties=penalties,additive.only=FALSE)
scan.sw.summary   <- summary(scan.sw, perms=scan.two.perms, alpha=0.2, pvalues=TRUE)

#Save the output files
saveRDS(scan.sw, file="scansw.rds", compressed=TRUE)
saveRDS(scan.sw.summary, file="scansw.summary.rds", compressed=TRUE)

#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

args = commandArgs(trailingOnly=TRUE)

#Set defaults
process <- 0
seed <- 54955149
n.perms <- 2
query_model <- "2011"
query_mtraits <- "total_berry_weight"
query_method <- "hk"

#Parse arguments, if provided
if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

query_mtraits.l.v <- unlist(strsplit(query_mtraits, split=","))

# loading libraries
library(qtl)

#Create a mmer list containing all the query_model.v plus fields for the process ID on the cluster, seed used, the number of permutations, a timestamp, system information (uname -a), and then the mmers for analysis
#NOTE: Rather than using a description field for mmers, I use the  colname to designate the description of the multivariate analysis.
operms <- vector("list",7)
names(operms) <- c("seed","n.perms","process","time.start", "time.end", "uname", "perms")
operms$seed      <- seed
operms$n.perms   <- n.perms
operms$process   <- process
operms$time.start <- Sys.time()
operms$uname     <- system("uname -a", intern=TRUE)

set.seed(seed)
#Read in the cross files for R/qtl and then calculate the genotype probabilities
#Generate a cross datastructure for running analyses on
cnjpop.cross <- readRDS(file="cross.rds")
#For the current analysis
#Derive the columns of interest for running the permutations on.
operms[["perms"]] <- scantwo(cnjpop.cross, n.perm=n.perms, pheno.col=query_mtraits.l.v, method=query_method, verbose=FALSE)
operms$time.end <- Sys.time()

rds.filename <- paste0("operms.2D.",as.character(process),".rds")
saveRDS(operms,file=rds.filename,compress=TRUE)

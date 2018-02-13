#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

# loading libraries
library(qtl)

args = commandArgs(trailingOnly=TRUE)

#Default command-line value
workflow="../../Workflows/1"

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

folder <- paste0(workflow,"/perms2D")
dirs <- list.dirs(path=folder,full.names=FALSE,recursive=FALSE)

#Derive all models from the paths
for ( i in 1:length(dirs) ) {
    dir <- dirs[[i]]
    subdirs <- list.dirs(path=paste0(folder,'/',dirs[[i]]),full.names=F)
    #The first entry is the current path, so skip it.
    subdirs <- subdirs[2:length(subdirs)]
    subdirs_fpath <- list.dirs(path=paste0(folder,'/',dirs[[i]]),full.names=T)
    #The first entry is the current path, so skip it.
    subdirs_fpath <- subdirs_fpath[2:length(subdirs_fpath)]
    operms.sub <- readRDS(paste0(subdirs_fpath[[1]],'/operms.2D.',as.character(subdirs[[1]]),'.rds'))
    operms <- operms.sub[["perms"]]
    for ( j in 2:length(subdirs) ) {
        operms.sub <- readRDS(paste0(subdirs_fpath[[j]],'/operms.2D.',as.character(subdirs[[j]]),'.rds'))
        operms <- c(operms,operms.sub[["perms"]])
    }
    saveRDS(operms,paste0(workflow,'/traits/',dir,'/operms.2D.rds'),compress=TRUE)
}

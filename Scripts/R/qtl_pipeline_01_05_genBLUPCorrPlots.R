#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script generates the correlation plots derived from means of two Vorsa populations along with the correlations of 
# both populations combined..

# loading libraries
source('./usefulFunctions.R')
library(tidyverse)

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

blups.collated.wide.tb <- readRDS(paste0(workflow,"/traits/blups_collated.wide.rds"))
blup_raw_cors <- blups.collated.wide.tb %>% group_by(model,trait) %>% select(blup,raw) %>% summarize(brcorr=cor(blup,raw,use="pairwise.complete.obs",method="spearman"))
write_csv(blup_raw_cors, file=paste0(workflow,"/traits/blups_raw_corrs.csv"))

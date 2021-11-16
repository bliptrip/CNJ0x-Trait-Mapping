#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

# loading libraries
source('./usefulFunctions.R')

library(dplyr)

#Read the phenotypic data in that was output by analyze_remove_outliers.R

cnjpop.pheno.df <- readRDSw('cnjpop.noout.df.rds')
cnjpop.pheno.df$upright_length = as.numeric(cnjpop.pheno.df$upright_length)
cnjpop.pheno.df$secondary_growth = as.numeric(cnjpop.pheno.df$secondary_growth)
cnjpop.pheno.df$dry_wt_leaves = as.numeric(cnjpop.pheno.df$dry_wt_leaves)

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
    #If pheno.file is specified on command-line, then specify it here.
    if( !is.na(pheno.file) ) {
        cnjpop.pheno.df <- readRDS(pheno.file)
    }
}

extract_row <- function(accession) {
    return(gsub("R([[:digit:]]+)C([[:digit:]]+)","\\1",accession,ignore.case=T))
}

extract_col <- function(accession) {
    return(gsub("R([[:digit:]]+)C([[:digit:]]+)","\\2",accession,ignore.case=T))
}

#Setup which phenotypes are invalid/missing by filling in with 'NA'
cnjpop.pheno.df.length_width_weight.zeros.idx <- which(cnjpop.pheno.df$berry_length == 0 | cnjpop.pheno.df$berry_width == 0 | cnjpop.pheno.df$berry_weight == 0)
cnjpop.pheno.df$berry_length[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df$berry_width[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df$berry_weight[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df[cnjpop.pheno.df.length_width_weight.zeros.idx,]
cnjpop.pheno.df.num_seeds.zeros.idx <- which(cnjpop.pheno.df$num_seeds == 0)
cnjpop.pheno.df$num_seeds[cnjpop.pheno.df.num_seeds.zeros.idx] <- NA
cnjpop.pheno.df.upright_length.zeros.idx <- which(cnjpop.pheno.df$upright_length == 0)
cnjpop.pheno.df$upright_length[cnjpop.pheno.df.upright_length.zeros.idx] <- NA
cnjpop.pheno.df.secondary_growth.zeros.idx <- which(cnjpop.pheno.df$secondary_growth == 0)
cnjpop.pheno.df$secondary_growth[cnjpop.pheno.df.secondary_growth.zeros.idx] <- NA
cnjpop.pheno.df.dry_wt_leaves.zeros.idx <- which(cnjpop.pheno.df$dry_wt_leaves == 0)
cnjpop.pheno.df$dry_wt_leaves[cnjpop.pheno.df.dry_wt_leaves.zeros.idx] <- NA
#Add row/column designators to indicate the location                            
cnjpop.pheno.df <- cnjpop.pheno.df %>%
                      mutate(row=extract_row(accession),
                             column=extract_col(accession),
							 accession_name=gsub("[- ]","_",accession_name,fixed=FALSE)) #Replace spaces/hyphens with underscores 

#Calculate means of uprights
cnjpop.pheno.means.df <- cnjpop.pheno.df %>%
                            group_by(population,year,accession_name,accession,row,column) %>%
                            summarize(across(num_peds:chimera_umccLogY,mean,na.rm=T))

#Filter out parents for comparing parental values to children values
#cnjpop.pheno.parents.df <- cnjpop.pheno.df %>% filter(!grepl("CNJ0.*", cnjpop.pheno.df$accession_name))
#Now filter out progeny
#cnjpop.pheno.progeny.df <- cnjpop.pheno.df %>% filter(grepl("CNJ0.*", cnjpop.pheno.df$accession_name))

#Separate two populations for analysis
#Population 1						
cnjpop.pheno.p1.means.df <- cnjpop.pheno.means.df %>%
                                filter(population == 1)
#Remove the population column
cnjpop.pheno.p1.means.df <- subset(cnjpop.pheno.p1.means.df, select=-c(population))                                

#Population 2
cnjpop.pheno.p2.means.df <- cnjpop.pheno.means.df %>%
                                filter(population == 2)
#Remove the population column
cnjpop.pheno.p2.means.df <- subset(cnjpop.pheno.p2.means.df, select=-c(population))                                

#Write the means to an output csv file.

write.csv(cnjpop.pheno.means.df, file=pheno_dpath2fpath("Data-combined-collated.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p1.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj04.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p2.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj02.means.csv"), row.names=FALSE)

save.image(paste0(workflow,"/.RData.1_means"))
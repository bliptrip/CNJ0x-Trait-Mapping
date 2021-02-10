#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

# loading libraries
source('./usefulFunctions.R')

library(dplyr)

#Read the phenotypic data in that was output by analyze_remove_outliers.R

cnjpop.pheno.df <- readRDSw('cnjpop.noout.df.rds')

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
#Add row/column designators to indicate the location                            
cnjpop.pheno.df <- cnjpop.pheno.df %>%
                      mutate(row=extract_row(accession),
                             column=extract_col(accession),
							 accession_name=gsub(" ","_",accession_name)) #Remove spaces from accession names as dataframe columns can't have them

#Calculate means of uprights
cnjpop.pheno.means.df <- cnjpop.pheno.df %>%
                            group_by(population,year,accession_name,accession,row,column) %>%
                            summarize(berry_length=mean(berry_length, na.rm=T), 
                                      berry_width=mean(berry_width, na.rm=T),
                                      berry_weight=mean(berry_weight, na.rm=T),
                                      num_seeds=mean(num_seeds, na.rm=T),
                                      num_peds=mean(num_peds, na.rm=T),
                                      num_berries=mean(num_berries, na.rm=T),
                                      total_berry_weight=mean(total_berry_weight, na.rm=T))

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

############### NOTE: Maybe consider removing the 'complete gentoype' code below when considering across each year itself. #####################

#Determine which genotypes are represented across all three years
#p1
#cnjpop.pheno.p1.means.table <- table(cnjpop.pheno.p1.means.df$accession_name) 
#Complete genotypes should apply for both full and means data
#cnjpop.pheno.p1.CompleteGenoTypes <- names(which(cnjpop.pheno.p1.means.table == 3))
#cnjpop.pheno.p1.means.df <- cnjpop.pheno.p1.means.df[cnjpop.pheno.p1.means.df$accession_name %in% cnjpop.pheno.p1.CompleteGenoTypes,]
#p2
#cnjpop.pheno.p2.means.table <- table(cnjpop.pheno.p2.means.df$accession_name) 
#Complete genotypes should apply for both full and means data
#cnjpop.pheno.p2.CompleteGenoTypes <- names(which(cnjpop.pheno.p2.means.table == 3))
#cnjpop.pheno.p2.means.df <- cnjpop.pheno.p2.means.df[cnjpop.pheno.p2.means.df$accession_name %in% cnjpop.pheno.p2.CompleteGenoTypes,]

#Write the means to an output csv file.
write.csv(cnjpop.pheno.means.df, file=pheno_dpath2fpath("Data-combined-collated.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p1.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj04.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p2.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj02.means.csv"), row.names=FALSE)

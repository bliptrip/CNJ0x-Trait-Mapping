#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

# loading libraries
source('./usefulFunctions.R')

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

#First thing first, only consider progeny
cnjpop.pheno.progeny.idx <- which(grepl("CNJ0.*", cnjpop.pheno.df$accession_name))
cnjpop.pheno.df <- cnjpop.pheno.df[cnjpop.pheno.progeny.idx,]

#Setup which phenotypes are invalid
cnjpop.pheno.df.length_width_weight.zeros.idx <- which(cnjpop.pheno.df$berry_length == 0 | cnjpop.pheno.df$berry_width == 0 | cnjpop.pheno.df$berry_weight == 0)
cnjpop.pheno.df$berry_length[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df$berry_width[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df$berry_weight[cnjpop.pheno.df.length_width_weight.zeros.idx] <- NA
cnjpop.pheno.df[cnjpop.pheno.df.length_width_weight.zeros.idx,]
cnjpop.pheno.df.num_seeds.zeros.idx <- which(cnjpop.pheno.df$num_seeds == 0)
cnjpop.pheno.df$num_seeds[cnjpop.pheno.df.num_seeds.zeros.idx] <- NA

#Compute means for each genotype x year across uprights
cnjpop.pheno.berry_length.means.df <- aggregate(formula=berry_length~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.berry_width.means.df <- aggregate(formula=berry_width~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.berry_weight.means.df <- aggregate(formula=berry_weight~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.num_seeds.means.df <- aggregate(formula=num_seeds~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.num_peds.means.df <- aggregate(formula=num_peds~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.num_berries.means.df <- aggregate(formula=num_berries~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)
cnjpop.pheno.total_berry_weight.means.df <- aggregate(formula=total_berry_weight~population+year+accession_name,data=cnjpop.pheno.df,FUN=mean,na.omit=T)


#Combine means into one dataframe
#cnjpop.pheno.means.df <- cbind(cnjpop.pheno.berry_length.means.df, berry_width=cnjpop.pheno.berry_width.means.df[,4], berry_weight=cnjpop.pheno.berry_weight.means.df[,4], num_seeds=cnjpop.pheno.num_seeds.means.df[,4], num_peds=cnjpop.pheno.num_peds.means.df[,4], num_berries=cnjpop.pheno.num_berries.means.df[,4], total_berry_weight=cnjpop.pheno.total_berry_weight.means.df[,4])
cnjpop.pheno.merge.df <- merge(cnjpop.pheno.berry_length.means.df, cnjpop.pheno.berry_width.means.df, by=c("population", "year", "accession_name"))
cnjpop.pheno.merge.df <- merge(cnjpop.pheno.merge.df, cnjpop.pheno.berry_weight.means.df, by=c("population", "year", "accession_name"))
cnjpop.pheno.merge.df <- merge(cnjpop.pheno.merge.df, cnjpop.pheno.num_seeds.means.df,by=c("population", "year", "accession_name"))
cnjpop.pheno.merge.df <- merge(cnjpop.pheno.merge.df, cnjpop.pheno.num_peds.means.df,by=c("population", "year", "accession_name"))
cnjpop.pheno.merge.df <- merge(cnjpop.pheno.merge.df, cnjpop.pheno.num_berries.means.df,by=c("population", "year", "accession_name"))
cnjpop.pheno.means.df <- merge(cnjpop.pheno.merge.df, cnjpop.pheno.total_berry_weight.means.df, by=c("population", "year", "accession_name"))
#colnames(cnjpop.pheno.means.df) <- c("population", "year", "accession_name", "berry_length", "berry_width", "berry_weight", "num_seeds", "num_peds", "num_berries", "total_berry_weight")

#Separate two populations for analysis
#Means
cnjpop.pheno.p1.means.df.idx <- which(cnjpop.pheno.means.df$population=="1")
cnjpop.pheno.p1.means.df <- cnjpop.pheno.means.df[cnjpop.pheno.p1.means.df.idx,-1] #remove the population field as no longer needed
cnjpop.pheno.p2.means.df.idx <- which(cnjpop.pheno.means.df$population=="2")
cnjpop.pheno.p2.means.df <- cnjpop.pheno.means.df[cnjpop.pheno.p2.means.df.idx,-1] #Remove the population field as no longer needed
#All values
cnjpop.pheno.p1.df.idx <- which(cnjpop.pheno.df$population=="1")
cnjpop.pheno.p1.df <- cnjpop.pheno.df[cnjpop.pheno.p1.df.idx,-1] #Remove the population field as no longer needed
cnjpop.pheno.p2.df.idx <- which(cnjpop.pheno.df$population=="2")
cnjpop.pheno.p2.df <- cnjpop.pheno.df[cnjpop.pheno.p2.df.idx,-1] #Remove the population field as no longer needed

############### NOTE: Maybe consider removing the 'complete gentoype' code below when considering across each year itself. #####################

#Determine which genotypes are represented across all three years
#p1
cnjpop.pheno.p1.means.table <- table(cnjpop.pheno.p1.means.df$accession_name) 
#Complete genotypes should apply for both full and means data
cnjpop.pheno.p1.CompleteGenoTypes <- names(which(cnjpop.pheno.p1.means.table == 3))
cnjpop.pheno.p1.df <- cnjpop.pheno.p1.df[cnjpop.pheno.p1.df$accession_name %in% cnjpop.pheno.p1.CompleteGenoTypes,]
cnjpop.pheno.p1.means.df <- cnjpop.pheno.p1.means.df[cnjpop.pheno.p1.means.df$accession_name %in% cnjpop.pheno.p1.CompleteGenoTypes,]
#p2
cnjpop.pheno.p2.means.table <- table(cnjpop.pheno.p2.means.df$accession_name) 
#Complete genotypes should apply for both full and means data
cnjpop.pheno.p2.CompleteGenoTypes <- names(which(cnjpop.pheno.p2.means.table == 3))
cnjpop.pheno.p2.df <- cnjpop.pheno.p2.df[cnjpop.pheno.p2.df$accession_name %in% cnjpop.pheno.p2.CompleteGenoTypes,]
cnjpop.pheno.p2.means.df <- cnjpop.pheno.p2.means.df[cnjpop.pheno.p2.means.df$accession_name %in% cnjpop.pheno.p2.CompleteGenoTypes,]

#Write the means to an output csv file.
write.csv(cnjpop.pheno.means.df, file=pheno_dpath2fpath("Data-combined-collated.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p1.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj04.means.csv"), row.names=FALSE)
write.csv(cnjpop.pheno.p2.means.df, file=pheno_dpath2fpath("Data-combined-collated.cnj02.means.csv"), row.names=FALSE)

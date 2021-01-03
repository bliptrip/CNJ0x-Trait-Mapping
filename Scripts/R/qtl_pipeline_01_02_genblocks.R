#!/usr/bin/env RScript
#Try to do two-stage geographical clustering on genotypes & phenotypes to determine any phenotype blocking effects.

# loading libraries
source('./usefulFunctions.R')

workflow <- "../../Workflows/1"

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

library(qtl)
library(ClustGeo)

#Load the configuration file for this workflow
source(paste0(workflow,"/configs/model.cfg"))

###### Stage 1: Cluster on genotypes using realized additive matrix as the 'distance' matrix. ###### 
#In this case, we are wanting to cluster based on heterogeneity of the genotypes, to make sure we 
#have representative genotype dispersion within our planting.  Subsequently, we'll average these 
#genotypic blocks phenotypes and then calculate a blocking effect at a higher level to find
#uniform blocks.
####################################################################################################

#Read in the phenotype file
pheno.means.df<-read.csv(file=pheno_dpath2fpath(pheno_file))

#Read in the marker maps
geno<-read.table(geno_rpath2fpath(geno_file),header=T,sep=',')  
rownames(geno)<-geno$X
geno<-geno[,-c(1:6)] #Remove the marker name, segregation pattern, phase, classification, position, and lg fields -- only keep the genotype calls
geno.num <-atcg1234(t(geno))

#Calculate the realized additive relationship matrix.  Note: This is D0, or the feature matrix
#distances, for ClustGeo.
D0  <- A.mat(geno.num)
D1  <- dist(pheno.means.df[,c("row","column")])

#Read in the model trait configuration file to determine how to model traits.
traits.df <- as.matrix(read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv")))



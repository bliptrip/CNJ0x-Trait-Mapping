#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script generates the correlation plots derived from means of two Vorsa populations along with the correlations of 
# both populations combined..

# loading libraries
source('./usefulFunctions.R')

library(lattice)

cnjpop.pheno.means.df <- read.csv(file=pheno_dpath2fpath("Data-combined-collated.means.csv"))
cnjpop.pheno.p1.means.df <- read.csv(file=pheno_dpath2fpath("Data-combined-collated.cnj04.means.csv"))
cnjpop.pheno.p2.means.df <- read.csv(file=pheno_dpath2fpath("Data-combined-collated.cnj02.means.csv"))

#Function: split_by_year()
#
# Purpose: To split the phenotypic dataframe into separate years for generating a phenotypic correlation matrix.
#
# Args: 
#     - pheno.df: The input dataframe.
#     - year.col: The name for the column specifying the year.
#     - vars.col: A vector containing the list of columns we care to split into separate years.
#     - by: The name of the column we want to merge on.  NOTE: This will usually be the genotype, cultivar, or accession name.
#
# Returns: A new dataframe containing columns named var.year for each variable by year combination, along organized by
#           the 'by' column (typically genotype, cultivar, or accession name).
#
# NOTE: The column specified by year.col should be a factor.
#
split_by_year <- function(pheno.df, year.col, vars.col, by) {
genos <- unique(pheno.df[,by])
split.pheno.df   <- data.frame(genos)
colnames(split.pheno.df) <- c(by)
for( var in vars.col ) {
    for( year in unique((pheno.df[,year.col])) ) {
        single_year.idx <- which(pheno.df[,year.col] == year) 
        newvar <- paste0(var,".",year)
        single_year.df <- pheno.df[single_year.idx,c(by,var)]
            colnames(single_year.df) <- c(by,newvar)
            split.pheno.df <- merge(split.pheno.df,single_year.df,by=by)
        }
    }
    return(split.pheno.df)
}

#This could potentially be an input from a configuration file, but in general, these initial scripts are specific to the project.
focal.cols <- c("berry_length", "berry_width", "berry_weight", "num_seeds", "num_peds", "num_berries", "total_berry_weight")

#Do a combined population assessment
cnjpop.pheno.combined.means.cor.mat <- cor(cnjpop.pheno.means.df[,focal.cols])
write.csvw(cnjpop.pheno.combined.means.cor.mat, "combinedpopulations_phenotype_correlations_global.csv")

#Generate the overall correlation matrix between variables across all years
cnjpop.pheno.p1.means.cor.mat <- cor(cnjpop.pheno.p1.means.df[,focal.cols])
write.csvw(cnjpop.pheno.p1.means.cor.mat, "CNJ04_phenotype_correlations_global.csv")
cnjpop.pheno.p2.means.cor.mat <- cor(cnjpop.pheno.p2.means.df[,focal.cols])
write.csvw(cnjpop.pheno.p2.means.cor.mat, "CNJ02_phenotype_correlations_global.csv")

#Generate the correlation matrix by with variables split year
cnjpop.pheno.p1.means.split.df <- split_by_year(cnjpop.pheno.p1.means.df, year.col="year", vars.col=focal.cols, by="accession_name")
cnjpop.pheno.p1.means.cor.split.mat <- cor(cnjpop.pheno.p1.means.split.df[,-1])
write.csvw(cnjpop.pheno.p1.means.cor.split.mat, "CNJ04_phenotype_correlations_by_year.csv")

cnjpop.pheno.p2.means.split.df <- split_by_year(cnjpop.pheno.p2.means.df, year.col="year", vars.col=focal.cols, by="accession_name")
cnjpop.pheno.p2.means.cor.split.mat <- cor(cnjpop.pheno.p2.means.split.df[,-1])
write.csvw(cnjpop.pheno.p2.means.cor.split.mat, "CNJ02_phenotype_correlations_by_year.csv")

#graph.cols <- colorRampPalette(c(rgb(0,0,1),rgb(1,0,0)))(100)
graph.cols <- colorRampPalette(c('white','black'))(100)
postscriptw(file="p12_phenotypes.wholeCor.eps", title="Phenotypic Correlations")
levelplot(cnjpop.pheno.combined.means.cor.mat,scales=list(x=list(rot=90)),main='Combined Populations', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p1.means.cor.mat,scales=list(x=list(rot=90)),main='CNJ04', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p1.means.cor.split.mat,scales=list(x=list(rot=90)),main='CNJ04: Year Split', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p2.means.cor.mat,scales=list(x=list(rot=90)),main='CNJ02', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p2.means.cor.split.mat,scales=list(x=list(rot=90)),main='CNJ02: Year Split', xlab='',ylab='',col.regions=graph.cols)
dev.off()

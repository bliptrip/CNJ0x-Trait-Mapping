#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script generates the correlation plots derived from means of two Vorsa populations along with the correlations of 
# both populations combined..

# loading libraries
source('./usefulFunctions.R')

library(Hmisc) #For rcorr -- does pairwise correlations with p-value stats included
library(RColorBrewer)
library(corrplot)
library(GGally)
library(ggthemes)
library(lattice)
library(plotly)

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
        newvar <- paste0(var,"y",substr(year,3,4))
        single_year.df <- pheno.df[single_year.idx,c(by,var)]
            colnames(single_year.df) <- c(by,newvar)
            split.pheno.df <- merge(split.pheno.df,single_year.df,by=by)
        }
    }
    return(split.pheno.df)
}



#This could potentially be an input from a configuration file, but in general, these initial scripts are specific to the project.
focal.cols <- c("berry_length", "berry_width", "berry_weight", "num_seeds", "num_peds", "num_berries", "total_berry_weight")
colnames(cnjpop.pheno.means.df) <- rename_traits(colnames(cnjpop.pheno.means.df), trait.abbrev.map.df)
colnames(cnjpop.pheno.p1.means.df) <- rename_traits(colnames(cnjpop.pheno.p1.means.df), trait.abbrev.map.df)
colnames(cnjpop.pheno.p2.means.df) <- rename_traits(colnames(cnjpop.pheno.p2.means.df), trait.abbrev.map.df)
focal.cols <- rename_traits(focal.cols, trait.abbrev.map.df)

#Do a combined population assessment
cnjpop.pheno.combined.means.cor.mat <- rcorr(as.matrix(cnjpop.pheno.means.df[,focal.cols]), type="pearson")
#write.csvw(cnjpop.pheno.combined.means.cor.mat, "combinedpopulations_phenotype_correlations_global.csv")

#Generate the overall correlation matrix between variables across all years
cnjpop.pheno.p1.means.sub.df <- cnjpop.pheno.p1.means.df[,focal.cols]
cnjpop.pheno.p1.means.cor.mat <- rcorr(as.matrix(cnjpop.pheno.p1.means.df[,focal.cols]), type="pearson")
#write.csvw(cnjpop.pheno.p1.means.cor.mat, "CNJ04_phenotype_correlations_global.csv")
cnjpop.pheno.p2.means.cor.mat <- rcorr(as.matrix(cnjpop.pheno.p2.means.df[,focal.cols]), type="pearson")
#write.csvw(cnjpop.pheno.p2.means.cor.mat, "CNJ02_phenotype_correlations_global.csv")

#Generate the correlation matrix by with variables split year
cnjpop.pheno.p1.means.split.df <- split_by_year(cnjpop.pheno.p1.means.df, year.col="year", vars.col=focal.cols, by="accession_name")
cnjpop.pheno.p1.means.cor.split.mat <- rcorr(as.matrix(cnjpop.pheno.p1.means.split.df[,-1]), type="pearson")
#write.csvw(cnjpop.pheno.p1.means.cor.split.mat, "CNJ04_phenotype_correlations_by_year.csv")
#cnjpop.pheno.p1.means.cor.split.mat <- read.csvw("CNJ04_phenotype_correlations_by_year.csv")


cnjpop.pheno.p2.means.sub.df <- cnjpop.pheno.p2.means.df[,focal.cols]
cnjpop.pheno.p2.means.split.df <- split_by_year(cnjpop.pheno.p2.means.df, year.col="year", vars.col=focal.cols, by="accession_name")
cnjpop.pheno.p2.means.cor.split.mat <- rcorr(as.matrix(cnjpop.pheno.p2.means.split.df[,-1]), type="pearson")
#write.csvw(cnjpop.pheno.p2.means.cor.split.mat, "CNJ02_phenotype_correlations_by_year.csv")
#cnjpop.pheno.p2.means.cor.split.mat <- read.csvw("CNJ02_phenotype_correlations_by_year.csv", row.names=1)

#graph.cols <- colorRampPalette(c(rgb(0,0,1),rgb(1,0,0)))(100)
graph.cols <- colorRampPalette(c('white','black'))(100)
postscriptw(file="p12_phenotypes.wholeCor.eps", title="Phenotypic Correlations")
levelplot(cnjpop.pheno.combined.means.cor.mat$r,scales=list(x=list(rot=90)),main='Combined Populations', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p1.means.cor.mat$r,scales=list(x=list(rot=90)),main='CNJ04', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p1.means.cor.split.mat$r,scales=list(x=list(rot=90)),main='CNJ04: Year Split', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p2.means.cor.mat$r,scales=list(x=list(rot=90)),main='CNJ02', xlab='',ylab='',col.regions=graph.cols)
levelplot(cnjpop.pheno.p2.means.cor.split.mat$r,scales=list(x=list(rot=90)),main='CNJ02: Year Split', xlab='',ylab='',col.regions=graph.cols)
dev.off()

#Use ggplot and plotly to develop interactive graphs of phenotypes, correlations b/w trait phenotypes, and heritabilities of traits.  It would also be interesting to generate
#correlations b/w trait BLUPs, but will forget about that for now.
#Plot correlations as a raster heatmap
g <- ggcorr(cnjpop.pheno.p2.means.split.df, label=FALSE, geom="circle", hjust=0.575, min_size=1, max_size=8, size=3, label_color=white) +
         theme(
            text = element_text(color="white"),
            axis.text.x = element_text(color="white"),
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            , panel.grid.major = element_blank() # get rid of major grid
            , panel.grid.minor = element_blank() # get rid of minor grid
            , legend.background = element_rect(fill = "transparent") # get rid of legend bg
            , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
         )
png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait_year.corr.png"), width="1365", height="1024", units="px", bg="transparent")
ggsave(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait.corr.png"), plot=g, device="png", background="transparent")
#png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait_year.corr.png"), width="1280", height="1024", units="px", bg="transparent")
#g
#dev.off()

#Also try using the corrplot package.
#Trait by year association
png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait_year.corrplot.png"), bg="transparent", width=1024, height=1024, units="px")
#corrplot(cnjpop.pheno.p2.means.cor.split.mat$r, p.mat=cnjpop.pheno.p2.means.cor.split.mat$P, method="square", order="hclust", addrect=4, insig="blank", tl.col="white", tl.srt=45, tl.cex=0.7, rect.col="white", pch="N", bg="transparent", addgrid.col="white")
#corrplot(cnjpop.pheno.p2.means.cor.split.mat$r, p.mat=cnjpop.pheno.p2.means.cor.split.mat$P, method="square", order="hclust", addrect=4, insig="blank", tl.col="white", tl.srt=45, tl.cex=0.7, rect.col="white", pch="N", bg="transparent", addgrid.col="white", col=brewer.pal(11,"Spectral"))
corrplot(cnjpop.pheno.p2.means.cor.split.mat$r, p.mat=cnjpop.pheno.p2.means.cor.split.mat$P, method="square", order="hclust", addrect=4, insig="blank", tl.col="black", tl.srt=45, tl.cex=0.8, bg="transparent")
dev.off()

#Trait association
png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait.corrplot.png"), bg="transparent", width=1024, height=1024, units="px")
#corrplot(cnjpop.pheno.p2.means.cor.split.mat$r, p.mat=cnjpop.pheno.p2.means.cor.split.mat$P, method="square", order="hclust", addrect=4, insig="blank", tl.col="white", tl.srt=45, tl.cex=0.7, rect.col="white", pch="N", bg="transparent", addgrid.col="white")
#corrplot(cnjpop.pheno.p2.means.cor.split.mat$r, p.mat=cnjpop.pheno.p2.means.cor.split.mat$P, method="square", order="hclust", addrect=4, insig="blank", tl.col="white", tl.srt=45, tl.cex=0.7, rect.col="white", pch="N", bg="transparent", addgrid.col="white", col=brewer.pal(11,"Spectral"))
corrplot(cnjpop.pheno.p2.means.cor.mat$r, p.mat=cnjpop.pheno.p2.means.cor.mat$P, method="square", order="hclust", addrect=3, insig="blank", tl.col="black", tl.srt=45, tl.cex=0.8, bg="transparent")
dev.off()


#All traits
colnames(cnjpop.pheno.p2.means.sub.df) <- trait.abbrev.map.df[colnames(cnjpop.pheno.p2.means.sub.df), 1]
#        theme_solarized_2(light = FALSE) +
g <- ggcorr(cnjpop.pheno.p2.means.sub.df, label=FALSE, geom="circle", hjust=0.575, min_size=1, max_size=16, size=5, label_color=white) +
         theme(
            text = element_text(color="white"),
            axis.text.x = element_text(color="white"),
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            , panel.grid.major = element_blank() # get rid of major grid
            , panel.grid.minor = element_blank() # get rid of minor grid
            , legend.background = element_rect(fill = "transparent") # get rid of legend bg
            , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
         )
png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait.corr.png"), bg="transparent", width=1024, height=1024, units="px")
g
dev.off()
#ggsave(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p12_phenotypes.trait.corr.png"), plot=g, device="png", background="transparent")

#Generate a scatterplot of the CNJ02 population's traits of interest
g <- ggpairs(cnjpop.pheno.p2.means.sub.df) +
         theme(
            text = element_text(color="white"),
            axis.text.x = element_text(color="white"),
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            , panel.grid.major = element_blank() # get rid of major grid
            , panel.grid.minor = element_blank() # get rid of minor grid
            , legend.background = element_rect(fill = "transparent") # get rid of legend bg
            , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
         )
png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"/p13_phenotypes.trait.pwscatter.png"), bg="transparent")
g
dev.off()

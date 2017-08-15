#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#Only install the following packages when running for the first time
install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("sommer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lattice"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#IRanges is a part of bioconductor package: 
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")
biocLite("GenomicRanges")

install.packages(c("intervals"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#The following may need to be run as a superuser for the first time on a macosx.  Also, one may manually need to install openmpi.
install.packages(c("Rmpi"), repos = "http://mirror.las.iastate.edu/CRAN/", configure.args="--with-Rmpi-include=/opt/openmpi/include --with-Rmpi-libpath=/opt/openmpi/lib --with-Rmpi-type=OPENMPI", dependencies=TRUE, verbose=TRUE)
install.packages(c("snow"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("doSNOW"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

# loading libraries
source('./usefulFunctions.R')


library(lme4)
library(qtl)
library(sommer)
library(IRanges)
library(GenomicRanges) 
library(intervals)
library(snow)
library(doSNOW)
library(lattice)
library(RColorBrewer)

require(openxlsx)

#convert_to_factors <- function(pheno.df) {
#    if( "population" %in% colnames(pheno.df) ) {
#        pheno.df$population <- as.factor(pheno.df$population)
#    }
#    pheno.df$year <- as.factor(pheno.df$year)
#
#    return(pheno.df)
#}

#Luis recommends running scan one on the dataset first

#Read the phenotypic data in that was output by analyze_remove_outliers.R
cnjpop.pheno.df <- readRDSw('cnjpop.df.noout.rds')
head(cnjpop.pheno.df)
str(cnjpop.pheno.df)

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

#Combine means into one dataframe
cnjpop.pheno.means.df <- cbind(cnjpop.pheno.berry_length.means.df, berry_width=cnjpop.pheno.berry_width.means.df[,4], berry_weight=cnjpop.pheno.berry_weight.means.df[,4], num_seeds=cnjpop.pheno.num_seeds.means.df[,4])
colnames(cnjpop.pheno.means.df) <- c("population", "year", "accession_name", "berry_length", "berry_width", "berry_weight", "num_seeds")

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

#Convert population, year, and accession_name to factors
#cnjpop.pheno.df <- convert_to_factors(cnjpop.pheno.df)
#cnjpop.pheno.means.df <- convert_to_factors(cnjpop.pheno.means.df)
#cnjpop.pheno.p1.df <- convert_to_factors(cnjpop.pheno.p1.df)
#cnjpop.pheno.p1.means.df <- convert_to_factors(cnjpop.pheno.p1.means.df)
#cnjpop.pheno.p2.df <- convert_to_factors(cnjpop.pheno.p2.df)
#cnjpop.pheno.p2.means.df <- convert_to_factors(cnjpop.pheno.p2.means.df)

#Now generate all pairwise comparisons for population 1
#Rather than build up a weird list and generate all pairwise combinations, use the dataframe merge
#function to build up a dataframe organized by genotypes (accession names) and each column has a year.variable
#name

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

focal.cols <- c("berry_length", "berry_width", "berry_weight", "num_seeds")

#Do a combined population assessment
cnjpop.pheno.combined.means.df <- rbind(cnjpop.pheno.p1.means.df,cnjpop.pheno.p2.means.df) #With incomplete genotypes removed
cnjpop.pheno.combined.means.cor.mat <- cor(cnjpop.pheno.combined.means.df[,focal.cols])
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

postscriptw(file="p12_phenotypes.wholeCor.eps", title="Phenotypic Correlations")
levelplot(cnjpop.pheno.combined.means.cor.mat,scales=list(x=list(rot=90)),main='Combined Populations', xlab='',ylab='',col.regions=colorRampPalette(c('white','black'))(100))
levelplot(cnjpop.pheno.p1.means.cor.mat,scales=list(x=list(rot=90)),main='CNJ04', xlab='',ylab='',col.regions=colorRampPalette(c('white','black'))(100))
levelplot(cnjpop.pheno.p1.means.cor.split.mat,scales=list(x=list(rot=90)),main='CNJ04: Year Split', xlab='',ylab='',col.regions=colorRampPalette(c('white','black'))(100))
levelplot(cnjpop.pheno.p2.means.cor.mat,scales=list(x=list(rot=90)),main='CNJ02', xlab='',ylab='',col.regions=colorRampPalette(c('white','black'))(100))
levelplot(cnjpop.pheno.p2.means.cor.split.mat,scales=list(x=list(rot=90)),main='CNJ02: Year Split', xlab='',ylab='',col.regions=colorRampPalette(c('white','black'))(100))
dev.off()

#From these correlation plots, the following is obvious:
# - Berry width, length, and weight are all highly correlated, while the number of seeds seems less correlated.
#     - Interestingly, width is more highly correlated to weight and so is length to weight, but width and length aren't necessarily as highly correlated as these
#     - The number of seeds seems to be most correlated to width for both populations, although it is still low overall.
# - The year to year correlation within some phenotypes is not as high as expected, indicating that year does not necessarily provide information on the value of a given phenotype.
#     - Exceptions: The berry length seems to have the highest year to year correlation (predicted more by the environment?)
#
# Luis pointed out that CNJ04 (smaller pop) has very low year-year correlation w/in traits, so it is not a good population to do assessment on as heritability is low.  He thus said
# that it is probably better to work on CNJ02 (higher within-trait correlations across years), and to first do multivariate analysis on highly correlated phenotypic variables
# (berry weight, length, and width) within years (years are independent) to allow them to correct for each other, then to see how year-year correlation stacks up after BLUPs calculated.
# 

#Luis mentioned that the only reason we really want to look at correlations at the phenotypic level is if we want to model them
#in multivariate analysis together such that they will correct for each other.
#
#The one catch to multivariate analysis is that if you want to get separate QTL peaks for the two separate variables that are being modeled as a response together,
#they will show the same peaks instead of being separate.
#
#However, if you want to find separate QTLs for the different traits, you need to model one trait in the covariates and treat the other trait as a response.


#
#He also talked about why we care to model as random effects: We need the variance components for heritability, and we need to be able to test that variance > 0.

#Luis says that it is best to use 'lattice' package and levelplot() function.
#cnjpop.pheno.df$population <- as.numeric(cnjpop.pheno.df$population)
#cnjpop.pheno.df.cormat <- round(cor(cnjpop.pheno.df[,c("year","population","berry_length","berry_width","berry_weight","num_seeds")],use='complete.obs'),3)

## separating data by population and year
#cnjpop.pheno.mat <- matrix(data=cnjpop.pheno.df,nrow=length(levels(cnjpop.pheno.df$population)),ncol=length(levels(cnjpop.pheno.df$year)))

#Begin assessing the data - Only assess population 2 first (CNJ02*)

#Read in the marker maps
p1.geno<-read.table(geno_rpath2fpath("CNJ04_AllASMapData.csv"),header=T,sep=',')  
p2.geno<-read.table(geno_rpath2fpath("CNJ02_AllASMapData.csv"),header=T,sep=',')  
rownames(p1.geno)<-p1.geno$X
rownames(p2.geno)<-p2.geno$X
p1.geno<-p1.geno[,-c(1:6)] #Remove the marker name, segregation pattern, phase, classification, position, and lg fields -- only keep the genotype calls
p2.geno<-p2.geno[,-c(1:6)] #Remove the marker name, segregation pattern, phase, classification, position, and lg fields -- only keep the genotype calls
p1.geno.num <-atcg1234(t(p1.geno))
p2.geno.num <-atcg1234(t(p2.geno))

#A note here: Luis calculated his additive genetic effects using the means of his phenotypic values.  I'd like to both try the means and the actual values (will need to duplicate rows for each duplicated genotype in realized additive matrix)
#and see what differences I can find in heritability of one versus the other.
#NOTE: Only doing means first, as I noticed that trying to do all data will take a tremendous amount of memory and computing resources on my computer.
#

mixed_model_analyze <- function(pheno, geno, header.cols, tag=NA, additional.model=NA) {
    analyses.l <- vector("list", 2)
    #Random genotype effects incidence matrix (per sommer documentation)
    Zg      <- model.matrix(~accession_name-1,pheno); colnames(Zg) <- gsub("accession_name","",colnames(Zg))
    A       <- A.mat(geno)
    if( !is.na(additional.model) ) {
        ETA.A <- additional.model
        ETA.A$accession_name <- list(Z=Zg,K=A)
    } else {
        ETA.A   <- list(accession_name=list(Z=Zg,K=A)) 
    }

    header              <- header.cols
    multivariate.traits <- c("berry_length","berry_width","berry_weight")
    repo    <- mmer(Y=pheno[,multivariate.traits],Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=FALSE)
    vcov    <- repo$var.comp
    h2      <- (diag(length(multivariate.traits))*(vcov$accession_name/Reduce('+',vcov))) %*% c(1,1,1)
    analyses.l[[1]] <- list(pheno=pheno[,c(header,multivariate.traits)], geno=geno, A=A, blups=repo$u.hat$accession_name, vcov=vcov, h2=h2, traits=multivariate.traits, repo=repo)
    univariate.trait <- "num_seeds"
    repo    <- mmer(Y=pheno[,univariate.trait],Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=FALSE)
    blups   <- repo$u.hat$accession_name
    colnames(blups) <- univariate.trait #Necessary for subsequent indexing to work on single-column datasets -- consider emailing Eduardo about this.
    vcov    <- repo$var.comp
    h2      <- vcov['Var(accession_name)','component']/sum(vcov[,'component'])
    #gcor     <- cor(as.matrix(blups),use='complete.obs')
    #rownames(gcor) <- rownames(geno)
    #colnames(gcor) <- colnames(geno)
    #l$univariate <- list(pheno=pheno[,c(header,univariate.trait)], geno=geno, A=A, blups=blups, vcov=vcov, h2=h2, gcor=gcor)
    analyses.l[[2]] <- list(pheno=pheno[,c(header,univariate.trait)], geno=geno, A=A, blups=blups, vcov=vcov, h2=h2, traits=univariate.trait)

    if( !is.na(tag) ) {
        l <- list(description=tag, analyses=analyses.l)
    } else {
        l <- list(analyses=analyses.l)
    }

    return(l)
}

#Separate the phenotypes by year
#Overwriting cnjpop.pheno.p2.df (was filled in with means)
intersect.geno.pheno <- intersect(unique(cnjpop.pheno.p2.means.df[,"accession_name"]), rownames(p2.geno.num))

cnjpop.mmer.p2 <- vector("list",length(unique((cnjpop.pheno.p2.means.df$year)))+1) #Add one for all years
idx <- 1
for (year in unique(cnjpop.pheno.p2.means.df$year)) {
    #Technically, since I'm analyzing the years separately here, I don't need to use 'complete genotypes' across all three years.  This is more important when
    #modeling the year as a covariate.  Think about changing this.
    include.idx                   <- which( (year == cnjpop.pheno.p2.means.df$year) & (cnjpop.pheno.p2.means.df$accession_name %in% cnjpop.pheno.p2.CompleteGenoTypes) & (cnjpop.pheno.p2.means.df$accession_name %in% intersect.geno.pheno) )

    #First do multivariate analysis on length, width, and weight, since they have a high level of correlation between them.
    pheno   <- cnjpop.pheno.p2.means.df[include.idx,]
    geno    <- p2.geno.num[as.character(pheno$accession_name),]
    pheno$accession_name <- as.character(pheno$accession_name)

    header <- c("year", "accession_name")
    cnjpop.mmer.p2[[idx]] <- mixed_model_analyze(pheno, geno, header, tag=year)

    idx <- idx + 1
}

#Now do a all-years analysis of traits
include.idx                   <- which( (cnjpop.pheno.p2.means.df$accession_name %in% cnjpop.pheno.p2.CompleteGenoTypes) & (cnjpop.pheno.p2.means.df$accession_name %in% intersect.geno.pheno) )
pheno   <- cnjpop.pheno.p2.means.df[include.idx,]
geno    <- p2.geno.num[unique(pheno$accession_name),]
header <- c("year", "accession_name")
#The following only generates the correct matrix if years are converted to factors first.
pheno$year <- as.factor(pheno$year)
Zy <- model.matrix(~year-1,pheno); colnames(Zy) <- gsub("year","",colnames(Zy))
years <- levels(pheno$year)
Ky <- diag(length(years))
colnames(Ky) <- years
rownames(Ky) <- years
cnjpop.mmer.p2[[idx]] <- mixed_model_analyze(pheno, geno, header, tag="all-years", additional.model=list(year=list(Z=Zy,K=Ky)))

#Save the mixed model analysis results to an R-data file for later analysis
saveRDSw(cnjpop.mmer.p2,'cnjpop.mmer2.p2.rds', compress=T)

#Read the mixed model analysis results from the R-data file for later analysis -- Start from here if we want to save time
cnjpop.mmer.p2 <- readRDSw('cnjpop.mmer2.p2.rds')


#Convert consensus map to bins
superMap.df<-read.table(geno_rpath2fpath('consensusMapAll2.csv'),header=T,sep=',')
superMap.df<-superMap.df[,c('marker','LG','consensus')]
superMap.df<-superMap.df[order(superMap.df[,2],superMap.df[,3]),]
superMap.df$binID<-NA
superMap.bin.df<-numeric()
for (LG in unique(superMap.df$LG)){
  f<-which(superMap.df$LG==LG)
  mybins<-unique(superMap.df$consensus[f])
  mybinsID<-paste0('bin_',LG,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap.df$consensus[f]==mybins[j])
    superMap.df$binID[f[f2]]<-mybinsID[j]      
    superMap.bin.df<-rbind(superMap.bin.df,superMap.df[f[f2[1]],])
  }
}
rownames(superMap.bin.df)<-superMap.bin.df$marker


## genetic analysis
geno<-read.table(geno_rpath2fpath('CNJ02_AllASMapData.csv'),header=T,sep=',')  


matrixK<-matrix(NA,nrow=nrow(geno),ncol=ncol(geno)-6)
g<-which(geno[,7:ncol(geno)]=='ac')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='ad')
matrixK[g]<-3  
g<-which(geno[,7:ncol(geno)]=='bc')
matrixK[g]<-2
g<-which(geno[,7:ncol(geno)]=='bd')
matrixK[g]<-4
g<-which(geno[,7:ncol(geno)]=='ll')
matrixK[g]<-5
g<-which(geno[,7:ncol(geno)]=='lm')
matrixK[g]<-6
g<-which(geno[,7:ncol(geno)]=='nn')
matrixK[g]<-7
g<-which(geno[,7:ncol(geno)]=='np')
matrixK[g]<-8
g<-which(geno[,7:ncol(geno)]=='hh')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='hk')
matrixK[g]<-10
g<-which(geno[,7:ncol(geno)]=='kk')
matrixK[g]<-4

matrixK[1:10,1:10]
gData<-t(matrixK)
rownames(gData)<-colnames(geno)[7:ncol(geno)]
colnames(gData)<-geno$X
gData[1:10,1:10]


#Generate the QTL cross files for reading into R/qtl
for (mmer in cnjpop.mmer.p2) {
    print(paste0("Population Analysis Set: ",mmer$description))
    for(analysis in mmer$analyses) {
        traits <- analysis$traits
        geno.intersect<-intersect(rownames(analysis$blups),rownames(gData))
        y<-as.data.frame(analysis$blups[geno.intersect,traits]) #Convert to data.frame to deal with issues in setting colnames in univariate analysis (blups aren't a data.frame in this case)
        colnames(y)<-traits

        gData.sub <- gData[geno.intersect,]
        geno.intersect.sub<-intersect(colnames(gData.sub),superMap.bin.df$marker)
        gData.sub<-gData.sub[,geno.intersect.sub]
        superMap.sub<-superMap.bin.df[geno.intersect.sub,]

        gData.sub<-rbind(superMap.sub$LG,superMap.sub$consensus,gData.sub)

        colnames(gData.sub)<-superMap.sub$binID

        gData.sub<-data.frame(rbind('','',y),gData.sub)
        rownames(gData.sub)[1:2]<-c('chr','pos')

        qtlfile=paste0('qtl/',mmer$description,"__",paste(traits,collapse="__"),"__QTL.csv")
        write.csv(file=geno_dpath2fpath(qtlfile),gData.sub,row.names = FALSE)
    }
}

#Read in the cross files for R/qtl and then calculate the genotype probabilities
for (i in 1:length(cnjpop.mmer.p2)) {
    mmer <- cnjpop.mmer.p2[[i]]
    for(j in 1:length(mmer)) {
        analysis <- mmer$analyses[[j]]
        #Have to index in main structure as R doesn't keep lists as references but instead makes copies!
        traits <- analysis$traits
        cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls <- vector("list",length(traits))
        qtlfile=paste0('qtl/',mmer$description,"__",paste(analysis$traits,collapse="__"),"__QTL.csv")
        cnjpop.cross <- (read.cross(format = "csv", file=geno_dpath2fpath(qtlfile),genotypes = NULL))
        cnjpop.mmer.p2[[i]]$analyses[[j]]$cross <- calc.genoprob(cnjpop.cross,step=0,map.function="kosambi") 
        cnjpop.mmer.p2[[i]]$analyses[[j]]$scanone.out.ehk <- scanone(cnjpop.mmer.p2[[i]]$analyses[[j]]$cross,method='ehk',pheno.col=1:length(traits))
        #To get statistical thresholds
        #operms <- scanone(cnjpop.mmer.p2[[i]]$analyses[[j]]$cross, n.perm=1000, verbose=FALSE)
        #cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.hk.sum <- summary(cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.hk, perms=operms, alpha=0.1, pvalues=TRUE)
        #cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.ehk.sum <- summary(cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.ehk, perms=operms, alpha=0.1, pvalues=TRUE)
    }
}

#Generate the plots -- assume symmetrical data
col.years.pal <- brewer.pal(3, 'Dark2')
mmer <- cnjpop.mmer.p2[[1]]
for( j in 1:length(mmer$analyses) ) {
    analysis <- mmer$analyses[[j]]
    quartzw(file=paste0("p13_scanone_QTLS_",paste(analysis$traits,collapse="__"),".pdf"), title="Whole Genome QTL Plots")
    old.par <- par(mfrow=c(2,length(analysis$traits))) #2 rows for 2 graphs per trait: 1 for all years combined, 1 for all-years in model
    scanones.year.l <- list(cnjpop.mmer.p2[[1]]$analyses[[j]]$scanone.out.ehk,cnjpop.mmer.p2[[2]]$analyses[[j]]$scanone.out.ehk,cnjpop.mmer.p2[[3]]$analyses[[j]]$scanone.out.ehk)
    years.v <- c(cnjpop.mmer.p2[[1]]$description, cnjpop.mmer.p2[[2]]$description, cnjpop.mmer.p2[[3]]$description)
    #First combine plots over all years
    for( k in 1:length(analysis$traits) ) {
        trait <- analysis$traits[[k]]
        do.call.list <- scanones.year.l
        do.call.list$main <- paste0("Trait ",analysis$traits[[k]],collapse="")
        do.call.list$col <-  col.years.pal[1:3]
        do.call.list$lodcolumn <- k
        do.call("plot", do.call.list)
        legend("topleft", legend=years.v, fill=col.years.pal[1:3], col=col.years.pal[1:3])
    }
    #Then show all-years plot.
    for( k in 1:length(analysis$traits) ) {
        trait <- analysis$traits[[k]]
        plot(cnjpop.mmer.p2[[4]]$analyses[[j]]$scanone.out.ehk, main=paste0("Trait (", cnjpop.mmer.p2[[4]]$description, ") ", analysis$traits[[k]],collapse=""),lodcolumn=k)
    }
    par(old.par)
    dev.off()
}

for (i in 1:length(cnjpop.mmer.p2)) {
    for(j in 1:length(mmer)) {
        analysis <- mmer$analyses[[j]]
        #Have to index in main structure as R doesn't keep lists as references but instead makes copies!
        traits <- analysis$traits
        cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls <- vector("list",length(traits))
        qtlfile=paste0('qtl/',mmer$description,"__",paste(analysis$traits,collapse="__"),"__QTL.csv")
        cnjpop.cross <- (read.cross(format = "csv", file=geno_dpath2fpath(qtlfile),genotypes = NULL))
        cnjpop.mmer.p2[[i]]$analyses[[j]]$cross <- calc.genoprob(cnjpop.cross,step=0,map.function="kosambi") 
        cnjpop.mmer.p2[[i]]$analyses[[j]]$scanone.out.ehk <- scanone(cnjpop.mmer.p2[[i]]$analyses[[j]]$cross,method='ehk',pheno.col=1:length(traits))
        #To get statistical thresholds
        #operms <- scanone(cnjpop.mmer.p2[[i]]$analyses[[j]]$cross, n.perm=1000, verbose=FALSE)
        #cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.hk.sum <- summary(cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.hk, perms=operms, alpha=0.1, pvalues=TRUE)
        #cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.ehk.sum <- summary(cnjpop.mmer.p2[[i]]$analyses[[j]]$trait_qtls[[k]]$scanone.out.ehk, perms=operms, alpha=0.1, pvalues=TRUE)
    }
}

cnjpop.mmer.p2[[1]]$analyses[[1]]$trait_qtls[[1]]$scanone.out.ehk
#Resave the mixed model analysis results to an R-data file for later analysis
saveRDSw(cnjpop.mmer.p2,'cnjpop.mmer2.p2.rds', compress=T)

#Print chromosome 2, 10, and 11 QTL plots for trait berry_length, width, and weight (num_seeds looked like noise on output)
postscriptw(file="p14_scanone_chr2_10_11_length_width_weight_QTLS.eps", title="Chromosome 2, 10, and 11 for berry_length, width, and weight QTL Plots")
for (i in 1:length(cnjpop.mmer.p2)) {
    mmer <- cnjpop.mmer.p2[[i]]
    analysis <- mmer$analyses[[1]]
    traits <- analysis$traits
    trait_qtls <- analysis$trait_qtls
    traits.num <- length(traits)
    for( k in 1:traits.num ) {
        plot(trait_qtls[[k]]$scanone.out.hk, trait_qtls[[k]]$scanone.out.ehk, chr=c(2,10,11), main=paste0(mmer$description, ": HK Regression on Trait ", traits[[k]]), col=c("blue","red"))
        legend("top", legend=c("HK", "EHK"), fill=c("blue", "red"), col=c("blue","red"))
    }
}
dev.off()

sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties,additive.only=FALSE)
  return(x1)
}

#For vaccininum server, set cores to 32
n.cores <- 32

cl <- snow::makeCluster(n.cores, spec="MPI")
registerDoSNOW(cl)
for (i in 1:5){
    for (j in 1:2){
        superOut <- foreach(ph = 1:ncol(phenoData[[i]][[j]]$cross$pheno),.packages = 'qtl') %dopar% perroQTL(phenoData[[i]][[j]]$cross,ph,max.qtl=10,penalties=sixo)
        phenoData[[i]][[j]]$scantwo<-superOut
    }
    print(i)
}

snow::stopCluster(cl)

for (tri in 1:5){
    for (ii in 1:2){
        modelRes<-vector("list", 4) 
        for (tri2 in 1:ncol(phenoData[[tri]][[ii]]$cross$pheno)){
            rob1<-phenoData[[tri]][[ii]]$scantwo[[tri2]]
            if (length(rob1)==0){

            }else{
                superQTL<-makeqtl(phenoData[[tri]][[ii]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
                md1<-fitqtl(phenoData[[tri]][[ii]]$cross,pheno.col = tri2,superQTL,formula = formula(rob1),get.ests=F)
                x<-summary(md1)
                modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                        drop=as.data.frame(x$result.drop),
                        model=rob1,
                        inter=numeric())
                for (j in 1:length(rob1$chr)){
                    mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
                    modelRes[[tri2]]$inter[j]<-mylod
                }
            }
        } 
        phenoData[[tri]][[ii]]$qtldata<-modelRes
    }
}

# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

jaguar1<-data.frame(trait=numeric(),year=numeric(),date=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:5){
  for (i2 in 1:2){
  t1<-mytraits[i]
  t2<-as.character(mymonths[i2])
  for (j in 1:ncol(phenoData[[i]][[i2]]$cross$pheno)){
    t3<-as.character(myyears[j])
    nQTL<-length(phenoData[[i]][[i2]]$qtldata[[j]]$model$chr)
    if (nQTL==0){

      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]][[i2]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]][[i2]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]][[i2]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]][[i2]]$qtldata[[j]]$drop))
      if (nQTL==1){
        xK<-phenoData[[i]][[i2]]$qtldata[[j]]$variance
      }else{
      x2<-phenoData[[i]][[i2]]$qtldata[[j]]$drop[f2,]
      xK<-x2$'%var'
    }

      jaguar1[cont,]<-c(t1,t3,t2,mychr,mypos,mymarker,xK,phenoData[[i]][[i2]]$qtldata[[j]]$variance,phenoData[[i]][[i2]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)

jaguar1$date2<-paste0(jaguar1$date,jaguar1$year)
mydates<-as.character(unique(jaguar1$date2))


# circos


# ideogram
ideo1<-numeric()
for (i in 1:12){
  ideo1<-c(ideo1,paste0('chr - c',i,' chr',i,' 0 ',max(superMap3$consensus[which(superMap3$LG==i)])*1000," vvdgrey"))
}
write(ideo1, file = '/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogram.txt')




mycols<-matrix(paste0(hueGen(30,0,340),'_a2'),nrow=6,byrow=F)
mycols2<-matrix(hueGen(30,0,340),nrow=6,byrow=F)


for (ii in 1:5){

mark1<-matrix(NA,nrow=length(which(jaguar1$trait==mytraits[ii])),ncol=1)
cont<-1
for (i in 1:nrow(jaguar1)){
  if (jaguar1$trait[i]==mytraits[ii]){
   mark1[cont,]<-paste0('c',jaguar1$chr[i]," ",round(jaguar1$position[i]*1000,0),' ',round(jaguar1$position[i]*1000,0)+1,' ',which(jaguar1$date2[i]==mydates),' color=',mycols[which(jaguar1$date2[i]==mydates),ii],',stroke_color=',mycols2[which(jaguar1$date2[i]==mydates),ii],',glyph_size=',round(jaguar1$marker.variance[i]*3,0))
  cont<-cont+1
}
}

filename<-paste0('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat',ii,'.txt')
write(mark1, file = filename)

mark2<-matrix(NA,nrow=length(which(jaguar1$trait==mytraits[ii])),ncol=1)

cont<-1
for (i in 1:nrow(jaguar1)){
  if (jaguar1$trait[i]==mytraits[ii]){
   minInt<-round(jaguar1$position[i]+(1000*(jaguar1$position[i]-.5*jaguar1$interval[i])),0)
   maxInt<-round(jaguar1$position[i]+(1000*(jaguar1$position[i]+.5*jaguar1$interval[i])),0)
   mark2[cont,]<-paste0('c',jaguar1$chr[i]," ",minInt,' ',maxInt,' ',which(jaguar1$date2[i]==mydates),' color=',mycols2[which(jaguar1$date2[i]==mydates),ii])
   cont<-cont+1
}
}


filename<-paste0('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line',ii,'.txt')
write(mark2, file = filename)


}


back1<-c('<backgrounds>',
          '<background>',
          'color=vdgrey',
          'y0=0',
          '</background>',
          '<background>',
          'y0=3.5',
          'y1=7',
          'color=dgrey',
          '</background>',
          '</backgrounds>')

back2<-c('<backgrounds>',
          '<background>',
          'y0=0',
          'y1=7',
          'color=vdgrey',
          '</background>',
          '</backgrounds>')



rmin<-seq(.3,.85,length.out=5)
rmax<-rmin-(.01-diff(rmin[1:2]))

rmin1<-rmin
rmax1<-rmin+.7*diff(rmin[1:2])

rmin2<-rmax1+0.01
rmax2<-rmax


rmin1<-round(rmin1,3)
rmin2<-round(rmin2,3)
rmax1<-round(rmax1,3)
rmax1<-round(rmax1,3)

circosStarter<-c('<<include etc/colors_fonts_patterns.conf>>',
                  '<<include /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ticks.conf>>',
                  'karyotype = /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogram.txt',
                  'chromosomes_units=1000',
                  'chromosomes_display_default=yes',
                  '<<include /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogramLalo.conf>>',
                  '<image>',
                  'radius* = 1500p',
                  '<<include etc/image.conf>>',
                  '</image>',
                  'chromosomes_radius=c1:1r,c2:1r,c3:1r,c4:1r,c5:1r,c6:1r,c7:1r,c8:1r,c9:1r,c10:1r,c11:1r,c12:1r',
                  'chromosomes_order=c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12',
                  '<<include etc/housekeeping.conf>>',
                  '<plots>',
                  

                  '<plot>',
                  paste0('r0=',rmin1[1],'r'),
                  paste0('r1=',rmax1[1],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat1.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line1.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[1],'r'),
                  paste0('r1=',rmax2[1],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                  '<plot>',
                  paste0('r0=',rmin1[2],'r'),
                  paste0('r1=',rmax1[2],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat2.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line2.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[2],'r'),
                  paste0('r1=',rmax2[2],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                                    '<plot>',
                  paste0('r0=',rmin1[3],'r'),
                  paste0('r1=',rmax1[3],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat3.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line3.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[3],'r'),
                  paste0('r1=',rmax2[3],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                  '<plot>',
                  paste0('r0=',rmin1[4],'r'),
                  paste0('r1=',rmax1[4],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat4.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line4.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[4],'r'),
                  paste0('r1=',rmax2[4],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',


                  '<plot>',
                  paste0('r0=',rmin1[5],'r'),
                  paste0('r1=',rmax1[5],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat5.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line5.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[5],'r'),
                  paste0('r1=',rmax2[5],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',




                  '</plots>')


write(circosStarter, file = '/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/circos.conf')


cmd <- paste("perl /Users/luisdiaz/Documents/chamba/BIOINFO_TOOLS/circos-0.67-7/bin/circos -conf ~/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/circos.conf")
system(cmd)




saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataCNJ02.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1CNJ02.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3CNJ02.rds')


## a better exploration of color

tmp1<-scantwo(phenoData[[2]]$sep$cross,pheno.col=1,chr=3)
plot(clean(tmp1))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[2]]$oct$cross,mname1=find.marker(phenoData[[2]]$oct$cross,3,i),draw=F)
print(i)
}

mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)







####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################'
####################################################################################################################################################################################'
############################################################################     P I C T U R E S      ##############################################################################'
####################################################################################################################################################################################'
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES





digitalCN<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markCNJ_WISC.csv',header=TRUE,sep=',')



mypops<-gregexpr('CNJ02',digitalCN$id3)
tt<-sapply(mypops, "[[", 1)

f02<-which(tt==1)
f04<-which(tt==-1)


mlx<-lm(digitalCN$g.mean[f02]~digitalCN$id3[f02])
yy02<-unique(predict(mlx))
yyRes02<-tapply(digitalCN$g.mean[f02],digitalCN$id3[f02],sd)
yyRes02<-yyRes02[1:166]

names(yy02)<-names(yyRes02)
cor(yy02,yyRes02,use='complete.obs')

digData02<-data.frame(color=yy02,colovar=yyRes02)
f<-which(digData02[,1]>0.51)
digData02<-digData02[-f,]
gt1<-gsub(pattern = "-",replacement = "_",rownames(digData02))
gt2<-gsub(pattern = ".JPG",replacement = "",gt1)
rownames(digData02)<-gt2


mlx<-lm(digitalCN$g.mean[f04]~digitalCN$id3[f04])
yy04<-unique(predict(mlx))
yyRes04<-tapply(digitalCN$g.mean[f04],digitalCN$id3[f04],sd)
yyRes04<-yyRes04[167:235]

names(yy04)<-names(yyRes04)
cor(yy04,yyRes04,use='complete.obs')
digData04<-data.frame(color=yy04,colovar=yyRes04)
gt1<-gsub(pattern = "-",replacement = "_",rownames(digData04))
gt2<-gsub(pattern = ".JPG",replacement = "",gt1)
rownames(digData04)<-gt2



superMap<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/consensusMapAll2.csv',header=T,sep=',')
superMap<-superMap[,c('marker','LG','consensus')]
superMap<-superMap[order(superMap[,2],superMap[,3]),]
superMap$binID<-NA
superMap2<-numeric()
for (i in 1:12){
  f<-which(superMap$LG==i)
  mybins<-unique(superMap$consensus[f])
  mybinsID<-paste0('bin_',i,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap$consensus[f]==mybins[j])
    superMap$binID[f[f2]]<-mybinsID[j]      
    superMap2<-rbind(superMap2,superMap[f[f2[1]],])
  }
}

rownames(superMap2)<-superMap2$marker



# cnj02
geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/CNJ02_AllASMapData.csv',header=T,sep=',')  



matrixK<-matrix(NA,nrow=nrow(geno),ncol=ncol(geno)-6)
g<-which(geno[,7:ncol(geno)]=='ac')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='ad')
matrixK[g]<-3  
g<-which(geno[,7:ncol(geno)]=='bc')
matrixK[g]<-2
g<-which(geno[,7:ncol(geno)]=='bd')
matrixK[g]<-4
g<-which(geno[,7:ncol(geno)]=='ll')
matrixK[g]<-5
g<-which(geno[,7:ncol(geno)]=='lm')
matrixK[g]<-6
g<-which(geno[,7:ncol(geno)]=='nn')
matrixK[g]<-7
g<-which(geno[,7:ncol(geno)]=='np')
matrixK[g]<-8
g<-which(geno[,7:ncol(geno)]=='hh')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='hk')
matrixK[g]<-10
g<-which(geno[,7:ncol(geno)]=='kk')
matrixK[g]<-4

matrixK[1:10,1:10]
gData<-t(matrixK)
rownames(gData)<-colnames(geno)[7:ncol(geno)]
colnames(gData)<-geno$X
gData[1:10,1:10]



f<-intersect(rownames(digData02),rownames(gData))
y2<-digData02[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]


write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))



phenoData<-list()
phenoData[[1]]<-list(cross=calc.genoprob(superF,step=0,map.function="kosambi"),blups=digData02)

phenoData[[1]]$blups<-phenoData[[1]]$blups*-1
#### cnj04

geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/CNJ04_AllASMapData.csv',header=T,sep=',')  



matrixK<-matrix(NA,nrow=nrow(geno),ncol=ncol(geno)-6)
g<-which(geno[,7:ncol(geno)]=='ac')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='ad')
matrixK[g]<-3  
g<-which(geno[,7:ncol(geno)]=='bc')
matrixK[g]<-2
g<-which(geno[,7:ncol(geno)]=='bd')
matrixK[g]<-4
g<-which(geno[,7:ncol(geno)]=='ll')
matrixK[g]<-5
g<-which(geno[,7:ncol(geno)]=='lm')
matrixK[g]<-6
g<-which(geno[,7:ncol(geno)]=='nn')
matrixK[g]<-7
g<-which(geno[,7:ncol(geno)]=='np')
matrixK[g]<-8
g<-which(geno[,7:ncol(geno)]=='hh')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='hk')
matrixK[g]<-10
g<-which(geno[,7:ncol(geno)]=='kk')
matrixK[g]<-4

matrixK[1:10,1:10]
gData<-t(matrixK)
rownames(gData)<-colnames(geno)[7:ncol(geno)]
colnames(gData)<-geno$X
gData[1:10,1:10]



f<-intersect(rownames(digData04),rownames(gData))
y2<-digData04[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]


write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))



phenoData[[2]]<-list(cross=calc.genoprob(superF,step=0,map.function="kosambi"),blups=digData04)


phenoData[[2]]$blups<-phenoData[[2]]$blups*-1





sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties)
  return(x1)
}



n.cores <- 6

cl <- snow::makeCluster(n.cores)
registerDoSNOW(cl)
for (i in 1:2){
superOut <- foreach(ph = 1:2,.packages = 'qtl') %dopar% perroQTL(phenoData[[i]]$cross,ph,max.qtl=5,penalties=sixo)
phenoData[[i]]$scantwo<-superOut
print(i)
}

snow::stopCluster(cl)







for (tri in 1:2){
  modelRes<-vector("list", 4) 
  for (tri2 in 1:2){
    rob1<-phenoData[[tri]]$scantwo[[tri2]]
    if (length(rob1)==0){
      
    }else{
      superQTL<-makeqtl(phenoData[[tri]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
      md1<-fitqtl(phenoData[[tri]]$cross,pheno.col = tri,superQTL,formula = formula(rob1),get.ests=F)
      x<-summary(md1)
      modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                      drop=as.data.frame(x$result.drop),
                      model=rob1,
                      inter=numeric())
      for (j in 1:length(rob1$chr)){
        mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
        modelRes[[tri2]]$inter[j]<-mylod
      }
    }
  }
  phenoData[[tri]]$qtldata<-modelRes
}


phenoData[[1]]$trait<-'cnj02'
phenoData[[2]]$trait<-'cnj04'

# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

myt1<-colnames(phenoData[[1]]$blups)
jaguar1<-data.frame(population=numeric(),trait=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:2){
  t1<-phenoData[[i]]$trait
  for (j in 1:2){
    nQTL<-length(phenoData[[i]]$qtldata[[j]]$model$chr)
    if (nQTL==1){
     jaguar1[cont,]<-c(t1,myt1[j],as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr),phenoData[[i]]$qtldata[[j]]$model$pos,paste0(x$binID[idx],'_',x$marker[idx]),rep(phenoData[[i]]$qtldata[[j]]$variance,2),phenoData[[i]]$qtldata[[j]]$inter)
      cont<-cont+1
      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]]$qtldata[[j]]$drop))
      x2<-phenoData[[i]]$qtldata[[j]]$drop[f2,]
      jaguar1[cont,]<-c(t1,myt1[j],mychr,mypos,mymarker,x2$'%var',phenoData[[i]]$qtldata[[j]]$variance,phenoData[[i]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)


saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_bothPop.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1Pictures_bothPop.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3_bothPop.rds')




## a better exploration of color

tmp1<-scantwo(phenoData[[2]]$cross,pheno.col=1,chr=3)
plot(clean(tmp1))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[2]]$cross,mname1=find.marker(phenoData[[2]]$cross,3,i),draw=F)
print(i)
}

mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(0.35,.46),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)





##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES


GRYG1<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2014.csv',header=TRUE,sep=',')
GRYG1$year<-'a'
GRYG2<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2015.csv',header=TRUE,sep=',')
GRYG2$year<-'b'
colnames(GRYG2)<-colnames(GRYG1)
GRYG3<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2016.csv',header=TRUE,sep=',')
GRYG3$year<-'c'
colnames(GRYG3)<-colnames(GRYG1)



GRYG1$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG1$id3)))
GRYG2$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG2$id3)))
GRYG3$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG3$id3)))

phenoData<-list()

phenoData[[1]]<-list(rawData=GRYG1)
phenoData[[2]]<-list(rawData=GRYG2)
phenoData[[3]]<-list(rawData=GRYG3)





for (i in 1:3){
  mdl1<-lmer(b.mean.1~(1|id3),data=phenoData[[i]]$rawData)
  blup1<-ranef(mdl1)$id3
  myColVar<-tapply(phenoData[[i]]$rawData$b.mean.1,list(myID=phenoData[[i]]$rawData$id3),sd)
  f<-intersect(names(myColVar),rownames(myColVar))
  blupsX<-as.data.frame(cbind(blup1[f,],myColVar[f]))
  colnames(blupsX)<-c('color','colorvar')
  phenoData[[i]]$blups<-blupsX
}

pairs(phenoData[[3]]$blups)

f<-intersect(rownames(phenoData[[1]]$blups),rownames(phenoData[[2]]$blups))
f<-intersect(f,rownames(phenoData[[3]]$blups))



for (i in 1:3){
  phenoData[[i]]$blups2<-phenoData[[i]]$blups[f,]
  rownames(phenoData[[i]]$blups2)<-paste0('P',rownames(phenoData[[i]]$blups2))
}


# in the following code I will combine years into a single list element. 

phenoDataX<-phenoData


# doing multivariate




geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/Gryg_AllASMapData.csv',header=T,sep=',')  
rownames(geno)<-geno$X
geno<-geno[,-c(1:6)]
geno2<-atcg1234(t(geno))


phenoData<-list()
mytraits<-c('color','colorvar')
myyears<-c('2014','2015','2016')
for (i in 1:2){

phenoA<-cbind(phenoDataX[[1]]$blups2[,i],phenoDataX[[2]]$blups2[,i],phenoDataX[[3]]$blups2[,i])

rownames(phenoA)<-rownames(phenoDataX[[1]]$blups2)


f<-intersect(rownames(phenoA),rownames(geno2))
phenoA2<-phenoA[f,]
geno3<-geno2[f,]
A<-A.mat(geno3)


Za <- diag(nrow(A)) 
ETA.A <- list(add=list(Z=Za,K=A)) 
repo1<-mmer(Y=phenoA2,Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=TRUE)
blups1<-repo1$u.hat$add
var1<-repo1$var.comp


phenoData[[i]]<-list(blups=as.matrix(blups1),vcov=var1,trait=mytraits[i])

}


for (i in 1:2){
  phenoData[[i]]$blups<-phenoData[[i]]$blups*-1
}


# calculating h2 and correlations

for (i in 1:2){
  
    gv<-diag(phenoData[[i]]$vcov$add)
    rv<-diag(phenoData[[i]]$vcov$Residual)
    phenoData[[i]]$h2<-gv/(gv+rv)
    x<-cor(phenoData[[i]]$blups,use='complete.obs')
    phenoData[[i]]$gcor<-x
    colnames(phenoData[[i]]$gcor)<-myyears
    rownames(phenoData[[i]]$gcor)<-myyears
    names(phenoData[[i]]$h2)<-myyears
    rownames(phenoData[[i]]$blups)<-rownames(phenoA2)
    colnames(phenoData[[i]]$blups)<-myyears
}




superMap<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/consensusMapAll2.csv',header=T,sep=',')
superMap<-superMap[,c('marker','LG','consensus')]
superMap<-superMap[order(superMap[,2],superMap[,3]),]
superMap$binID<-NA
superMap2<-numeric()
for (i in 1:12){
  f<-which(superMap$LG==i)
  mybins<-unique(superMap$consensus[f])
  mybinsID<-paste0('bin_',i,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap$consensus[f]==mybins[j])
    superMap$binID[f[f2]]<-mybinsID[j]      
    superMap2<-rbind(superMap2,superMap[f[f2[1]],])
  }
}

rownames(superMap2)<-superMap2$marker

## genetic analysis
geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/Gryg_AllASMapData.csv',header=T,sep=',')  

matrixK<-matrix(NA,nrow=nrow(geno),ncol=ncol(geno)-6)
g<-which(geno[,7:ncol(geno)]=='ac')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='ad')
matrixK[g]<-3  
g<-which(geno[,7:ncol(geno)]=='bc')
matrixK[g]<-2
g<-which(geno[,7:ncol(geno)]=='bd')
matrixK[g]<-4
g<-which(geno[,7:ncol(geno)]=='ll')
matrixK[g]<-5
g<-which(geno[,7:ncol(geno)]=='lm')
matrixK[g]<-6
g<-which(geno[,7:ncol(geno)]=='nn')
matrixK[g]<-7
g<-which(geno[,7:ncol(geno)]=='np')
matrixK[g]<-8
g<-which(geno[,7:ncol(geno)]=='hh')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='hk')
matrixK[g]<-10
g<-which(geno[,7:ncol(geno)]=='kk')
matrixK[g]<-4

matrixK[1:10,1:10]
gData<-t(matrixK)
rownames(gData)<-colnames(geno)[7:ncol(geno)]
colnames(gData)<-geno$X
gData[1:10,1:10]

for (ii in 1:2){


f<-intersect(rownames(phenoData[[ii]]$blups),rownames(gData))
y2<-phenoData[[ii]]$blups[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]

#write.csv(file='/crdata4/luis4/QTLcolorServer_copy/testQTL',gData5,row.names = FALSE)
write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))


phenoData[[ii]]$cross<-numeric()
phenoData[[ii]]$cross <- calc.genoprob(superF,step=0,map.function="kosambi") 

}




sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties)
  return(x1)
}



n.cores <- 6

cl <- snow::makeCluster(n.cores)
registerDoSNOW(cl)
for (i in 1:2){
superOut <- foreach(ph = 1:3,.packages = 'qtl') %dopar% perroQTL(phenoData[[i]]$cross,ph,max.qtl=10,penalties=sixo)
phenoData[[i]]$scantwo<-superOut
print(i)
}

snow::stopCluster(cl)





for (tri in 1:2){
  modelRes<-vector("list", 4) 
  for (tri2 in 1:3){
    rob1<-phenoData[[tri]]$scantwo[[tri2]]
    if (length(rob1)==0){
      
    }else{
      superQTL<-makeqtl(phenoData[[tri]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
      md1<-fitqtl(phenoData[[tri]]$cross,pheno.col = tri2,superQTL,formula = formula(rob1),get.ests=F)
      x<-summary(md1)
      modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                      drop=as.data.frame(x$result.drop),
                      model=rob1,
                      inter=numeric())
      for (j in 1:length(rob1$chr)){
        mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
        modelRes[[tri2]]$inter[j]<-mylod
      }
    }
  }
  phenoData[[tri]]$qtldata<-modelRes
}



# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

jaguar1<-data.frame(trait=numeric(),date=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:2){
  t1<-phenoData[[i]]$trait
  for (j in 1:3){
    nQTL<-length(phenoData[[i]]$qtldata[[j]]$model$chr)
    if (nQTL==1){
     jaguar1[cont,]<-c(t1,myt1[j],as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr),phenoData[[i]]$qtldata[[j]]$model$pos,paste0(x$binID[idx],'_',x$marker[idx]),rep(phenoData[[i]]$qtldata[[j]]$variance,2),phenoData[[i]]$qtldata[[j]]$inter)
      cont<-cont+1
      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]]$qtldata[[j]]$drop))
      x2<-phenoData[[i]]$qtldata[[j]]$drop[f2,]
      jaguar1[cont,]<-c(t1,myyears[j],mychr,mypos,mymarker,x2$'%var',phenoData[[i]]$qtldata[[j]]$variance,phenoData[[i]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)


saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_GRYG_3years.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1Pictures_GRYG_3years.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3_GRYG_3years.rds')



## a better exploration of color

phenoData<-readRDS('~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_GRYG_3years.rds')

tmp1<-scantwo(phenoData[[1]]$cross,pheno.col=1,chr=3)

postscript(file="~/Documents/chamba/PHENOTYPING/QTLcolorPaper/scantwogryg.eps", 
            width=10, 
            height=10, 
            horizontal=TRUE) 
plot(clean(tmp1))
dev.off()


plot(scanone(phenoData[[1]]$cross,pheno.col=1,chr=3))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[1]]$cross,mname1=find.marker(phenoData[[1]]$cross,3,i),draw=F)
print(i)
}
postscript(file="~/Documents/chamba/PHENOTYPING/QTLcolorPaper/effectplotGRYG.eps", 
            width=20, 
            height=10, 
            horizontal=TRUE) 
mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(-.03,.03),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)

dev.off()




#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts


# loading libraries
source('./usefulFunctions.R')


workflow <- get0("workflow", ifnotfound="../../Workflows/1")
P1_Name <- get0("P1_Name", ifnotfound="Mullica_Queen")
P2_Name <- get0("P2_Name", ifnotfound="Crimson_Queen")

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

library(sommer)
library(qtl)
library(dplyr)
library(rlist)
library(tidyverse)

#Load the configuration file for this workflow
source(paste0(workflow,"/configs/model.cfg"))

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
geno<-read.csv(geno_rpath2fpath(geno_file),header=T)
rownames(geno)<-geno$X

generateParentalMarker <- function(segregations,parent) {
	subexpr <- paste0("gsub('<([ablmn]{2})x([cdlnp]{2})>','",'\\\\',parent,"',segregations,fixed=FALSE)")
	return(eval(parse(text=subexpr)))
}

#Calculate and add-in the parental markers for calculating relationship matrix and generating parental BLUP values
geno.p <- geno
geno.p[,P1_Name] = generateParentalMarker(geno.p$Segregation,parent=1)
geno.p[,P2_Name] = generateParentalMarker(geno.p$Segregation,parent=2)
rownames(geno.p) <- geno.p$X

geno <- geno.p %>% dplyr::select(-c('X','Segregation','Phase','Classification','Position','LG'))
#geno.num is input to sommer's A.mat() function in order to calculate the realized additive matrix, aiding in calculation of breeding values (BLUPS)
#
#

#Make a copy of geno for passing to atcg1234() and ultimately A.mat()
#NOTE: I'm removing this b/c I think I did this in error -- The goal is to convert ab x cd bi-allelic SNPs to an hk x hk form, thus allowing atcg1234() to do its work.
geno.amat <- as.matrix(geno)

#Always filter scaffolds here, as this section is for calculating the Additive relationship matrix for calculating BLUPS, which is distinct
#from the encoding for the R/qtl cross file.
#geno.scaffolds.idx <- which(grepl("^scaffold_", rownames(geno.amat), ignore.case=T))
#Not sure I can really convert ab x cd segregation patterns into hk-style segregation - although I'm confused, as I thought GBS outputs are biallelic!
#The reason I doubt myself in converting abxcd to hkxhk is b/c you can't map hk x hk segregation-type markers -- you don't know which parent the h or k came from, so you can't definitively
#know in child if it is a recombination event or not
#ac.i <- which(geno.amat[geno.scaffolds.idx,] == "ac")
#bc.i <- which(geno.amat[geno.scaffolds.idx,] == "bc")
#ad.i <- which(geno.amat[geno.scaffolds.idx,] == "ad")
#bd.i <- which(geno.amat[geno.scaffolds.idx,] == "bd")
#geno.amat[geno.scaffolds.idx,][ac.i] <- "hh"
#geno.amat[geno.scaffolds.idx,][bc.i] <- "hk"
#geno.amat[geno.scaffolds.idx,][ad.i] <- "hk"
#geno.amat[geno.scaffolds.idx,][bd.i] <- "kk"

geno.num <- atcg1234(t(geno.amat))$M

#Read in the phenotype file
pheno.means.df <- read.csv(file=pheno_dpath2fpath(pheno_file))
pheno.means.df$rowf <- as.factor(pheno.means.df$row) #Needed for modeling row effects
pheno.means.df$columnf <- as.factor(pheno.means.df$column) #Needed for modeling column effects
pheno.means.df$year <- as.factor(paste0("20",as.character(pheno.means.df$year)))
#pheno.means.df$year <- as.factor(pheno.means.df$year) #Needed for modeling column effects

superMap.df       <- read.csv(geno_rpath2fpath(geno_consensus_file),header=T)[,c("marker","LG","consensus")] %>%
                        mutate(binID = generate_bin_id(LG,consensus))
rownames(superMap.df) <- superMap.df$marker
superMap.bin.df   <- condense_map_bins(superMap.df)
rownames(superMap.bin.df) <- superMap.bin.df$marker
saveRDS(superMap.bin.df, file=geno_rpath2fpath(paste0(geno_consensus_file,".rds")), compress=T)

#The following is necessary to convert the mapqtl codes in a four-way cross to those that are defined
#according to the R/qtl read.cross() specifications.
#Remove parents, as they should not be in the cross file
geno.no.p <- geno %>% dplyr::select(-c(P1_Name,P2_Name))
rownames(geno.no.p) <- rownames(geno)
matrixK<-matrix(NA,nrow=nrow(geno.no.p),ncol=ncol(geno.no.p))
g<-which(geno.no.p=='ac')
matrixK[g]<-1
g<-which(geno.no.p=='ad')
matrixK[g]<-3  
g<-which(geno.no.p=='bc')
matrixK[g]<-2
g<-which(geno.no.p=='bd')
matrixK[g]<-4
g<-which(geno.no.p=='ll')
matrixK[g]<-5
g<-which(geno.no.p=='lm')
matrixK[g]<-6
g<-which(geno.no.p=='nn')
matrixK[g]<-7
g<-which(geno.no.p=='np')
matrixK[g]<-8
g<-which(geno.no.p=='hh')
matrixK[g]<-1
g<-which(geno.no.p=='hk')
matrixK[g]<-10
g<-which(geno.no.p=='kk')
matrixK[g]<-4

gData<-t(matrixK)
rownames(gData)<-colnames(geno.no.p)
colnames(gData)<-rownames(geno.no.p)
saveRDS(gData, file=paste0(workflow,"/traits/rqtl.gdata.rds"), compress=TRUE)

#Need to rethink this, as I'm certain I'm calculating the heritability incorrectly in more complicated analyses
generate_h2_formula <- function(sigma, genotype_id_string) {
    num_var_components <- nrow(sigma)
    genotype_var_idx   <- which(grepl(genotype_id_string, rownames(sigma)))
    h2_denom           <- vector("character",num_var_components)
    for( i in 1:nrow(sigma) ) {
        h2_denom[i]    <- paste0("V",i)
    }
    h2_string          <- paste0("h2 ~ V",genotype_var_idx,"/(",paste0(h2_denom,collapse=" + "),")")
    return(h2_string)
}


#A note here: Luis calculated his additive genetic effects using the means of his phenotypic values.  I'd like to both try the means and the actual values (will need to duplicate rows for each duplicated genotype in realized additive matrix)
#and see what differences I can find in heritability of one versus the other.
#NOTE: Only doing means first, as I noticed that trying to do all data will take a tremendous amount of memory and computing resources on my computer.
#
mixed_model_analyze <- function(trait.cfg, pheno, geno) {
    pheno$id    <- as.factor(pheno$accession_name)
    geno        <- geno[unique(pheno$accession_name),]
    trait       <- trait.cfg$trait
    A           <- A.mat(geno)

    random <- as.formula(trait.cfg$random)
    yuyur <- strsplit(as.character(random[2]), split = "[+-]")[[1]]
    randomtermss <- apply(data.frame(yuyur),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    
    #Generate a grid of all possible random model terms, but always include first term -- Define the first term as always being the genotype
    randomtermss_combos = expand.grid(c(TRUE,rep(list(c(FALSE,TRUE)),length(randomtermss)-1)))
    AICs                = rep(Inf, nrow(randomtermss_combos)) # All the AIC terms -- the model with minimum AIC is the winner
    LogLikelihoods      = rep(-Inf, nrow(randomtermss_combos)) # All the LogLikelihoods terms
    randomfs            = rep("", nrow(randomtermss_combos))
    models              = list()
    print(paste0("Running trait \"", trait, "\" and model \"",trait.cfg$model,"\" analysis."))
    for( i in 1:nrow(randomtermss_combos) ) {
        randomtermss_include = unlist(randomtermss_combos[i,])
        randomtermss_current = randomtermss[randomtermss_include]
        randomf              = paste("~",paste(randomtermss_current,collapse = " + ")) #Convert back to formula string
        randomfs[i]          = randomf
        if( !all(is.na(pheno[,trait.cfg$label_short])) ) {
            mmer.expr <- paste0(c("mmer(fixed=", trait.cfg$fixed, ", random=", randomf, ", rcov=", trait.cfg$rcov, ", data=pheno"), collapse="")
            if( !is.empty(trait.cfg["mmer_args"]) ) {
                mmer.expr <- paste0(mmer.expr,",",trait.cfg["mmer_args"])
            }
            mmer.expr <- paste0(mmer.expr,", date.warning=FALSE)")
            print(mmer.expr)
            model <- eval(parse(text=mmer.expr))
            if( !is.null(model) && !is.null(model$AIC) && !is.null(model$monitor) ) {
                AICs[i]             = model$AIC
                LogLikelihoods[i]   = model$monitor[1,ncol(model$monitor)]
            }
            models = list.append(models,model)
        }
    }
    
    min_index = which.min(AICs)
    if( AICs[min_index] != Inf ) #Inf indicates no model could be fit
    {
        model       = models[[min_index]]
        model$A.mat = A #Store the additive relationship matrix so it does not need to be calculated later when doing a drop1 style anova to compare two models.
    } else {
        model = list() #Empty list indicates optimal model could not be found and fit
    }
    return(list(bestmodel=model, stats=tibble(randomfs=randomfs, AICs=AICs, LogLikelihoods=LogLikelihoods)))
}

#Read in the model trait configuration file to determine how to model traits.
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"))
for( i in 1:length(traits.df[,1]) ) {
    pheno.filtered.df    <- pheno.means.df %>% filter(accession_name %in% rownames(geno.num))
    trait.cfg <- traits.df[i,]
    if ( trait_is_unmasked(trait.cfg) ) {
        #Make a trait folder if it doesn't exist and add output files, such as blups, heritability, and cross files to it.
        trait            <- trait.cfg$trait
        trait_subfolder  <- paste0(c(trait.cfg$model,trait),collapse="--")
        trait_subfolder_fpath <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        dir.create(trait_subfolder_fpath, showWarnings = FALSE)
        if( !is.empty(trait.cfg$filter) ) {
            pheno.filtered.df <- eval(parse(text=paste0("pheno.filtered.df %>% filter(", trait.cfg$filter, ")")))
        }
        response    <- mixed_model_analyze(trait.cfg, pheno.filtered.df, geno.num)
        write_csv(response$stats, file=paste0(trait_subfolder_fpath,"/all_fitted_model_stats.csv"))
        model       <- response$bestmodel
        if( !is.null(model) && !is.null(model$U) ) { #Indicates that the model was singular, and failed to find a solution
            model.sum <- summary(model)
            saveRDS(model, file=paste0(trait_subfolder_fpath,"/mmer.rds"), compress=TRUE)
            blups   <- model$U[["u:id"]][[trait]]
            vcov    <- model.sum$varcomp
            randomtermss <- trimws(strsplit(as.character(model$call$random)[[2]],split='[+]')[[1]])
            rcovtermss <- trimws(strsplit(as.character(model$call$rcov)[[2]],split='[+]')[[1]])
            rownames(vcov) <- c(randomtermss,rcovtermss)
            write.csv(vcov,file=paste0(trait_subfolder_fpath,"/vcov.csv"))
        } else {
            cat(paste0("Best model for trait ", trait, " is singular or failed to find a solution.\n"))
            saveRDS(model, file=paste0(trait_subfolder_fpath,"/mmer.rds"), compress=TRUE)
        }
    }
}

save.image(paste0(workflow,"/.RData.01_genmodeleeffects"))

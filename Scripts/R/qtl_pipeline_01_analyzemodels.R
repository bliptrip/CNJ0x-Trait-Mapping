#year_var=numeric(), This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#!/usr/bin/env RScript

#Analyzes previously fitted models by generating anova tables using a drop1 term likelihood ratio approach to determine significance of an effect,
#and based on whether g:y interaction effects are significant or not (where applicable), chooses a better model and generates relevant cross
#files for r/qtl using selected model BLUPs.
#


# loading libraries
source('./usefulFunctions.R')

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

library(ggplot2)
library(qtl)
library(rlist)
library(sommer)
library(tidyverse)

source(paste0(workflow,'/configs/model.cfg'))

#Loop through previously generated models and gather summary stats on them, plot residuals, etc.  Consider comparing across years and also for all-years.
#Things to assess/plot:
# 1) Model significances.
# 2) Total model variances for BLUPs, residuals, and other covariates.
# 3) Plot of residuals versus index to see if residuals are homo- or heteroscedastic
# 4) Compare residuals across three years.
# 5) Look at year-to-year correlation of BLUPs.  Are they tight for 1st-order differences?  What about 2nd-order differences?
# 6) How about spearman-rank correlations.  What is the makeup of this compared to Pearson's correlation (should be similar).
# 7) Are the years significantly different for the three years?  What about pairwise comparisons of years?


readModelsCB  <- function(trait.cfg, trait.path, models.l) {
    trait               <- trait.cfg$trait
    model.map.l         <- models.l$model_map 
    model.collated.df   <- models.l$model_collated_table
    model               <- readRDS(file=paste0(trait.path,"/mmer.rds"))
    cat(paste0("\nCalculating pvalues for Trait: ",trait.cfg$trait," Model: ",trait.cfg$model,"\n"), fill=FALSE)
    append.pointer(model.map.l, model)
    model_idx           <- length(model.map.l$value)
    var.comp <- summary(model)$varcomp
    var.i <- grep("u:id", rownames(var.comp)) #Return variance component with genotype identifier
    units.i <- grep("units", rownames(var.comp))
    append.pointer(model.collated.df, c(trait.cfg$model, trait.cfg$trait, model_idx, var.comp$VarComp[var.i], var.comp$VarComp[units.i]))
    #Interesting components of model
    # $var.comp[$year, accession_name, ...]
    # $u.hat - Estimated BLUPs for RE covariates (years & genotypic values)
    # $Var.u.hat - Variance estimates for RE covariates (years & genotypic values)
    # $fitted.y
    # $fitted.u
    # $LL
    # $AIC
    # $BIC
    # $ZETA - the original Z-matrices specified for RE (year & genotypic values)
    #AMNOTE: Why is genomic relationship matrix of individual w/ itself not 1?  It looks like this might be intrinsic to the way additive relationship matrices
    # are calculated (XX'/c).
    
    
    #Use a 'drop1' type approach to determine significance of model terms
    #Get our A matrix, since it will be necessary.
    A <- model$A.mat
    #Gather fixed terms
    fixed <- model$call$fixed
    yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
    response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    fixedtermss <- apply(data.frame(yuyuf),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fixedtermss <- fixedtermss[which(fixedtermss != "-1")]
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
    #Gather random terms
    random <- model$call$random
    yuyur <- strsplit(as.character(random)[2], split = "[+-]")[[1]]
    randomtermss <- apply(data.frame(yuyur),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    randomtermss    <- randomtermss[which(randomtermss != "-1")]
    randomtermss    <- randomtermss[which(randomtermss != "1")]
    rcov <- model$call$rcov
    yuyuu <- strsplit(as.character(rcov[2]), split = "[+-]")[[1]]
    rcovtermss <- apply(data.frame(yuyuu),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    allterms          <- c(fixedtermss,randomtermss,rcovtermss)
    allterms.size     <- length(allterms)
    anova.combined    <- tibble(dropterm=character(allterms.size), #Term that is dropped from the model
                                Chisq=character(allterms.size),
                                ChiDf=character(allterms.size),
                                PrChisq=numeric(allterms.size),
                                PrChisqInfo=character(allterms.size))
    #rownames(anova.combined) <- as.character(allterms)
    #Now for each fixed term, drop it from model and calculate prob under assumption of chisqr distro for likelihood ratio
    i = 1
    #Now for each random term, drop it from model and calculate prob under assumption of chisqr distro for likelihood ratio
    k = 1
    for (j in seq_along(randomtermss)) {
        usef <- setdiff(randomtermss,randomtermss[k])
        randomf <- paste("~",paste(usef,collapse = " + "))
        if( length(usef) > 0 ) {
                mmer.expr <- paste0(c("mmer(fixed=", fixed, ", random=", randomf, ", rcov=", rcov, ", data=model$dataOriginal"), collapse="")
        } else {
                mmer.expr <- paste0(c("mmer(fixed=", fixed, ", rcov=", rcov, ", data=model$dataOriginal"), collapse="")
        }
        if( !is.empty(trait.cfg["mmer_args"]) ) {
                mmer.expr <- paste0(mmer.expr,",",trait.cfg["mmer_args"])
        }
        mmer.expr <- paste0(mmer.expr,", date.warning=FALSE)")
        cat(mmer.expr)
        model.drop1 <- try(eval(parse(text=mmer.expr)))
        if( !is.null(model.drop1) && !is_bare_list(model.drop1) ) { #Indicates that the model was singular, and failed to find a solution
            anova.cols    <- c("Chisq", "ChiDf", "PrChisq", "PrChisqInfo")
            anova.drop1 <- anova(model, model.drop1)[2,]
            prchisq <- strsplit(anova.drop1$PrChisq, split = " ")[[1]]
            anova.drop1$PrChisq <- as.numeric(prchisq[1])
            anova.drop1$PrChisqInfo <- prchisq[2]
            anova.combined[k+(i-1),anova.cols] = anova.drop1[,anova.cols]
            anova.combined[k+(i-1),"dropterm"] = randomtermss[k]
            k = k + 1
        }
    }    
    vcov <- read_csv(file=paste0(trait.path,'/vcov.csv'),show_col_types=FALSE)
    vcov$PrNorm <- pnorm(abs(vcov$Zratio),lower.tail=FALSE)
    vcov$PrNormInfo <- unlist(map(vcov$PrNorm,siginfo))
    vcov$dropterm <- rownames(vcov)
    vcov <- as_tibble(vcov)
    vcov.cols <- c("dropterm","Zratio","PrNorm","PrNormInfo")
    anova.combined <- anova.combined %>%
                        filter(dropterm != "") %>%
                        mutate(PrChisqInfo = ifelse(is.na(PrChisqInfo),"NS",PrChisqInfo)) %>%
                        left_join(vcov %>% dplyr::select(vcov.cols), by="dropterm")
    #Where applicable, add Zratio values for random effects (as found in vcov matrix) and pnorm values for these
    saveRDS(anova.combined, file=paste0(trait.path,"/anova.rds"), compress=TRUE)
    write.csv(anova.combined, file=paste0(trait.path,"/anova.csv"), row.names=TRUE)
}

model.map.l.p       <- newPointer(list())
model.collated.df   <- data.frame(model=character(),traits=character(),model_idx=integer(), bv_var=numeric(), res_var=numeric(), stringsAsFactors=FALSE)
model.collated.df.p <- newPointer(model.collated.df)
loopThruTraits(workflow, readModelsCB, list(model_map=model.map.l.p, model_collated_table=model.collated.df.p))

#Now that we've generated anova tables, choose the optimal model for all-years based on significance values of any g:y effects -- if not significant, then choose
#the model without this effect, as I've found that sommer's mmer() function will give invalid g:y variance values (negative values) under some situations, and if I remove
#this from the modeling term, the variance terms look more reasonable
#Use the selected model to generate relevant cross file for QTL analysis and BLUP plot generation/etc.
selectModelCB <- function(trait.cfg, trait.path, selectedModels) {
    trait          <- trait.cfg$trait
    model_label      <- trait.cfg$model
    cat(paste0("\nChecking significance of g:y interaction effects: ",trait.cfg$trait," Model: ",trait.cfg$model,"\n"), fill=FALSE)
    anova.drop1.selectedmodel.label <- NA
    anova.drop1.selectedmodel.label <- "full" #Default to using full model
    if( model_label == 'all-years' ) { #Interaction effects only apply for 'all-years' model - year-by-year model doesn't have an interaction effect
        anova.file <- paste0(trait.path,"/anova.csv")
        if( file.exists(anova.file) ) {
            anova.int    <- read_csv(anova.file, show_col_types=FALSE) %>% filter(dropterm == "id:year") #Filter out row with g:y interaction effects
            if( (nrow(anova.int) > 0) && (anova.int$PrChisq > int_effects_alpha) ) { #g:y interaction effect is not significant -- choose model that drops this term
                anova.drop1.selectedmodel.label <- "id:year"
            }
        }
    }
    append.pointer(selectedModels, c(model_label,trait,anova.drop1.selectedmodel.label))
}
models.selected.df   <- data.frame(model=character(),trait=character(),selected_model=character())
models.selected.df.p <- newPointer(models.selected.df)
loopThruTraits(workflow, selectModelCB, models.selected.df.p)
models.selected.df <- models.selected.df.p$value
write.csv(models.selected.df, file=paste0(workflow,"/traits/selectedModels.csv"), row.names=F)

generateCrossFilesCB <- function(trait.cfg,trait.path,args.l) {
    gData             <- args.l[[1]]
    map.df             <- args.l[[2]]
    trait_name         <- trait.cfg$label_short
    model_label         <- trait.cfg$model
    model   <- readRDS(file=paste0(trait.path,"/mmer.rds"))
    blups   <- model$U[["u:id"]][[trait_name]]
    cat(paste0("\nGenerating cross file for: ", trait.cfg$trait, " Model: ", trait.cfg$model, "\n"), fill=FALSE)
    generate_cross_file(trait.cfg, blups, gData, map.df, trait.path)
}

gData.rqtl        <- readRDS(file=paste0(workflow,"/traits/rqtl.gdata.rds"))
superMap.df       <- read_csv(geno_rpath2fpath(geno_consensus_file), show_col_types=FALSE)[,c("marker","LG","consensus")] %>%
                        mutate(binID = generate_bin_id(LG,consensus))
rownames(superMap.df) <- superMap.df$marker
loopThruTraits(workflow,generateCrossFilesCB,list(gData.rqtl,superMap.df,models.selected.df.p))

save.image(paste0(workflow,"/.RData.01_analyzemodels"))
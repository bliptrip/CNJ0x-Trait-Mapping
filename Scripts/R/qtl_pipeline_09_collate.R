#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
library(jsonlite)
library(qtl)
library(tidyverse)

# loading libraries
source('./usefulFunctions.R')

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}


source(paste0(workflow,"/configs/model.cfg"))

#Read in the consensus map file, modified in qtl_pipeline_01_genmodels.R
supermap.bin.df <- readRDS(file=geno_rpath2fpath(paste0(geno_consensus_file,".rds")))

#Function to deconstruct a QTL formula into it's parts
deconstruct_qtl_formula <- function(qtl_formula) {
    qtls <- trimws(unlist(strsplit(gsub("y[ ]*~[ ]*(.+)","\\1",qtl_formula),"\\+")))
    return(qtls)
}

#Function that iterates through a set of model qtls and fills in the collated datatable
generate_qtl_collate <- function(cross, qtl, qtl_model, method, model, trait, loopArgs) {
    collated.df.p <- loopArgs$qtls

    qtl.model.terms  <- deconstruct_qtl_formula(formula(qtl))
    qtl.num          <- length(qtl.model.terms)
    #Now add in anova p-values for genotype effects (BLUPs), and if applicable, for GxE
	cat(paste0("Collated QTLs for trait: ",trait,", model: ", model, ", method: ",method,",", "formula: ",formula(qtl),".\n"))
    anova.df <- read.csv(file=paste0(workflow,'/traits/',model,'--',trait,'/anova.csv'), header=TRUE, row.names=2)
    GLRpvalue    <- anova.df['vs(id, Gu = A)','PrChisq']
    GZRpvalue    <- anova.df['vs(id, Gu = A)','PrNorm']
    if( model == "all-years" ) {
        GxYLRpvalue <- anova.df['id:year','PrChisq']
        GxYZRpvalue <- anova.df['id:year','PrNorm']
    } else {
        GxYLRpvalue <- NA
        GxYZRpvalue <- NA
    }
    for (qtl.model.term in qtl.model.terms) {
        #Initialize some of the entries to NA for each step.
        mychr2        <- NA
        mypos2        <- NA
        mymarker      <- NA
        qtl_lodint    <- NA

        #Is this an QTL:QTL pairwise interaction or just a single additive QTL effect
        qtl.model.term.v   <- unlist(strsplit(qtl.model.term,":"))
        qtl_idx1           <- which(qtl$altname == qtl.model.term.v[1])
        mychr              <- as.numeric(qtl$chr[qtl_idx1])
        mypos              <- qtl$pos[qtl_idx1]
        #single-additive effect
        if( length(qtl.model.term.v) == 1 ) {
            qtl_lodint    <- diff(lodint(qtl,qtl.index=qtl_idx1)[c(1,3),2])
            f             <- which(supermap.bin.df$LG==mychr)
            x             <- supermap.bin.df[f,]
            idx           <- which.min(abs(x$consensus-mypos))
            mymarker      <- paste0(x$binID[idx],'_',x$marker[idx])
            f2            <- intersect(qtl$name[qtl_idx1],rownames(qtl_model$result.drop))
            qtl_model.var <- qtl_model$result.full["Model","%var"]
            if (qtl.num==1) {
                qtl.var   <- qtl_model$result.full["Model","%var"]
                qtl.lod   <- qtl_model$result.full["Model","LOD"]
                qtl.p     <- qtl_model$result.full["Model","Pvalue(F)"]
            } else {
                qtl.var   <- qtl_model$result.drop[f2,"%var"]
                qtl.lod   <- qtl_model$result.drop[f2,"LOD"]
                qtl.p     <- qtl_model$result.drop[f2,"Pvalue(F)"]
            }
        } else { #pairwise effect
            qtl_idx2      <- which(qtl$altname == qtl.model.term.v[2])
            mychr2        <- as.numeric(qtl$chr[qtl_idx2])
            mypos2        <- qtl$pos[qtl_idx2]
            qtl_model.var <- qtl_model$result.full["Model","%var"]
            qtl_name      <- paste0(qtl$name[qtl_idx1],":",qtl$name[qtl_idx2]);
            f2            <- intersect(qtl_name,rownames(qtl_model$result.drop))
            if (qtl.num==1) {
                qtl.var   <- qtl_model$result.full["Model","%var"]
                qtl.lod   <- qtl_model$result.full["Model","LOD"]
                qtl.p     <- qtl_model$result.full["Model","Pvalue(F)"]
            } else {
                qtl.var   <- qtl_model$result.drop[f2,"%var"]
                qtl.lod   <- qtl_model$result.drop[f2,"LOD"]
                qtl.p     <- qtl_model$result.drop[f2,"Pvalue(F)"]
            }
        }
        #TODO: Fill in interaction qtls
        append.pointer(collated.df.p, c(method, model, trait, mychr, mypos, mychr2, mypos2, mymarker, qtl.lod, qtl.var, qtl.p, qtl_model.var, qtl_lodint, GLRpvalue, GxYLRpvalue, GZRpvalue, GxYZRpvalue))
    }
}


#Generate a single dataframe with qtl information in it, and save to a csv file.
#index into qtl.collated.df
qtl.collated.df   <- data.frame(method=character(),model=character(),trait=character(),chr=numeric(),position=numeric(),chr2=numeric(),position2=numeric(),nearest.marker=numeric(),qtl.lod=numeric(),marker.variance=numeric(),qtl.pvalue=numeric(),model.variance=numeric(),interval=numeric(),GLRpvalue=numeric(),GxYLRpvalue=numeric(),GZRpvalue=numeric(),GxYZRpvalue=numeric(),stringsAsFactors=FALSE)
qtl.collated.df.p <- newPointer(qtl.collated.df)

collateQtlCB      <- function(trait.cfg, trait.path, loopArgs) {
    model <- as.character(trait.cfg$model)
    #Read in the cross file
    cross_path = paste0(trait.path,'/cross.rds')
    if( file.exists(cross_path) ) {
        cross <- readRDS(file=cross_path)
        trait <- trait.cfg$trait
        if( any(cross$pheno[trait] != 0) ) {
            trait_subsubfolder_fpath <- file.path(trait.path, trait)
            #scanone-derived results
            sofile <- paste0(trait_subsubfolder_fpath,'/scanone.qtl.rds')
            somdfile <- paste0(trait_subsubfolder_fpath,'/scanone.md.rds')
            if( file.exists(sofile) && file.exists(somdfile) ) {
                scan.one.qtl <- readRDS(file=sofile)
                scan.one.md  <- readRDS(file=somdfile)
                generate_qtl_collate(cross, scan.one.qtl, scan.one.md, "scanone", model, trait, loopArgs)
            }
            #stepwiseqtl-derived results
            swfile <- paste0(trait_subsubfolder_fpath,'/scansw.rds')
            swmdfile <- paste0(trait_subsubfolder_fpath,'/scansw.md.rds')
            if( file.exists(swfile) && file.exists(swmdfile) ) {
                scan.sw.qtl <- readRDS(file=swfile)
                scan.sw.md  <- readRDS(file=swmdfile)
                generate_qtl_collate(cross, scan.sw.qtl, scan.sw.md, "stepwiseqtl", model, trait, loopArgs)
            }
        }
    }
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, collateQtlCB, loopArgs=list(qtls=qtl.collated.df.p))

qtl.collated.df <- qtl.collated.df.p$value
write.csv(qtl.collated.df, file=paste0(workflow,'/traits/qtl_collated.csv'), row.names=F)

#Construct a set of names associated with the effect list, just in case the collated file changes order, etc.  This allows
#us to quickly lookup the qtl effects
qtl.collated.names.df <- qtl.collated.df %>% mutate(model_trait=paste(model,trait,sep="--"))
qtl.collated.names.df <- qtl.collated.names.df %>% mutate(id=paste(method,model_trait,trait,chr,position,chr2,position2,sep="/"))

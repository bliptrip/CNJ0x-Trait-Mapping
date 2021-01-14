#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
library(qtl)
library(dplyr)
library(jsonlite)

# loading libraries
source('./usefulFunctions.R')

workflow <- "../../Workflows/1"

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

#Function to strip out uncessary 'bin_<chr>.<pos>cM.' prefix from the effect names.
strip_effect_prefix <- function(effect_names) {
    return(gsub("^.*\\.([[:alpha:]]+)","\\1",effect_names))
}

#Additive-only effects
strip_effect_prefix_1way <- function(effects) {
    names(effects$Means) <- strip_effect_prefix(names(effects$Means))
    names(effects$SEs) <- strip_effect_prefix(names(effects$SEs))
    return(effects)
}

#Pairwise effects -- different structure for the effects
strip_effect_prefix_2way <- function(effects) {
    rownames(effects$Means)  <- strip_effect_prefix(rownames(effects$Means)) 
    colnames(effects$Means)  <- strip_effect_prefix(colnames(effects$Means)) 
    rownames(effects$SEs)  <- strip_effect_prefix(rownames(effects$SEs)) 
    colnames(effects$SEs)  <- strip_effect_prefix(colnames(effects$SEs)) 
    return(effects)
}


#Function that iterates through a set of model qtls and fills in the collated datatable
generate_qtl_collate <- function(cross, qtl, qtl_model, method, model, trait, loopArgs) {
    collated.df.p <- loopArgs$qtls
    effects.l.p   <- loopArgs$effects

    qtl.model.terms  <- deconstruct_qtl_formula(formula(qtl))
    qtl.num          <- length(qtl.model.terms)
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
            #Derive the effects for the QTL
            #m.effects         <- effectplot(cross, pheno.col=trait, mname1=paste0("bin_",mychr,".",mypos,"cM"), draw=FALSE)
            m.effects         <- effectplot(cross, pheno.col=trait, mname1=paste0(mychr,"@",mypos), draw=FALSE)
            #Make the names of the genotypes more manageable to read
            m.effects         <- strip_effect_prefix_1way(m.effects)
            #Convert Means/SEs to list so that the toJSON function will include the name of each effect (associative array)
            m.effects$Means   <- as.list(m.effects$Means)
            m.effects$SEs     <- as.list(m.effects$SEs)
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
            #Interaction effects
            #m.effects         <- effectplot(cross, pheno.col=trait, mname1=paste0("bin_",mychr,".",mypos,"cM"), mname2=paste0("bin_",mychr2,".",mypos2,"cM"), draw=FALSE)
            m.effects         <- effectplot(cross, pheno.col=trait, mname1=paste0(mychr,"@",mypos), mname2=paste0(mychr2,"@",mypos2), draw=FALSE)
            m.effects         <- strip_effect_prefix_2way(m.effects)
            m.effects$Means   <- data.frame(m.effects$Means)
            m.effects$SEs     <- data.frame(m.effects$SEs)
        }
        #TODO: Fill in interaction qtls
        append.pointer(collated.df.p, c(method, model, trait, mychr, mypos, mychr2, mypos2, mymarker, qtl.lod, qtl.var, qtl.p, qtl_model.var, qtl_lodint))
        append.pointer(effects.l.p, m.effects)
    }
}


#Generate a single dataframe with qtl information in it, and save to a csv file.
#index into qtl.collated.df
qtl.collated.df   <- data.frame(method=character(),model=character(),trait=character(),chr=numeric(),position=numeric(),chr2=numeric(),position2=numeric(),nearest.marker=numeric(),qtl.lod=numeric(),marker.variance=numeric(),qtl.pvalue=numeric(),model.variance=numeric(),interval=numeric(),stringsAsFactors=FALSE)
qtl.collated.df.p <- newPointer(qtl.collated.df)
qtl.effects.l.p   <- newPointer(list())

collateQtlCB      <- function(trait.cfg, trait.path, loopArgs) {
    model <- as.character(trait.cfg$model)
    #Read in the cross file
    cross <- readRDS(file=paste0(trait.path,'/cross.rds'))
    trait <- trait.cfg$trait
    trait_subsubfolder_fpath <- file.path(trait.path, trait)
    #scanone-derived results
    scan.one.qtl <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scanone.qtl.rds'))
    scan.one.md  <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scanone.md.rds'))
    generate_qtl_collate(cross, scan.one.qtl, scan.one.md, "scanone", model, trait, loopArgs)
    #stepwiseqtl-derived results
    scan.sw.qtl <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scansw.rds'))
    scan.sw.md  <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scansw.md.rds'))
    generate_qtl_collate(cross, scan.sw.qtl, scan.sw.md, "stepwiseqtl", model, trait, loopArgs)
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, collateQtlCB, loopArgs=list(qtls=qtl.collated.df.p, effects=qtl.effects.l.p))

browser()
qtl.collated.df <- qtl.collated.df.p$value
write.csv(qtl.collated.df, file=paste0(workflow,'/traits/qtl_collated.csv'), row.names=F)

#Construct a set of names associated with the effect list, just in case the collated file changes order, etc.  This allows
#us to quickly lookup the qtl effects
qtl.collated.names.df <- qtl.collated.df %>% mutate(model_trait=paste(model,trait,sep="--"))
qtl.collated.names.df <- qtl.collated.names.df %>% mutate(id=paste(method,model_trait,trait,chr,position,chr2,position2,sep="/"))
qtl.effects.l   <- qtl.effects.l.p$value
names(qtl.effects.l) <- qtl.collated.names.df$id
saveRDS(qtl.effects.l, file=paste0(workflow,'/traits/effects_collated.rds'), compress=T)
qtl.effects.json  <- toJSON(qtl.effects.l, dataframe="columns", auto_unbox=F, pretty=T)
write(qtl.effects.json, paste0(workflow,'/effects_collated.json'))

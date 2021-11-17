#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
#install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

library(dplyr)
library(qtl)
library(jsonlite)

workflow <- get0("workflow", ifnotfound="../../Workflows/1")
#Which circos traits to render.  This file is in similar format to model-traits.cfg.csv.  All it needs is the following columns: trait, label, mask.
#Any mask==TRUE fields means these fields aren't rendered in the plot.
#

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source('./usefulFunctions.R')
source(paste0(workflow,"/configs/model.cfg"))

#Need the supermap file to get marker info.
supermap.bin.df <- readRDS(file=geno_rpath2fpath(paste0(geno_consensus_file,".rds")))

#Extract and refactor LOD profiles for circos plots
exportLODProfile <- function(scan.obj, qtl_type, model, trait) {
    lodprofs         <- attr(scan.obj, "lodprofile")
    lodprofs.organized.l.p <-  newPointer(list())
    num_lod_entries  <- length(unlist(lodprofs))/ncol(lodprofs[[1]])
    model.v          <- vector("character", num_lod_entries)
    trait.v          <- rep(trait, num_lod_entries)
    chr.v            <- vector("numeric", num_lod_entries)
    lod.v            <- vector("numeric", num_lod_entries)
    position.v       <- vector("numeric", num_lod_entries)
    current_index    <- 1
    for( lodprof in lodprofs ) {
        lodprof.len                         <- nrow(lodprof)
        end_index                           <- current_index+lodprof.len-1
        model.v[current_index:end_index]    <- rep(model, nrow(lodprof))
        chr.v[current_index:end_index]      <- as.numeric(as.character(lodprof$chr))
        lod.v[current_index:end_index]      <- lodprof$lod
        position.v[current_index:end_index] <- lodprof$pos
        current_index                       <- end_index + 1
    }
    tibble(chr=chr.v, lod=lod.v, position=position.v) %>%
           group_by(chr, position) %>%
           summarize(lod=max(lod)) %>% #Due to overlapping positions b/w different detected QTLs with different LOD scores (how?), I've decided to take the max lod score of either where they overlap to make the graphs look continuous.
           ungroup() %>%
           group_by(chr) %>%
           group_walk( ~ {
             chr = .y$chr
             assign.pointer(lodprofs.organized.l.p, .x, paste0("vm",chr))
           })
    return(lodprofs.organized.l.p$value)
}

exportLODProfiles <- function(trait, data) {
    lod_profile_models.l.p      <-  newPointer(list())
    data %>% 
      group_by(model) %>%
      group_walk( ~ {
          model                     <- .y$model
          lodprofile_method.l       <- list()
          trait_subfolder           <- paste0(model,"--",trait);
          trait_subfolder_fpath     <- file.path(paste0(workflow,"/traits"), trait_subfolder)
          trait_subsubfolder_fpath  <- file.path(paste0(trait_subfolder_fpath, "/", trait))
          scan.sw                   <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.rds'))
          lodprofile_method.l$stepwiseqtl <- exportLODProfile(scan.sw, "stepwiseqtl", model, trait)
          scan.one.qtl <- readRDS(paste0(trait_subsubfolder_fpath,'/scanone.qtl.rds'))
          lodprofile_method.l$scanone <-  exportLODProfile(scan.one.qtl, "scanone", model, trait)
          assign.pointer(lod_profile_models.l.p, lodprofile_method.l, model)
      })
    return(lod_profile_models.l.p$value)
}

#Consider looping through all QTL's and only keep those that are consistent across all years?
# finding QTL in at least two years
qtl.collated.df <- read.csv(file=paste0(workflow,"/traits/qtl_collated.consensus.csv"),header=T)
#Periods in variable names aren't allowed in javascript, as they have special meaning.  Replace all periods with underscores in the column names.
colnames(qtl.collated.df) <- gsub(".", "_", colnames(qtl.collated.df), fixed=T)
write.csv(qtl.collated.df, file=paste0(workflow,"/",circosfile2path("qtl_collated.consensus.csv")))

lod_profile_traits.l.p <- newPointer(list())
qtl.collated.df %>%
    group_by(trait) %>%
    group_walk( ~ {
                    assign.pointer(lod_profile_traits.l.p, exportLODProfiles(.y$trait, .x), .y$trait)
                  } )

lod_filepath <- paste0(workflow,"/traits/lod_profiles.json")
write_json(lod_profile_traits.l.p$value, lod_filepath, auto_unbox=T, pretty=T)

# Generate Karyotype File for Cranberry using consensus map
karyotype<-numeric()
k.chrs <- 1:12
k.ids <- paste0("vm",k.chrs)
k.labs <-  paste0("LG",k.chrs)
k.cols <- rep("rgb(150,150,150)", length(k.chrs))
k.lens <- (arrange(supermap.bin.df, LG) %>% group_by(LG) %>% summarize(max_consensus=max(consensus)*1000))$max_consensus
karyotype.df <- data.frame(id=k.ids, label=k.labs, color=k.cols, len=k.lens)
karyotype.json <- toJSON(karyotype.df, pretty=T)
write(karyotype.json, file=paste0(workflow,"/",circosfile2path("karyotype.json")))

#Save image for reloading later if desired
save.image(".RData.11_gencircos")

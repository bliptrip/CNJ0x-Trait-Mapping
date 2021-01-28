#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
library(qtl)

args = commandArgs(trailingOnly=TRUE)

workflow="../../Workflows/1"
if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source(paste0(workflow,"/configs/model.cfg"))
source('./usefulFunctions.R')

traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)

#Loop over all mmers, trait groups, subgroups, and perform makeqtl() and fitqtl().
fitQtlCB <- function(trait.cfg, trait.path, funArgs) {
    cross <- readRDS(file=paste0(trait.path,"/cross.rds"))
    trait <- trait.cfg$trait
    if( any(cross$pheno[trait] != 0) ) {
		scan.two.perms <- readRDS(paste0(trait.path, "/operms.2D.rds"))
		trait_subsubfolder_fpath <- file.path(trait.path, trait);
		#scanone build model, fit qtls
		scan.one.qtl <- readRDS(file=paste0(trait_subsubfolder_fpath, "/scanone.qtl.rds"))
		print(paste0("Running fitqtl() (scanone) with model: ",trait.cfg$model, " | trait: ", trait))
		scan.one.md  <- fitqtl(cross, pheno.col=c(trait), qtl=scan.one.qtl, formula = formula(scan.one.qtl), method=qtl_method, get.ests=TRUE)
		saveRDS(scan.one.md, file=paste0(trait_subsubfolder_fpath,'/scanone.md.rds'),compress=T)
		#stepwiseqtl model fitting
		scan.sw <- readRDS(paste0(trait_subsubfolder_fpath, "/scansw.rds"))
		print(paste0("Running fitqtl() (stepwiseqtl) with model: ",trait.cfg$model, " | trait: ", trait))
		scan.sw.md  <- fitqtl(cross, pheno.col=c(trait), qtl=scan.sw, formula = formula(scan.sw), method=qtl_method, get.ests=TRUE)
		saveRDS(scan.sw.md, file=paste0(trait_subsubfolder_fpath,'/scansw.md.rds'),compress=T)
	}
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, fitQtlCB)

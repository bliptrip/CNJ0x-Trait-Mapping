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

#Loop over all mmers, trait groups, subgroups, and perform makeqtl() and fitqtl().
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
for( i in 1:length(traits.df[,1]) ) {
    trait.cfg       <- traits.df[i,]
    if ( is.na(trait.cfg$mask) || (trait.cfg$mask != "TRUE") ) {
        traits           <- unlist(strsplit(trait.cfg$mtraits,","))
        trait.names      <- paste0(traits,collapse="__")
        trait_subfolder  <- paste0(c(trait.cfg$model,trait.names),collapse="--")
        trait_subfolder_fpath <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        cross <- read.cross(format='csv', file=paste0(trait_subfolder_fpath,"/cross.csv"), genotypes=NULL)
        cross <- calc.genoprob(cross, step=0, map.function="kosambi")
        scan.two.perms <- readRDS(paste0(trait_subfolder_fpath, "/operms.2D.rds"))
        for( j in 1:length(traits) ) {
            trait <- traits[j]
            trait_subsubfolder_fpath <- file.path(trait_subfolder_fpath, trait)
            scan.sw <- readRDS(paste0(trait_subsubfolder_fpath, "/scansw.rds"))
            print(paste0("Running makeqtl() and fitqtl() with model: ",trait.cfg$model," | traits: ",trait.names, " | trait: ", trait))
            scan.sw.qtl <- makeqtl(cross, chr=scan.sw$chr, pos=scan.sw$pos, what='prob')
            scan.sw.md  <- fitqtl(cross, pheno.col=c(trait), scan.sw.qtl, formula = formula(scan.sw), get.ests=F)
            scan.sw.summary <-summary(scan.sw.md)
            scan.sw.mdres <-list(variance=scan.sw.summary$result.full[1,5], drop=as.data.frame(scan.sw.summary$result.drop), model=scan.sw, inter=numeric())
            for (k in 1:length(scan.sw$chr)){
                mylod<-diff(lodint(scan.sw,qtl.index=k)[c(1,3),2])
                scan.sw.mdres$inter[k]<-mylod
            }
            saveRDS(scan.sw.qtl, file=paste0(trait_subsubfolder_fpath,'/scansw.qtl.rds'),compress=T)
            saveRDS(scan.sw.md, file=paste0(trait_subsubfolder_fpath,'/scansw.md.rds'),compress=T)
            saveRDS(scan.sw.summary, file=paste0(trait_subsubfolder_fpath,'/scansw.summary.rds'),compress=T)
            saveRDS(scan.sw.mdres , file=paste0(trait_subsubfolder_fpath,'/scansw.mdres.rds'),compress=T)
        }
    }
}



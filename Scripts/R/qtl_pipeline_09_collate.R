#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
library(qtl)

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

#Generate a single dataframe with qtl information in it, and save to a csv file.
qtl.collated.df <-data.frame(model=character(),year=numeric(),mtraits=character(),trait=character(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric(),stringsAsFactors=FALSE)
#index into qtl.collated.df
qtl.i <- 1
#Loop over all mmers, trait groups, subgroups, and perform makeqtl() and fitqtl().
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
for( i in 1:length(traits.df[,1]) ) {
    trait.cfg       <- traits.df[i,]
    if ( is.na(trait.cfg$mask) || (trait.cfg$mask != "TRUE") ) {
        model                    <- as.character(trait.cfg$model)
        print(model)
        year                     <- as.numeric(trait.cfg$year)
        print(year)
        traits                   <- unlist(strsplit(trait.cfg$mtraits,","))
        trait.names              <- paste0(traits,collapse="__")
        print(trait.names)
        trait_subfolder          <- paste0(c(trait.cfg$model,trait.names),collapse="--")
        trait_subfolder_fpath    <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        #Read in the model result file
        for( trait in traits ) {
            trait_subsubfolder_fpath <- file.path(trait_subfolder_fpath, trait)
            scan.sw.mdres  <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scansw.mdres.rds'))
            qtl.num  <- length(scan.sw.mdres$model$chr)
            if( qtl.num > 0 ) {
                for (k in 1:qtl.num) {
                    mypos     <- scan.sw.mdres$model$pos[k]
                    mychr     <- as.numeric(scan.sw.mdres$model$chr[k])
                    f         <- which(supermap.bin.df$LG==mychr)
                    x         <- supermap.bin.df[f,]
                    idx       <- which.min(abs(x$consensus-mypos))
                    mymarker  <- paste0(x$binID[idx],'_',x$marker[idx])
                    f2        <- intersect(scan.sw.mdres$model$name[k],rownames(scan.sw.mdres$drop))
                    if (qtl.num==1) {
                        x.var <- scan.sw.mdres$variance
                    } else {
                        x2    <- scan.sw.mdres$drop[f2,]
                        x.var <- x2$'%var'
                    }
                    qtl.collated.df[qtl.i,] <- c(model, year, trait.names, trait, mychr, mypos, mymarker, x.var, scan.sw.mdres$variance, scan.sw.mdres$inter[k])
                    qtl.i                   <- qtl.i + 1
                }
            }
        }
    }
}

write.csv(qtl.collated.df, file=paste0(workflow,'/traits/qtl_collated.csv'), row.names=F)

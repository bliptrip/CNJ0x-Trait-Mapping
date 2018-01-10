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
#index into qtl.collated.df
qtl.collated.df <- data.frame(model=character(),year=numeric(),mtraits=character(),trait=character(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric(),stringsAsFactors=FALSE)
qtl.collated.df.p <- newPointer(qtl.collated.df)
collate_qtl     <- function(trait.cfg, trait.names, traits, trait.path, funArgs) {
        model <- as.character(trait.cfg$model)
        year  <- as.numeric(trait.cfg$year)
        for( trait in traits ) {
                trait_subsubfolder_fpath <- file.path(trait.path, trait)
                scan.sw.mdres  <- readRDS(file=paste0(trait_subsubfolder_fpath,'/scansw.mdres.rds'))
                qtl.num  <- length(scan.sw.mdres$model$chr)
                if( qtl.num > 0 ) {
                        for (i in 1:qtl.num) {
                                mypos     <- scan.sw.mdres$model$pos[i]
                                mychr     <- as.numeric(scan.sw.mdres$model$chr[i])
                                f         <- which(supermap.bin.df$LG==mychr)
                                x         <- supermap.bin.df[f,]
                                idx       <- which.min(abs(x$consensus-mypos))
                                mymarker  <- paste0(x$binID[idx],'_',x$marker[idx])
                                f2        <- intersect(scan.sw.mdres$model$name[i],rownames(scan.sw.mdres$drop))
                                if (qtl.num==1) {
                                        x.var <- scan.sw.mdres$variance
                                } else {
                                        x2    <- scan.sw.mdres$drop[f2,]
                                        x.var <- x2$'%var'
                                }
                                append.pointer(funArgs, c(model, year, trait.names, trait, mychr, mypos, mymarker, x.var, scan.sw.mdres$variance, scan.sw.mdres$inter[i]))
                        }
                }
        }
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, collate_qtl, qtl.collated.df.p)

qtl.collated.df <- qtl.collated.df.p$value
write.csv(qtl.collated.df, file=paste0(workflow,'/traits/qtl_collated.csv'), row.names=F)

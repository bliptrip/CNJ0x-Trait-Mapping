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

pens.df     <- data.frame(model=character(), mtraits=character(), trait=character(), pen_main=integer(), pen_heavy=integer(), pen_light=integer(), stringsAsFactors=F)
pens.df.p   <- newPointer(pens.df)

#Loop over all mmers, trait groups, subgroups, and perform stepwiseqtl()
constructPensTable <- function(trait.cfg, trait.names, traits, trait.path, funArgs) {
	traits.len = length(traits)
    for( j in 1:traits.len ) {
        trait <- traits[j]
        trait_subsubfolder_fpath = paste0(trait.path, '/', trait) 
        scan.sw <- readRDS(file=paste0(trait_subsubfolder_fpath, "/scansw.rds"))
        pens <- attr(scan.sw, "penalties")
        append.pointer(funArgs,c(trait.cfg$model, trait.cfg$mtraits, trait, pens[["main"]], pens[["heavy"]], pens[["light"]]))
    }
}

#Loop through all legitimate traits and build collated qtl file.
loopThruTraits(workflow, constructPensTable, pens.df.p)
pens.df <- pens.df.p$value
write.csv(pens.df, file=paste0(workflow,'/traits/pens_collated.csv'), row.names=F)

#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
    #Maybe get seed from system time
    process=0
    seed=54955149
    n.perms=2
    query_mmers="2011"
    #query_traits="berry_length,berry_width,berry_weight"
    query_traits="total_berry_weight"
    #query_subtrait="berry_weight"
    query_subtraits="total_berry_weight"
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

query_mmers.v <- unlist(strsplit(query_mmers, split=","))
query_traits.l.v <- unlist(strsplit(query_traits, split=";")) #trait groupings are split by ';', subtraits within by ','
query_subtraits.l.v <- unlist(strsplit(query_subtraits, split=";")) #subtrait groupings are split by ';', query subtraits within by ','

# loading libraries
library(qtl)

cnjpop.mmer.p2 <- readRDS('cnjpop.mmer2.p2.rds')
#Create a mmer list containing all the query_mmers.v plus fields for the process ID on the cluster, seed used, the number of permutations, a timestamp, system information (uname -a), and then the mmers for analysis
#NOTE: Rather than using a description field for mmers, I use the  colname to designate the description of the multivariate analysis.
cnjpop.operms.p2 <- vector("list",length(query_mmers.v)+5)
names(cnjpop.operms.p2) <- c("seed","n.perms","process","timestamp","uname", query_mmers.v)
cnjpop.operms.p2$seed      <- seed
cnjpop.operms.p2$n.perms   <- n.perms
cnjpop.operms.p2$process   <- process
cnjpop.operms.p2$timestamp <- Sys.time()
cnjpop.operms.p2$uname     <- system("uname -a", intern=TRUE)

set.seed(seed)
#Read in the cross files for R/qtl and then calculate the genotype probabilities
for (i in 1:length(cnjpop.mmer.p2)) {
    mmer <- cnjpop.mmer.p2[[i]]
    if (mmer$description %in% query_mmers.v ) {
        cnjpop.operms.p2[[as.character(mmer$description)]] <- vector("list", length(query_traits.l.v))
	k <- 1
	trait.names <- vector()
        for(j in 1:length(mmer$analyses)) {
            analysis <- mmer$analyses[[j]]
            #For the current analysis, see if it is listed in the query_traits.l.v
            analysis.traits <- paste0(analysis$traits,collapse=",")
            query_traits.l.idx <- which( query_traits.l.v %in% analysis.traits )
            if( length(query_traits.l.idx) == 1 ) {
                #Generate a cross datastructure for running analyses on
                qtlfile=paste0('qtl/',mmer$description,"__",paste(analysis$traits,collapse="__"),"__QTL.csv")
                cnjpop.cross <- (read.cross(format = "csv", file=qtlfile, genotypes = NULL))
                cnjpop.cross <- calc.genoprob(cnjpop.cross,step=0,map.function="kosambi") 
                #For the current analysis
                #Derive the columns of interest for running the permutations on.
                #First, splitup the subraits
                query_subtraits.v <- unlist(strsplit(query_subtraits.l.v[query_traits.l.idx], ","))
                #Determine which of the subtraits are in the global traits
                query_subtraits.v.phenos <- query_subtraits.v[which(query_subtraits.v %in% analysis$traits)]
                cnjpop.operms.p2[[as.character(mmer$description)]][[k]] <- scanone(cnjpop.cross, n.perm=n.perms, pheno.col=query_subtraits.v.phenos, verbose=FALSE)
		if( length(query_subtraits.v.phenos) == 1 ) {
		   colnames(cnjpop.operms.p2[[as.character(mmer$description)]][[k]]) <- query_subtraits.v.phenos
		}
		trait.names[k] <- paste0(c(gsub(",","_",analysis$traits),gsub(",","_",query_subtraits.v.phenos)),collapse="___")
		k <- k + 1
	    }
        }
	names(cnjpop.operms.p2[[as.character(mmer$description)]]) <- trait.names
    }
}

rds.filename <- paste0("cnjpop.operms.p2.",as.character(process),".rds")
saveRDS(cnjpop.operms.p2,rds.filename,compress=TRUE)

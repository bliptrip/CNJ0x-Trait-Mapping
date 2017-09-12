#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
    total_perms=1000
    num_clusters_per_trait=1
    #Maybe get seed from system time
    query_mmers="2011"
    #query_traits="berry_length,berry_width,berry_weight"
    query_traits="total_berry_weight"
    #query_subtrait="berry_weight"
    query_subtraits="total_berry_weight"
    folder='runperms'
    output='cnjpop.operms.p2.rds'
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

query_mmers.v <- unlist(strsplit(query_mmers, split=";"))
query_traits.l.v <- unlist(strsplit(query_traits, split=";")) #trait groupings are split by ';', subtraits within by ','
query_subtraits.l.v <- unlist(strsplit(query_subtraits, split=";")) #subtrait groupings are split by ';', query subtraits within by ','

# loading libraries
library(qtl)

#Create a mmer list containing all the query_mmers.v plus fields for the process ID on the cluster, seed used, the number of permutations, a timestamp, system information (uname -a), and then the mmers for analysis
#NOTE: Rather than using a description field for mmers, I use the  colname to designate the description of the multivariate analysis.
cnjpop.operms.p2 <- vector("list",length(query_mmers.v))
names(cnjpop.operms.p2) <- query_mmers.v

#Read in the cross files for R/qtl and then calculate the genotype probabilities
for (mmer.name in names(cnjpop.operms.p2)) {
    cnjpop.operms.p2[[mmer.name]] <- vector("list",length(query_traits.l.v)) 
    for ( i in 1:length(query_traits.l.v) ) {
        trait.names <- paste0(unlist(strsplit(query_traits.l.v[i], ",")),collapse="__")
        query_subtraits.v <- unlist(strsplit(query_subtraits.l.v[i], ","))
        cnjpop.operms.p2[[mmer.name]][[trait.names]] <- vector("list",length(query_subtraits.v))
        names(cnjpop.operms.p2[[mmer.name]][[trait.names]]) <- query_subtraits.v
        for( query_subtrait in query_subtraits.v ) {
            masterfold <- paste0(c(mmer.name,trait.names,query_subtrait),collapse="--")
            masterfold <- paste0(folder,"/",masterfold)
            glob <- Sys.glob(paste0(masterfold,"/[0-9]*"))
            cnjpop.operms.p2[[mmer.name]][[trait.names]][[query_subtrait]] <- vector(mode="list",length=length(glob))
            for( j in 1:length(glob) ) {
                operms.file = paste0(masterfold,"/",as.character(j-1),"/","cnjpop.operms.p2.",as.character(j-1),".rds")
                operms <- readRDS(operms.file)
                #This could be a bit cleaner -- just access element 1 as there is only one subtrait per trait per folder in each analysis (not combining at this point)
                cnjpop.operms.p2[[mmer.name]][[trait.names]][[query_subtrait]][[j]] <- operms[[mmer.name]][[1]]
            }
        }
    }
}

saveRDS(cnjpop.operms.p2,output,compress=TRUE)

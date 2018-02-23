#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.

library(dplyr)
library(factoextra)

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

#source(paste0(workflow,"/configs/model.cfg"))

qtl.collated.df <- read.csv(file=paste0(workflow,'/traits/qtl_collated.csv'), head=TRUE)
#Only look at those qtls derived from the specified method, and moreover, remove interaction qtls and focus on additive-only ones.
qtl.collated.df <- qtl.collated.df %>%
                        filter(is.na(chr2) & is.na(position2)) %>%
                        arrange(method,mtraits,trait,chr,model)

#Group by trait and chromosome, and then apply hierarchical agglomerative clustering to see if can create consensus for some of the markers
qtl.collated.grouped.df <- qtl.collated.df %>% group_by(method,mtraits,trait,chr)

generate_consensus <- function(model, mtraits, trait, chr, position) {
    #Initialize the consensus position to the current positions
    consensus_positions <- position
    if( length(position) > 1 ) {
        group.dist        <- dist(position)
        group.hclust      <- hclust(group.dist)
        clust             <- cutree(group.hclust, h = 5) #Cut at height 5
        for( i in 1:length(clust) ) {
            f <- which(clust == i)
            position.mean <- mean(position[f])
            consensus_positions[f] <- position.mean
        }
    }
    return(consensus_positions)
}

qtl.collated.consensus.df <- qtl.collated.grouped.df %>% mutate(position.consensus = generate_consensus(model, mtraits, trait, chr, position))
write.csv(qtl.collated.consensus.df, file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), row.names=FALSE)

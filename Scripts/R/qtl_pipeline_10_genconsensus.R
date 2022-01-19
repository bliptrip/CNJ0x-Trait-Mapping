library(dplyr)
library(factoextra)

# loading libraries
source('./usefulFunctions.R')

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

#source(paste0(workflow,"/configs/model.cfg"))

qtl.collated.orig.df <- read.csv(file=paste0(workflow,'/traits/qtl_collated.csv'), head=TRUE)
#Only look at those qtls derived from the specified method, and moreover, remove interaction qtls and focus on additive-only ones.
qtl.collated.df <- qtl.collated.orig.df %>%
                        filter(is.na(chr2) & is.na(position2)) %>%
                        arrange(method,trait,chr,model)

#Group by trait and chromosome, and then apply hierarchical agglomerative clustering to see if can create consensus for some of the markers
qtl.collated.grouped.df <- qtl.collated.df %>% group_by(method,trait,chr)

generate_consensus <- function(model, trait, chr, position) {
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

qtl.collated.consensus.df <- qtl.collated.grouped.df %>% ungroup() %>% mutate(position.consensus = generate_consensus(model, trait, chr, position))
qtl.collated.epistatics.df <- qtl.collated.orig.df %>% 
                                    filter(!is.na(chr2) & !is.na(position2)) %>%
                                    mutate(position.consensus = NA)

qtl.collated.both.df <- rbind(qtl.collated.consensus.df, qtl.collated.epistatics.df) %>%
                            arrange(method,trait,chr,model)
write.csv(qtl.collated.both.df, file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), row.names=FALSE)

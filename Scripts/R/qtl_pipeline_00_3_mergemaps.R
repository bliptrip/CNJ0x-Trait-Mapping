#!/usr/bin/env RScript

#This script is intended to merge uniparental maps into an integrated map using LPMerge utility.
#NOTE: These maps were provided by Luis Diaz-Garcia when constructing a new Cranberry genome
# from PacBio Data.  He used the GBS data taken from https://academic.oup.com/g3journal/article/7/4/1177/6031795
# to call new SNPs and generate linkage maps for assisting in genome construction (ihttps://www.frontiersin.org/articles/10.3389/fpls.2021.633310/full)

source('./usefulFunctions.R')
library(LPmerge)
library(tidyverse)

max_interval <- get0('max_interval', ifnotfound=10) #Maximum interval to test up to in LPmerge function call

args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

map_names <- c("gmap_cnj02-p1.csv","gmap_cnj02-p2.csv")
maps <- list()
for (i in 1:length(map_names)) {
    filename  <- geno_rpath2fpath(map_names[i])
    maps[[i]] <- as.data.frame(read_csv(filename))
}

lpmerge_results <- list()
for( j in 1:12 ) {
    lg_maps <- list()
    for (i in 1:length(map_names)) {
        lg_maps[[i]] <- as.data.frame(maps[[i]] %>% 
            filter(chr == j) %>% 
            select(marker,pos))
    }
    lpmerge_results[[j]] <- LPmerge(lg_maps, max.interval = 1:max_interval)
}

#Find the best of each 
lpmerge_best <- lapply(lpmerge_results, function(x) { x[[which.min(unlist(lapply(x, function(y) { attr(y,"RMSE.mean") + attr(y,"RMSE.sd") })))]] })
for( i in 1:12 ) {
    lpmerge_best[[i]] <- lpmerge_best[[i]] %>% mutate(LG=i) #Re-inject the chromosome number before merging dataset
}

lpmerge_map <- reduce(lpmerge_best, function(a, n) { a <- a %>% bind_rows(n) } ) %>% select(marker, consensus, LG)
write_csv(lpmerge_map,geno_rpath2fpath("gmap.cnj02.consensus.csv"))

#Now merge genotype calls datasets
geno_names <- c("geno_cnj02-p1.csv","geno_cnj02-p2.csv")
genos <- list()
for (i in 1:length(geno_names)) {
    filename  <- geno_rpath2fpath(geno_names[i])
    genos[[i]] <- as.data.frame(read_csv(filename))
    genos[[i]] <- genos[[i]][order(genos[[i]]$id),]
}


markers_common <- maps[[1]]$marker
markers_all <- maps[[1]]$marker
for (i in 2:length(map_names)) {
    markers_common <- base::intersect(markers_common,maps[[i]]$marker)
    markers_all <- base::union(markers_all,maps[[i]]$marker)
}
markers_not_common <- markers_all[!(markers_all %in% markers_common)]
mismatch_idx <- genos[[1]][,c("id",markers_common)]!=genos[[2]][,c("id",markers_common)]
genos_common_1_mismatch <- genos[[1]][,c("id",markers_common)]
row.names(genos_common_1_mismatch) <- genos_common_1_mismatch[,"id"]
genos_common_1_mismatch[!mismatch_idx] <- NA
genos_common_1_mismatch <- genos_common_1_mismatch[,markers_common]
genos_common_2_mismatch <- genos[[2]][,c("id",markers_common)]
row.names(genos_common_2_mismatch) <- genos_common_2_mismatch[,"id"]
genos_common_2_mismatch[!mismatch_idx] <- NA
genos_common_2_mismatch <- genos_common_2_mismatch[,markers_common]
genos_common_mismatch = array(NA, dim=dim(genos_common_1_mismatch))
for( i in 1:nrow(genos_common_1_mismatch) ) {
    for( j in 1:ncol(genos_common_1_mismatch) ) {
        if( !is.na(genos_common_1_mismatch[i,j]) ) {
            genos_common_mismatch[i,j] <- paste0(genos_common_1_mismatch[i,j],",",genos_common_2_mismatch[i,j])
        }
    }
}

genos_common_mismatch <- merge(genos_common_1_mismatch,genos_common_2_mismatch)
genos_common_mismatch <- paste0(genos_common_1_mismatch, genos_common_2_mismatch, sep=",")
genos_common_fixed <- genos[[1]][,c("id",markers_common)]
genos_common_fixed[mismatch_idx] <- "AB" #Fix mismatches
genos_not_common <- lapply(genos, function(x) { x %>% select(!markers_common) })
genos.merged <- reduce(genos_not_common, function(a,n) { a <- a %>% inner_join(n, by="id") })
genos.all.merged <- genos.merged %>% inner_join(genos_common_fixed, by="id")
write_csv(genos.all.merged,geno_rpath2fpath("geno_cnj02.consensus.csv"))


#
#
##Extract out shared markers, which will have the AB x AB segregation type
#
#maps_binned <- list()
#for (i in 1:length(map_names)) {
#    map <- maps[[i]]
#    maps_binned[[i]] <- map %>%
#                            mutate(pos = round(pos, digits=3)) %>% #LPMerge rounds to 3 digits, and we need to bin at this level and then remove redundancies
#                            group_by(chr, pos) %>%
#                            summarize(bin_marker = ifelse(n() > 1, 
#                                                      ifelse(any(marker %in% markers_common), 
#                                                             marker[marker %in% markers_common][1], 
#                                                             marker[1]), 
#                                                      marker[1]),
#                                      .groups="drop")
#}


save.image(paste0(workflow,"/.RData.1_mergemaps"))

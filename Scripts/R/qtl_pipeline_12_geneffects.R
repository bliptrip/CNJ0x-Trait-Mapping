#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
#install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

library(jsonlite)
library(qtl)
library(tidyverse)

#Defaults
workflow="../../Workflows/1"
qtl_method   <- "stepwiseqtl" #Which method to filter
num_top_qtls <- 2 #Number of the top QTLs to calculate effect sizes for

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source('./usefulFunctions.R')
source(paste0(workflow,"/configs/model.cfg"))

extract_effects <- function(method, model, trait, chr, position) {
    print(paste0("method = ",method, ", model = ", model, "trait = ", trait))
    qtl.mname <- paste0(chr,'@',format(round(as.numeric(position),digits=1), nsmall=1))
    cross <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/cross.rds'))
    qtl   <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/',trait,'/', ifelse(method == 'scanone', 'scanone.qtl', 'scansw'), '.rds'))
    #The qtl$prob contains the list of significant QTLs and the probability of a given genotype at the QTL.  I would like to show a boxplot of
    #blup values at the different genotypes for each QTL, but since the genotype is a mixed distribution at each QTL, I will only include genotypes with a higher
    #than, say, 95% probability of being a given genotype, and choose that as the representative genotype (assigning the entire BLUP and/or trait to that genotype for organization)
    #Determine which index corresponds to QTL
    qtl.i <- which(qtl$name == qtl.mname)
    print(paste0("qtl.mname = ", qtl.mname, ", qtl.i = ",qtl.i))
    print(paste0("qtl$name = ", qtl$name))
    qtl.p <- data.frame(qtl$prob[[qtl.i]])
    qtl.p$id.i <- rownames(qtl.p)
    qtl.p.nest <- qtl.p %>% 
                    mutate(blup=cross$pheno[,trait]) %>% 
                    pivot_longer(!c(id.i,blup),names_to='genotype', values_to='probability') %>% 
                    filter(probability > 0.95) %>%
                    mutate(genotype=paste0(genotype,'.blups')) %>%
                    group_by(genotype) %>%
                    nest() %>%
                    pivot_wider(names_from="genotype", values_from="data")
    m.effects <- effectplot(cross, pheno.col=trait, mname1=qtl.mname, draw=FALSE)
    return(list(m.effects, qtl.p.nest))
}

#For each qtl in the collated file, use it's position and consensus position to calculate the effects.  Store this information in the collated file?
generate_collated_effects <- function(qtl.collated.tb,position_col,num_top_qtls) {
    qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                    arrange(method,trait,model,desc(marker.variance)) %>%
                                    group_by(method,trait,model) %>%
                                    mutate(top_qtls = c(rep(TRUE,num_top_qtls),rep(FALSE,length(marker.variance)-num_top_qtls))) %>%
                                    filter(top_qtls == TRUE)
    effects.tb <- eval(parse(text=paste0("qtl.collated.filtered.tb %>% select(method, model, trait, chr, ", position_col, ")")))
    effects.tb <- effects.tb %>%
                    mutate(effectMarker="",
                           AC=0,
                           BC=0,
                           AD=0,
                           BD=0,
                           ACSE=0,
                           BCSE=0,
                           ADSE=0,
                           BDSE=0,
                           AC.blups=list(1),
                           AD.blups=list(1),
                           BC.blups=list(1),
                           BD.blups=list(1))
    for( i in 1:nrow(effects.tb) ) {
        e.tb <- effects.tb[i,]
        effs.l <- extract_effects(e.tb$method, e.tb$model, e.tb$trait, e.tb$chr, e.tb[position_col])
        effs <- effs.l[[1]]
        blups <- effs.l[[2]]
        re <- "(.+)\\.((AC)|(BC)|(AD)|(BD))$"
        marker <- unique(gsub(re, "\\1", names(effs$Means), fixed=FALSE))
        cm <- gsub(re, "\\2", names(effs$Means), fixed=FALSE)
        names(effs$Means) <- cm
        cs  <- gsub(re, "\\2SE", names(effs$SEs), fixed=FALSE)
        names(effs$SEs) <- cs
        e.tb$effectMarker <- marker
        e.tb[c(cm, cs)] <- t(as.numeric(cbind(effs$Means, effs$SE)))
        e.tb$AC.blups <- blups$AC.blups
        e.tb$BC.blups <- blups$BC.blups
        e.tb$AD.blups <- blups$AD.blups
        e.tb$BD.blups <- blups$BD.blups
        effects.tb[i,] = e.tb
    }
    effects.collated.tb  <- qtl.collated.filtered.tb %>%
                            left_join(effects.tb, by=c("method","trait","model","chr",position_col))
    return(effects.collated.tb)
}

generate_effects_plot <- function(e.tb, methods=NULL, models=NULL, traits=NULL, num_top_qtls) {
    if( methods != NULL ) {
        e.tb <- filter(e.tb, method %in% methods)
    }
    if( models != NULL ) {
        e.tb <- filter(e.tb, model %in% methods)
    }
    if( traits != NULL ) {
        e.tb <- filter(e.tb, trait %in% traits)
    }
    e.tb <- e.tb %>% mutate(top_qtls = c(rep(TRUE,num_top_qtls),rep(FALSE,length(marker.variance)-num_top_qtls))) %>%
                     filter(top_qtls == TRUE)
}
qtl.tb  <- read_csv(file=paste0(workflow,'/traits/qtl_collated.csv'), col_names=TRUE)
effs.tb <- generate_collated_effects(qtl.tb, "position", num_top_qtls)
saveRDS(effs.tb,file=paste0(workflow,'/traits/effects_collated.rds'), compress=TRUE)
#write.csv(effs.tb, file=paste0(workflow,'/traits/effects_collated.csv'), row.names=FALSE)

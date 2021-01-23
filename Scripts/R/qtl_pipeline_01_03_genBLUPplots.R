#!/usr/bin/env RScript

#NOTE: This particular script generates informative data plots for the BLUPs: Things like pairwise correlation plots, boxplots, etc.
# Think of combining BLUP data into a table format similar to the input phenotype table, and then can display data in a way similar
# to how this was done for the original phenotypes.
#
# GGally is our friend here!

# loading libraries
library(GGally)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(tidyverse)

source('./usefulFunctions.R')

workflow <- "../../Workflows/1"

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

collateBLUPs <- function(trait.cfg, trait.path, models.tb.p) {
    trait         <- trait.cfg$trait
    model         <- readRDS(file=paste0(trait.path,"/mmer.rds"))
    u             <- model$U[["u:id"]]
    blups         <- u[[trait]]
    blupsZ        <- scale(blups)[,1]
    if( trait.cfg$model == "all-years" ) {
        pheno <- model$data[,c("id",trait)] %>% group_by(id) %>% summarize(mean=mean(.data[[trait]]))
        pheno <- pheno[["mean"]]
    } else {
        pheno         <- model$data[[trait]]
    }
    phenoZ        <- scale(pheno)[,1]
    ids           <- names(blups)
    #print(paste0("Model: ",trait.cfg$model, ", Trait: ", trait, ", Num blups: ", length(blups), " Num phenos: ", length(pheno)))
    model_label   <- ifelse(trait.cfg$model == 'all-years', 'All Years', trait.cfg$model)
    models.tb.p$value  <- models.tb.p$value %>% 
                            add_row(id=factor(ids), model=factor(trait.cfg$model, labels=model_label), trait=factor(trait, labels=factor(trait.cfg$label)), value=blups, valueZ=blupsZ, type=factor("BLUPs"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short)) %>% 
                            add_row(id=factor(ids), model=factor(trait.cfg$model, labels=model_label), trait=factor(trait, labels=factor(trait.cfg$label)), value=pheno, valueZ=phenoZ, type=factor("Raw"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short)) 
}

traits.df   <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=T)
#Just build a 'long' version of collated BLUPs as tibble table, and then use tidyverse 'pivot_wider()' to flatten
model.collated.long.tb     <- tibble(id=factor(), model=factor(), trait=factor(), value=numeric(), valueZ=numeric(), type=factor(), label=factor(), label_short=factor())
model.collated.long.tb.p   <- newPointer(model.collated.long.tb)
loopThruTraits(workflow, collateBLUPs, loopArgs=model.collated.long.tb.p)
model.collated.long.tb     <- model.collated.long.tb.p$value
#Save long version, as it can be manipulated in many ways to produce meaningful plots/graphs
saveRDS(model.collated.long.tb, file=paste0(workflow,'/traits/blups_collated.long.rds'), compress=TRUE)

#Spread out
model.collated.tb <- spread(model.collated.long.tb, trait, blup)
write.csv(model.collated.tb, file=paste0(workflow,'/traits/blups_collated.csv'), row.names=F)

#Generate one plot per unique model

#Generate plot across only models with valid 'year' -- not 'all-years'
model.collated.years.tb <- model.collated.tb %>% filter(model != 'all-years')
full_traits <- c()
for (trait in unique(model.collated.long.tb$trait)) {
    if (!any(is.na(model.collated.years.tb[,trait]))) {
        full_traits <- c(full_traits, trait)
    }
}
g <- ggscatmat(model.collated.years.tb[,c('id','model',full_traits)], color="model", alpha=0.8) +
         theme_minimal() +
         theme(axis.title = element_blank(),
               axis.text.x  = element_text(size=8,angle=60),
               axis.text.y  = element_text(size=8),
               strip.text = element_text(face="bold",size=8),
               legend.text=element_text(face="bold",size=8),
               legend.title=element_text(face="bold",size=10))
png(filename=paste0(workflow,'/traits/plots/blups_collated.traits.years.png'), width=1280, height=960, bg="white")
g
dev.off()

#Generate a plot across all models.  NOTE: Need to remove any traits with any NA's as these will result in ggscatmat returning an error (unbalanced)
model.collated.allmodels.tb <- model.collated.tb
full_traits <- c()
for (trait in unique(model.collated.long.tb$trait)) {
    if (!any(is.na(model.collated.allmodels.tb[,trait]))) {
        full_traits <- c(full_traits, trait)
    }
}
g <- ggscatmat(model.collated.allmodels.tb[,c('id','model',full_traits)], color="model", alpha=0.8) +
         theme_minimal() +
         theme(axis.title = element_blank(),
               axis.text.x  = element_text(size=8,angle=60),
               axis.text.y  = element_text(size=8),
               strip.text = element_text(face="bold",size=10),
               legend.text=element_text(face="bold",size=8),
               legend.title=element_text(face="bold",size=10))
png(filename=paste0(workflow,'/traits/plots/blups_collated.traits.allmodels.png'), width=1280, height=960, bg="white")
g
dev.off()

model.collated.long.tb <- model.collated.long.tb %>% 
                            mutate(model_trait=as.factor(paste0(model,' ',trait)), model_type=as.factor(paste0(model,' ',type)))

model.collated.long.stats.tb <- model.collated.long.tb %>%
                                    group_by(model,trait,type) %>%
                                    summarize(min=min(value), minZ=min(valueZ), min_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.min(value)]),
                                              max=max(value), maxZ=max(valueZ), max_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.max(value)]))

model.collated.long.ylims.tb <- model.collated.long.tb %>%
                                    group_by(type) %>%
                                    summarize(min=min(value, na.rm=TRUE) - 2.5, minZ=min(valueZ, na.rm=TRUE)-1,
                                              max=max(value, na.rm=TRUE) + 2.5, maxZ=max(valueZ, na.rm=TRUE)+1) %>%
                                    gather("limit","value",-"type")
model.collated.long.ylims.tb$model <- model.collated.long.tb$model[1]

#Display boxplots of Raw Values + BLUPs without Z-score scaling
g <- model.collated.long.tb %>% 
        ggplot(aes(x=factor(model), y=value, fill=factor(trait))) +
        geom_blank(data=model.collated.long.ylims.tb %>% filter((limit == "min") | (limit == "max")), aes(x=factor(model), y=value)) +
        geom_boxplot(alpha=0.4) +
        geom_jitter(width=0.1, alpha=0.4) +
        geom_text(data=model.collated.long.stats.tb, aes(y=min, label=min_geno), size=8, fontface='bold', vjust=1.5, hjust=1.5, angle=45) +
        geom_text(data=model.collated.long.stats.tb, aes(y=max, label=max_geno), size=8, fontface='bold', vjust=-1.5, hjust=0, angle=45) +
        facet_grid(rows=vars(type), cols=vars(trait), scales="free_y") +
        guides(fill = 'none') +
        xlab("Year") + 
        ylab("Trait Values") +
        ggtitle(label="CNJ02 Population Raw Values and BLUPs", subtitle="Yield-Related Traits") +
        theme( axis.text.x = element_text(face="bold", size=32, angle = 60, hjust = 1),
               axis.text.y = element_text(face="bold", size=32),
               strip.text  = element_text(face="bold", size=36, angle = 30),
               axis.title  = element_text(face="bold",size=48),
               plot.title   = element_text(face="bold",size=52, hjust=0.5),
               plot.subtitle = element_text(size=48, hjust = 0.5))
png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.allmodels.png'), width=2560, height=1920, bg="white")
g
dev.off()

#Display boxplots of Raw Values + BLUPs with Z-score scaling
g <- model.collated.long.tb %>% 
        ggplot(aes(x=factor(model), y=valueZ, fill=factor(trait))) +
        geom_blank(data=model.collated.long.ylims.tb %>% filter((limit == "minZ") | (limit == "maxZ")), aes(x=factor(model), y=value)) +
        geom_boxplot(alpha=0.4) +
        geom_jitter(width=0.1, alpha=0.4) +
        geom_text(data=model.collated.long.stats.tb, aes(y=minZ, label=min_geno), size=8, fontface='bold', vjust=1.5, hjust=1.5, angle=45) +
        geom_text(data=model.collated.long.stats.tb, aes(y=maxZ, label=max_geno), size=8, fontface='bold', vjust=-1.5, hjust=0, angle=45) +
        facet_grid(rows=vars(type), cols=vars(trait), scales="free_y") +
        guides(fill = 'none') +
        xlab("Year") + 
        ylab("Trait Values") +
        ggtitle(label="CNJ02 Population Raw Values and BLUPs Z-Score Scaled", subtitle="Yield-Related Traits") +
        theme( axis.text.x = element_text(face="bold", size=32, angle = 60, hjust = 1),
               axis.text.y = element_text(face="bold", size=32),
               strip.text  = element_text(face="bold", size=36, angle = 30),
               axis.title  = element_text(face="bold",size=48),
               plot.title   = element_text(face="bold",size=52, hjust=0.5),
               plot.subtitle = element_text(size=48, hjust = 0.5))
png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.allmodels.zscore.png'), width=2560, height=1920, bg="white")
g
dev.off()

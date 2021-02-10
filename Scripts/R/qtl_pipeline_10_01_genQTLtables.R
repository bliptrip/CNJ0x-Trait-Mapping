#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(formattable)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(tidyverse)

source('./usefulFunctions.R')

#Defaults (can be overridden with command-line invokation
workflow     <- "../../Workflows/1"
qtl_method   <- "stepwiseqtl" #Which method to filter
num_top_qtls <- 2 #Number of top QTLs to show per trait


#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

trait.cfg.tb    <- read_csv(file=paste0(workflow,'/configs/model-traits.cfg.csv'), col_names=TRUE)
qtl.collated.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), col_names=TRUE)
consensusMapWithGenes.tb <- read_csv(file=geno_dpath2fpath("consensusMapAll2_withGenes.csv"), col_names=TRUE)

consensusMapWithGenes.filtered.tb <- consensusMapWithGenes.tb %>% 
                                        select(marker, gaccession, blast1, blast2) %>%
                                        mutate(blast1 = gsub("0", "",blast1)) %>%
                                        mutate(blast2 = gsub("0", "",blast2))

qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                filter(method == qtl_method) %>%
                                mutate(marker = gsub(".+cM_(.+)", "\\1",nearest.marker), model_name = model_to_name(trait.cfg.tb,model,trait), trait_name = trait_to_name(trait.cfg.tb,model,trait)) %>%
                                left_join(consensusMapWithGenes.filtered.tb, by="marker") %>%
                                arrange(trait,model,desc(marker.variance)) %>%
                                mutate(marker.variance=percent(marker.variance/100), 
                                       trait_repeat=(trait_name == c("",trait_name[-length(trait_name)])),
                                       model_repeat=(model_name == c("",model_name[-length(model_name)]))) %>%
                                group_by(trait_name, model_name) %>%
                                mutate(top_qtls = c(rep(TRUE,num_top_qtls),rep(FALSE,length(marker.variance)-num_top_qtls))) %>%
                                filter(top_qtls == TRUE)


qtl.collated.succinct.tb <- qtl.collated.filtered.tb %>%
                                select(trait_name,model_name,chr,position,marker.variance,marker,blast1,blast2,trait_repeat,model_repeat) %>%
                                mutate(marker = ifelse(is.na(marker)," ",marker)) %>%
                                mutate(blast1 = ifelse(is.na(blast1)," ",blast1)) %>%
                                mutate(blast2 = ifelse(is.na(blast2)," ",blast2)) %>%
                                mutate(position = color_bar("lightgreen")(round(position,digits=1)),
                                       marker.variance = color_bar("lightblue")(marker.variance),
                                       trait_name = ifelse(trait_repeat, "", cell_spec(trait_name, "html", bold=TRUE)),
                                       model_name = ifelse(model_repeat, "", cell_spec(model_name, "html", bold=TRUE)))

colnames(qtl.collated.succinct.tb) <- c("Trait", "Model", "Linkage Group", "Marker Location (cM)", "Variance Explained by QTL", "Nearest Marker", "Putative Function 1", "Putative Function 2", "trait_repeat", "model_repeat")

qtl.f <- formattable(qtl.collated.succinct.tb, formatters=list(
                        `Variance Explained by QTL` = color_bar("lightblue"),
                        `Trait` = formatter("span", style = ~ style(visibility=ifelse(trait_repeat,"hidden","visible"), font.weight = ifelse(trait_repeat,NA,"bold"))),
                        `Model` = formatter("span", style = ~ style(visibility=ifelse(model_repeat,"hidden","visible"), font.weight = ifelse(model_repeat,NA,"bold"))),
                        model_repeat = FALSE, #Used to hide this column
                        trait_repeat = FALSE),
                      format="html",
                      digits=3,
                      align=c('l', 'l', 'r', 'r', 'r', 'l', 'l', 'l')
                    )


qtl.collated.succinct.tb %>%
    select(-trait_repeat,-model_repeat) %>%
    kable("html", align='llrlrlll', escape = FALSE) %>%
    kable_paper("striped", full_width=FALSE) %>%
    column_spec(3, width = "1cm") %>%
    column_spec(4, width = "2.5cm") %>%
    column_spec(5, width = "3cm") %>%
    column_spec(7, width = "10cm") %>%
    column_spec(8, width = "8cm")

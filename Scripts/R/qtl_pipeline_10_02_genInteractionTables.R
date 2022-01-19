#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(flextable)
library(formattable)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(tidyverse)

source('./usefulFunctions.R')

#Defaults (can be overridden with command-line invocation)
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top interaction QTLs to show per trait
table.font_size <- get0("table.font_size", ifnotfound=10)
table.width     <- get0("table.width", ifnotfound=0.85)

qtl_scan_method <- "stepwiseqtl" #Only stepwiseqtl method valid for interactions

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
num_top_qtls <- as.numeric(num_top_qtls) #Convert to numeric if passed in through command-line

trait.cfg.tb    <- read_csv(file=paste0(workflow,'/configs/model-traits.cfg.csv'), col_names=TRUE)
qtl.collated.ints.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.csv'), col_names=TRUE) %>%
                            filter(!is.na(chr2) & !is.na(position2)) %>% #Only include interaction QTLs for now
                            group_by(trait) %>%
                            mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>% #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                            mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                            ungroup()

qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                filter(method == qtl_scan_method) %>%
                                mutate(marker = gsub(".+cM_(.+)", "\\1",nearest.marker), 
                                       model_name = model_to_name(trait.cfg.tb,model,trait), 
                                       trait_name = trait_to_name(trait.cfg.tb,model,trait)) %>%
                                arrange(trait_name,model_name,desc(marker.variance)) %>%
                                group_by(trait_name,model_name) %>%
                                mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker.variance))),rep(FALSE,length(marker.variance)-min(num_top_qtls,length(marker.variance))))) %>%
                                filter(top_qtls == TRUE) %>%
                                ungroup()

qtl.collated.succinct.tb <- qtl.collated.filtered.tb %>%
                                mutate(marker = ifelse(is.na(marker)," ",marker)) %>%
                                mutate(position = round.digits(position,2),
                                       position2 = round.digits(position2,2)
                                       qtl.lod = paste0(round.digits(qtl.lod,2),"$^{",unlist(map(qtl.pvalue,siginfo)),"}$"))


generateTableSignifSymbols <- function(tbl) {
    return( tbl %>%
                mutate(trait_name=paste0(trait_name,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$"),
                       model_name=paste0(model_name,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) )
}

generateTableRemoveRepeats <- function(tbl) {
    return( tbl %>%
        mutate(trait_repeat=(trait_name == c("",trait_name[-length(trait_name)])),
                model_repeat=(model_name == c("",model_name[-length(model_name)]))) %>%
        mutate(trait_name = ifelse(trait_repeat, "",trait_name),
               model_name = ifelse(model_repeat, "",model_name)) )
}

generateReducedTable <- function(tbl, caption=NULL) {
    qtable1 <- tbl %>%
        generateTableRemoveRepeats() %>%
        select(trait_name,chr,position,chr2,position2,qtl.lod,marker.variance) %>%
        mutate(marker.variance=color_bar("lightblue")(percent(marker.variance/100))) %>%
        rename("Trait"=trait_name,
               "LG1"=chr,
               "Position1 (cM)"=position,
               "LG2"=chr,
               "Position2 (cM)"=position,
               "pLOD"=qtl.lod,
               "Variance Explained by QTL"=marker.variance)
    if( is_html_output() ) {
        qtable2 <- qtable1 %>% kable("html", caption=caption, align='lrrrrrc', escape = FALSE, table.attr="id=\"kableTable\"")
    } else {
        qtable2 <- qtable1 %>% kable(align='lrrrrrc', caption=caption, escape = FALSE)
    }
    qtable3  <- qtable2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        column_spec(1, bold=TRUE) %>%
        column_spec(2, width = "0.5cm") %>%
        column_spec(3, width = "1.5cm") %>%
        column_spec(4, width = "0.5cm") %>%
        column_spec(5, width = "1.5cm") %>%
        column_spec(6, width = "1.5cm") %>%
        column_spec(7, width = "2cm")
    return(qtable3)
}

generateTable <- function(tbl, caption=NULL) {
    qtable1 <- tbl %>%
        generateTableSignifSymbols() %>%
        generateTableRemoveRepeats() %>%
        select(trait_name,model_name,chr,position,chr2,position2,qtl.lod,marker.variance) %>%
        mutate(marker.variance=color_bar("lightblue")(percent(marker.variance/100))) %>%
        rename('Trait[note]'=trait_name,
               'Model[note]'=model_name,
               "LG1"=chr,
               "Position1 (cM)"=position,
               "LG2"=chr,
               "Position2 (cM)"=position,
               "pLOD"=qtl.lod,
               "Variance Explained by QTL"=marker.variance)
    if( is_html_output() ) {
        qtable2 <- qtable1 %>% kable("html", caption=caption, align='llrrrrrc', escape = FALSE, table.attr="id=\"kableTable\"")
    } else {
        qtable2 <- qtable1 %>% kable(align='llrrrrrc', caption=caption, escape = FALSE)
    }
    qtable3  <- qtable2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        column_spec(1, bold=TRUE) %>%
        column_spec(2, bold=TRUE) %>%
        column_spec(3, width = "0.5cm") %>%
        column_spec(4, width = "1.5cm") %>%
        column_spec(5, width = "0.5cm") %>%
        column_spec(6, width = "1.5cm") %>%
        column_spec(7, width = "1.5cm") %>%
        column_spec(8, width = "2cm") %>%
        add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                    "Significance codes for Genotype Effects"))
    return(qtable3)
}

save.image(paste0(workflow,"/.RData.10_02.genInteractionTable.",num_top_qtls))
load(paste0(workflow,"/.RData.10_02.genInteractionTable.",num_top_qtls))

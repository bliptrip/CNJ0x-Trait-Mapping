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
qtl_scan_method      <- get0("qtl_scan_method", ifnotfound="stepwiseqtl") #Which method to filter
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top QTLs to show per trait
table.font_size <- get0("table.font_size", ifnotfound=10)
table.width <- get0("table.width", ifnotfound=0.85)

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
num_top_qtls <- as.numeric(num_top_qtls) #Convert to numeric if passed in through command-line

trait.cfg.tb    <- read_csv(file=paste0(workflow,'/configs/model-traits.cfg.csv'), col_names=TRUE)
qtl.collated.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.csv'), col_names=TRUE) %>%
                        filter(is.na(chr2) & is.na(position2)) %>% #Remove interaction effects for now
                        group_by(trait) %>%
                        mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>% #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                        mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                        ungroup()

consensusMapWithGenes.tb <- read_csv(file=geno_dpath2fpath("consensusMapAll2_withGenes.csv"), col_names=TRUE)
consensusMapWithGenes.filtered.tb <- consensusMapWithGenes.tb %>% 
                                        select(marker, gaccession, blast1, blast2) %>%
                                        mutate(blast1 = gsub("0", "",blast1)) %>%
                                        mutate(blast2 = gsub("0", "",blast2))
qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                filter(method == qtl_scan_method) %>%
                                mutate(marker = gsub(".+cM_(.+)", "\\1",nearest.marker), 
                                       model_name = model_to_name(trait.cfg.tb,model,trait), 
                                       trait_name = trait_to_name(trait.cfg.tb,model,trait)) %>%
                                inner_join(consensusMapWithGenes.filtered.tb, by="marker") %>%
                                arrange(trait_name,model_name,desc(marker.variance)) %>%
                                group_by(trait_name,model_name) %>%
                                mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker.variance))),rep(FALSE,length(marker.variance)-min(num_top_qtls,length(marker.variance))))) %>%
                                filter(top_qtls == TRUE) %>%
                                ungroup()

qtl.collated.succinct.tb <- qtl.collated.filtered.tb %>%
                                mutate(marker = ifelse(is.na(marker)," ",marker)) %>%
                                mutate(blast = ifelse((!is.na(blast1) & !is.na(blast2)),paste0(blast1,";",blast2),ifelse(is.na(blast1),ifelse(is.na(blast2)," ", blast2),blast1))) %>%
                                mutate(position = paste0(round.digits(position,2),"\u00b1",round.digits(interval/2,2)),
                                       qtl.lod = paste0(round.digits(qtl.lod,2),"$^{",unlist(map(qtl.pvalue,siginfo)),"}$"))


generateTableSignifSymbols <- function(tbl) {
    return( tbl %>%
                mutate(trait_name=paste0(trait_name,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$"),
                       model_name=paste0(model_name,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) )
}

#arrange(trait_name,model_name,desc(marker.variance)) %>%
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
        select(trait_name,chr,position,qtl.lod,marker.variance,marker,blast) %>%
        mutate(marker.variance=color_bar("lightblue")(percent(marker.variance/100))) %>%
        rename("Trait"=trait_name,
               "LG"=chr,
               "Marker Location ± 1.5LOD (cM)"=position,
               "pLOD"=qtl.lod,
               "Variance Explained by QTL"=marker.variance,
               "Nearest Marker"=marker,
               "Putative Function"=blast)
    if( is_html_output() ) {
        qtable2 <- qtable1 %>% kable("html", caption=caption, align='lrrrrll', escape = FALSE, table.attr="id=\"kableTable\"")
    } else {
        qtable2 <- qtable1 %>% kable(align='lrrrrll', caption=caption, escape = FALSE)
    }
    qtable3  <- qtable2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        column_spec(1, bold=TRUE) %>%
        column_spec(2, width = "1cm") %>%
        column_spec(3, width = "2cm") %>%
        column_spec(4, width = "2cm") %>%
        column_spec(7, width = "7cm")
    return(qtable3)
}

generateTable <- function(tbl, caption=NULL) {
    qtable1 <- tbl %>%
        generateTableSignifSymbols() %>%
        generateTableRemoveRepeats() %>%
        select(trait_name,model_name,chr,position,qtl.lod,marker.variance,marker,blast) %>%
        mutate(marker.variance=color_bar("lightblue")(percent(marker.variance/100))) %>%
        rename('Trait[note]'=trait_name,
               'Model[note]'=model_name,
               "LG"=chr,
               "Marker Location ± 1.5LOD (cM)"=position,
               "pLOD[note]"=qtl.lod,
               "Variance Explained by QTL"=marker.variance,
               "Nearest Marker"=marker,
               "Putative Function"=blast)
    if( is_html_output() ) {
        qtable2 <- qtable1 %>% kable("html", caption=caption, align='llrrrrll', escape = FALSE, table.attr="id=\"kableTable\"")
    } else {
        qtable2 <- qtable1 %>% kable(align='llrrrrll', caption=caption, escape = FALSE)
    }
    qtable3  <- qtable2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        row_spec(row=which(tbl$trait_name != "")-1, hline_after = TRUE) %>%
        column_spec(1, bold=TRUE) %>%
        column_spec(3, width = "1cm") %>%
        column_spec(4, width = "2.5cm") %>%
        column_spec(5, width = "2cm") %>%
        column_spec(8, width = "20cm") %>%
        add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                    "Significance codes for Genotype Effects",
                    "Significance codes from QTL pvalues"))
    return(qtable3)
}

generateFlexTable <- function(tbl, caption=NULL) {
        flxtable <- qtl.collated.succinct.tb %>%
        select(-trait_repeat,-model_repeat) %>%
        mutate(Trait=gsub("$","",Trait,fixed=TRUE),
                pLOD=gsub("$","",pLOD,fixed=TRUE),
                Model=gsub("$","",Model,fixed=TRUE)) %>%
        mutate(Trait=gsub(" ",'\\ ',Trait,fixed=TRUE),
                pLOD=gsub(" ",'\\ ',pLOD,fixed=TRUE),
                Model=gsub(" ",'\\ ',Model,fixed=TRUE)) %>%
        flextable() %>%
        bold(j=1, bold=TRUE, part="body") %>%
        align(j=c(1,2,7,8), align='left', part="body") %>%
        align(j=c(3,4,5,6), align='right', part="body") %>%
        mk_par(j = "Trait", part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        mk_par(j = "Model", part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        mk_par(j = "pLOD", part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        flextable::footnote(i = 1, 
                                j = ~ Trait + Model + pLOD,
                                value=as_paragraph(c(    "Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                                                        "Significance codes for Genotype Effects",
                                                        "Significance codes from QTL pvalues")),
                                ref_symbols=c("a","b","c"),
                                part = "header") %>%
        theme_booktabs() %>%
        hline(i=(which(qtl.collated.succinct.tb$Trait != "")-1)[-1], part="body") %>%
        set_table_properties(width=table.width)
        return(flxtable)
}

cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/'),mustWork=TRUE),'/qtl_collated.',qtl_scan_method,'.png'),"'")

save.image(paste0(workflow,"/.RData.10_01.genQTLtable.",qtl_scan_method,".",num_top_qtls))
load(paste0(workflow,"/.RData.10_01.genQTLtable.",qtl_scan_method,".",num_top_qtls))

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

#Defaults (can be overridden with command-line invocation)
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
qtl_method      <- get0("qtl_method", ifnotfound="stepwiseqtl") #Which method to filter
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top QTLs to show per trait

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
                                filter(method == qtl_method) %>%
                                mutate(marker = gsub(".+cM_(.+)", "\\1",nearest.marker), 
                                       model_name = model_to_name(trait.cfg.tb,model,trait), 
                                       trait_name = trait_to_name(trait.cfg.tb,model,trait)) %>%
                                inner_join(consensusMapWithGenes.filtered.tb, by="marker") %>%
                                arrange(trait_name,model_name,desc(marker.variance)) %>%
                                mutate(trait_name=paste0(trait_name,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$")) %>%
                                mutate(model_name=paste0(model_name,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) %>%
                                mutate(marker.variance=percent(marker.variance/100), 
                                       trait_repeat=(trait_name == c("",trait_name[-length(trait_name)])),
                                       model_repeat=(model_name == c("",model_name[-length(model_name)]))) %>%
                                group_by(trait_name,model_name) %>%
                                mutate(top_qtls = c(rep(TRUE,num_top_qtls),rep(FALSE,length(marker.variance)-num_top_qtls))) %>%
                                filter(top_qtls == TRUE)


qtl.collated.succinct.tb <- qtl.collated.filtered.tb %>%
                                mutate(marker = ifelse(is.na(marker)," ",marker)) %>%
                                mutate(blast = ifelse((!is.na(blast1) & !is.na(blast2)),paste0(blast1,";",blast2),ifelse(is.na(blast1),ifelse(is.na(blast2)," ", blast2),blast1))) %>%
                                mutate(marker.variance = color_bar("lightblue")(marker.variance),
                                       trait_name = ifelse(trait_repeat, "", cell_spec(trait_name, "html", bold=TRUE)),
                                       model_name = ifelse(model_repeat, "", cell_spec(model_name, "html", bold=TRUE)),
                                       position = paste0(round.digits(position,2),"\u00b1",round.digits(interval/2,2)),
                                       qtl.lod = paste0(round.digits(qtl.lod,2),"$^{",unlist(map(qtl.pvalue,siginfo)),"}$")
                                       ) %>%
                                select(trait_name,model_name,chr,position,qtl.lod,marker.variance,marker,blast,trait_repeat,model_repeat)

colnames(qtl.collated.succinct.tb) <- c("Trait[note]", "Model[note]", "Linkage Group", "Marker Location \u00b1 1.5LOD (cM)", "pLOD[note]", "Variance Explained by QTL", "Nearest Marker", "Putative Function", "trait_repeat", "model_repeat")

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
    kable("html", align='llrrrrll', escape = FALSE, table.attr="id=\"kableTable\"") %>%
    kable_paper("striped", full_width=FALSE) %>%
	row_spec(row=which(qtl.collated.succinct.tb[,"Trait[note]"] != ""), extra_css = "border-top: 1px solid #ddd") %>%
    column_spec(3, width = "1cm") %>%
    column_spec(4, width = "2.5cm") %>%
    column_spec(5, width = "2cm") %>%
    column_spec(8, width = "20cm") %>%
	add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
				   "Significance codes for Genotype Effects",
				   "Significance codes from QTL pvalues"))


cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/qtl_collated.',qtl_method,'.png'),mustWork=TRUE),"'"))

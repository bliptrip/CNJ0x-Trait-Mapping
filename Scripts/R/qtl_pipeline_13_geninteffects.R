#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(cowplot)
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


extract_int_effects <- function(model, trait, chr, position, chr2, position2) {
    #print(paste0("method = ",method, ", model = ", model, "trait = ", trait))
    qtl.mname <- paste0(chr,'@',round.digits(position,1))
    qtl2.mname <- paste0(chr2,'@',round.digits(position2,1))
    cross <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/cross.rds'))
    qtl   <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/',trait,'/scansw.rds'))
    #The qtl$prob contains the list of significant QTLs and the probability of a given genotype at the QTL.  I would like to show a boxplot of
    #blup values at the different genotypes for each QTL, but since the genotype is a mixed distribution at each QTL, I will only include genotypes with a higher
    #than, say, 95% probability of being a given genotype, and choose that as the representative genotype (assigning the entire BLUP and/or trait to that genotype for organization)
    #Determine which index corresponds to QTL
    qtl.i <- which(qtl$name == qtl.mname)
    qtl2.i <- which(qtl$name == qtl2.mname)
    qtl.p <- data.frame(qtl$prob[[qtl.i]])
    qtl2.p <- data.frame(qtl$prob[[qtl2.i]])
    qtl.p$id_i <- rownames(qtl.p)
    qtl2.p$id_i <- rownames(qtl2.p)
    qtl.p.genos <- qtl.p %>% 
                    mutate(blup=cross$pheno[,trait]) %>% 
                    pivot_longer(!c(id_i,blup),names_to='genotype_q1', values_to='probability') %>% 
                    filter(probability > 0.95) %>%
                    select(!probability)
    qtl2.p.genos <- qtl2.p %>%
                    pivot_longer(!c(id_i),names_to='genotype_q2', values_to='probability') %>% 
                    filter(probability > 0.95) %>%
                    select(!probability)
    qtl.both.p.genos <- qtl.p.genos %>%
                            inner_join(qtl2.p.genos, by="id_i") %>%
                            nest_by(genotype_q1,genotype_q2,.key="blups")
    return(qtl.both.p.genos)
}

#For each qtl in the collated file, use it's position and consensus position to calculate the effects.  Store this information in the collated file?
generate_collated_int_effects <- function(qtl.collated.tb) {
    qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                    arrange(method,trait,model,desc(marker_variance))
                                    
    effects.tb <- NULL
    for( i in 1:nrow(qtl.collated.filtered.tb) ) {
        e.tb <- qtl.collated.filtered.tb[i,]
        effs.tb <- bind_cols(e.tb,extract_int_effects(e.tb$model, e.tb$trait, e.tb$chr, e.tb$position, e.tb$chr2, e.tb$position2))
        if( is.null(effects.tb) ) {
            effects.tb = effs.tb
        } else {
            effects.tb <- bind_rows(effects.tb, effs.tb)
        }
    }
    return(effects.tb)
}
qtl.ints.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), col_names=TRUE) %>%
                        filter(method == "stepwiseqtl") %>%
                        filter(!is.na(chr2) & !is.na(position2)) %>% #Only include interaction QTLs for now
                        group_by(trait) %>%
                        mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>% #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                        mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                        ungroup()
colnames(qtl.ints.tb) <- gsub('.','_',colnames(qtl.ints.tb),fixed=TRUE)

effs.ints.tb <- generate_collated_int_effects(qtl.ints.tb) %>%
                    mutate(genotype_q1 = factor(genotype_q1,levels=c("AC","AD","BC","BD"),ordered=TRUE),
                           genotype_q2 = factor(genotype_q2,levels=c("AC","AD","BC","BD"),ordered=TRUE))

saveRDS(effs.ints.tb,file=paste0(workflow,'/traits/effects_interactions_collated.rds'), compress=TRUE)
effs.ints.tb <- readRDS(paste0(workflow,'/traits/effects_interactions_collated.rds')) #We can start here to load older state


effs.ints.tops.tb <- effs.ints.tb %>% 
                nest_by(method,trait,model,chr,position,chr2,position2,marker_variance) %>%
                group_by(method,trait,model) %>%
                arrange(method,trait,model,desc(marker_variance)) %>%
                mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker_variance))),rep(FALSE,length(marker_variance)-min(num_top_qtls,length(marker_variance))))) %>%
                filter(top_qtls == TRUE) %>%
                select(-top_qtls) %>%
                ungroup() %>%
                unnest(data) %>%
                unnest(blups)


effs.ints.plots.tb <- effs.ints.tops.tb %>% 
                        group_by(method,trait,model) %>%
                        mutate(min_blup = min(blup), max_blup=max(blup)) %>%
                        ungroup() %>%
                        group_by(method, trait, model, chr, position, chr2, position2) %>%
                        mutate(genotype = factor(paste0(genotype_q1,"×",genotype_q2), levels=c("AC×AC","AC×AD","AC×BC","AC×BD","AD×AC","AD×AD","AD×BC","AD×BD","BC×AC","BC×AD","BC×BC","BC×BD","BD×AC","BD×AD","BD×BC","BD×BD"))) %>%
                        do(plot0  = ggplot(., aes(x=genotype_q2, y=blup, fill=genotype)) +
                                            geom_blank(aes(y=min_blup)) +
                                            geom_blank(aes(y=max_blup)) +
                                            geom_boxplot(alpha=0.5) +
                                            scale_x_discrete(drop=FALSE) +
                                            geom_jitter(width=0.1, alpha=0.3) +
                                            geom_hline(mapping=aes(yintercept=0),color="grey30",linetype="dashed") +
                                            guides(fill = 'none') +
                                            coord_flip() +
                                            facet_grid(cols=vars(genotype_q1)) +
                                            theme(axis.title.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  axis.line   = element_line(color="black"),
                                                  axis.ticks = element_line(color="black"),
                                                  panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                  panel.background = element_rect(fill="transparent"),
                                                  plot.background = element_rect(fill="transparent")),
                           plot1 = ggplot(., aes(x=1, y=blup, group=genotype_q1, fill=genotype_q1)) +
                                            geom_blank(aes(y=min_blup)) +
                                            geom_blank(aes(y=max_blup)) +
                                            geom_boxplot(alpha=0.5) +
                                            scale_x_discrete(drop=FALSE) +
                                            geom_jitter(width=0.1, alpha=0.3) +
                                            geom_hline(mapping=aes(yintercept=0),color="grey30",linetype="dashed") +
                                            guides(fill = 'none') +
                                            coord_flip() +
                                            facet_grid(cols=vars(genotype_q1)) +
                                            theme(axis.title.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  axis.line   = element_line(color="black"),
                                                  axis.ticks = element_line(color="black"),
                                                  panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                  panel.background = element_rect(fill="transparent"),
                                                  plot.background = element_rect(fill="transparent")),
                           plot2 = ggplot(., aes(x=genotype_q2, y=blup, fill=genotype_q2)) +
                                            geom_blank(aes(y=min_blup)) +
                                            geom_blank(aes(y=max_blup)) +
                                            geom_boxplot(alpha=0.5) +
                                            scale_x_discrete(drop=FALSE) +
                                            geom_jitter(width=0.1, alpha=0.3) +
                                            geom_hline(mapping=aes(yintercept=0),color="grey30",linetype="dashed") +
                                            guides(fill = 'none') +
                                            coord_flip() +
                                            theme(axis.text.y = element_blank(),
                                                  axis.title.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  axis.line   = element_line(color="black"),
                                                  axis.ticks = element_line(color="black"),
                                                  panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                  panel.background = element_rect(fill="transparent"),
                                                  plot.background = element_rect(fill="transparent"))) %>%
                        mutate(plot_filename = paste0(normalizePath(paste0(workflow,'/traits/',model,'--',trait,'/',trait),mustWork=TRUE), '/effects_interactions_plot.blups.chr',chr,'_',round.digits(position,2),'cM--chr',chr2,'_',round.digits(position2,2),'cM.png')) %>%
                        group_walk(~ {
                                        plot = plot_grid(.x$plot0[[1]],.x$plot2[[1]],.x$plot1[[1]],nrow=2,ncol=2,rel_widths=c(1,.25))
                                        print(paste0("Saving ",.x$plot_filename))
                                        #png(filename=.x$plot_filename, width=640, height=320, bg="white")
                                        ggsave(filename=.x$plot_filename, plot = plot, device="png", bg="transparent", dpi=300, width=20, height=20, units="cm")
                                    })

effs.ints.all.tb <- effs.ints.plots.tb %>%
                        left_join(qtl.ints.tb, by=c("method","trait","model","chr","position","chr2","position2")) %>%
                        mutate(model_name = model_to_name(trait.cfg.tb,model,trait), 
                               trait_name = trait_to_name(trait.cfg.tb,model,trait),
                               position = round.digits(position,2),
                               position2 = round.digits(position2,2),
                               qtl_lod = paste0(round.digits(qtl_lod,2),"$^{",unlist(map(qtl_pvalue,siginfo)),"}$"))


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
generateReducedTable(effs.ints.all.tb)
generateReducedTable <- function(tbl, caption=NULL) {
    etbl1 <- tbl %>%
        generateTableRemoveRepeats() %>%
        select(trait_name,chr,position,chr2,position2,qtl_lod,marker_variance) %>%
        mutate(marker_variance=color_bar("lightblue")(percent(marker_variance/100))) %>%
        rename("Trait"=trait_name,
               "LG1"=chr,
               "Position1 (cM)"=position,
               "LG2"=chr,
               "Position2 (cM)"=position,
               "pLOD"=qtl_lod,
               "Variance Explained by QTL"=marker_variance) %>%
        mutate("Interaction Effect Boxplots[note]"="")
    if( is_html_output() ) {
        etbl2 <- etbl1 %>% kable("html", caption=caption, align='lrrrrrrc', escape = FALSE, table.attr="id=\"kableTable\"")
    } else {
        etbl2 <- etbl1 %>% kable(align='lrrrrrrc', caption=caption, escape = FALSE)
    }
    etbl3  <- etbl2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        column_spec(1, bold=TRUE) %>%
        column_spec(2, width = "0.5cm") %>%
        column_spec(3, width = "1.5cm") %>%
        column_spec(4, width = "0.5cm") %>%
        column_spec(5, width = "1.5cm") %>%
        column_spec(6, width = "1.5cm") %>%
        column_spec(7, width = "1cm") %>%
        column_spec(8, width = "20cm", image=spec_image(tbl$plot_filename,2560,2560)) %>%
        add_footnote(c("Fill in"))
    return(etbl3)
}

save.image(paste0(workflow,"/.RData.13.geninteffects",num_top_qtls))
load(paste0(workflow,"/.RData.10_02.geninteffects",num_top_qtls))

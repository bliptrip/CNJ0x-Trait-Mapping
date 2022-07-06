#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(cowplot)
library(flextable)
library(formattable)
library(knitr)
library(RColorBrewer)
library(tidyverse)


#Defaults (can be overridden with command-line invocation)
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top interaction QTLs to show per trait
table.width     <- get0("table.width", ifnotfound=0.85)
collate_effects <- get0("collate_effects", ifnotfound=TRUE)
reload_table_functions <- get0("reload_table_functions", ifnotfound=FALSE)

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
num_top_qtls <- as.numeric(num_top_qtls) #Convert to numeric if passed in through command-line
collate_effects <- as.logical(collate_effects) #Convert to logical if passed through command-line
reload_table_functions=as.logical(reload_table_functions)

if( !reload_table_functions ) {
    source('./usefulFunctions.R')
    source(paste0(workflow,"/configs/model.cfg"))
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

#Non-interactive QTLs
    qtl.nints.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), col_names=TRUE) %>%
                    filter(method == "stepwiseqtl") %>%
                    filter(is.na(chr2) & is.na(position2)) %>%
                    select(method,model,trait,chr,position,marker.variance)
    colnames(qtl.nints.tb) <- gsub('.','_',colnames(qtl.nints.tb),fixed=TRUE)
    qtl.nints.tb <- unique(qtl.nints.tb) #Needed in case of duplicates

#Interactive QTLs
    qtl.ints.tb <- read_csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), col_names=TRUE) %>%
                            filter(method == "stepwiseqtl") %>%
                            filter(!is.na(chr2) & !is.na(position2)) %>% #Only include interaction QTLs for now
                            group_by(trait) %>%
                            mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue)), #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                                GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                            ungroup()
    colnames(qtl.ints.tb) <- gsub('.','_',colnames(qtl.ints.tb),fixed=TRUE)
    qtl.ints.tb <- unique(qtl.ints.tb) #Needed in case of duplicates

    if( collate_effects ) {
        effs.ints.tb <- generate_collated_int_effects(qtl.ints.tb) %>%
                            mutate(genotype_q1 = factor(genotype_q1,levels=c("AC","AD","BC","BD"),ordered=TRUE),
                                genotype_q2 = factor(genotype_q2,levels=c("AC","AD","BC","BD"),ordered=TRUE)) %>%
                            nest_by(method,trait,model,chr,position,chr2,position2,marker_variance) %>%
                            left_join(qtl.nints.tb, by=c("method","trait","model","chr","position"), suffix=c("","1")) %>%
                            left_join(qtl.nints.tb %>% rename(chr2=chr,position2=position), by=c("method","trait","model","chr2","position2"), suffix=c("","2")) %>%
                            unnest(data) %>%
                            ungroup()

        saveRDS(effs.ints.tb,file=paste0(workflow,'/traits/effects_interactions_collated.rds'), compress=TRUE)
    } else {
        effs.ints.tb <- readRDS(paste0(workflow,'/traits/effects_interactions_collated.rds')) #We can start here to load older state
    }

    effs.ints.tops.tb <- effs.ints.tb %>% 
                    nest_by(method,trait,model,chr,position,chr2,position2,marker_variance) %>%
                    group_by(method,trait,model) %>%
                    arrange(method,trait,model,desc(marker_variance)) %>%
                    mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker_variance))),rep(FALSE,length(marker_variance)-min(num_top_qtls,length(marker_variance))))) %>%
                    filter(top_qtls == TRUE) %>%
                    select(-top_qtls) %>%
                    ungroup() %>%
                    unnest(data) %>%
                    unnest(blups) %>%
                    group_by(method,trait,model) %>%
                    mutate(min_blup = min(blup), max_blup=max(blup)) %>%
                    ungroup() %>%
                    group_by(method, trait, model, chr, position, chr2, position2, genotype_q1) %>%
                    mutate(mean_q1 = mean(blup), sd_q1=sd(blup)) %>%
                    ungroup() %>%
                    group_by(method, trait, model, chr, position, chr2, position2, genotype_q2) %>%
                    mutate(mean_q2 = mean(blup), sd_q2=sd(blup)) %>%
                    ungroup() %>%
                    group_by(method, trait, model, chr, position, chr2, position2, genotype_q1, genotype_q2) %>%
                    mutate(mean_q12 = mean(blup), sd_q12=sd(blup)) %>%
                    ungroup() %>%
                    mutate(mean_diff = (mean_q1+mean_q2)-mean_q12,
                        mean_sum = (mean_q1+mean_q2),
                        genotype = factor(paste0(genotype_q1,"×",genotype_q2), levels=c("AC×AC","AC×AD","AC×BC","AC×BD","AD×AC","AD×AD","AD×BC","AD×BD","BC×AC","BC×AD","BC×BC","BC×BD","BD×AC","BD×AD","BD×BC","BD×BD")))

    effs.ints.plots.tb <- effs.ints.tops.tb %>% 
                            group_by(method, trait, model, chr, position, chr2, position2) %>%
                            do(mad=max(abs(.$mean_diff)),
                            nmad = max(abs(.$mean_diff))/(.$max_blup[1]-.$min_blup[1]), #Normalized mean absolute difference
                            plot0  = ggplot(., aes(y=genotype_q2, x=blup, fill=genotype)) +
                                                geom_boxplot(alpha=0.5) +
                                                geom_jitter(width=0.1, alpha=0.3) +
                                                geom_point(aes(x=mean_sum), shape="diamond", size=3, alpha=0.7) +
                                                geom_point(aes(x=mean_q12), shape="circle", size=3, alpha=0.7) +
                                                geom_blank(aes(x=min_blup)) +
                                                geom_blank(aes(x=max_blup)) +
                                                scale_y_discrete(drop=FALSE) +
                                                geom_vline(mapping=aes(xintercept=0),color="grey30",linetype="dashed") +
                                                guides(fill = 'none') +
                                                facet_grid(cols=vars(genotype_q1)) +
                                                theme(axis.title.x = element_blank(),
                                                    axis.title.y = element_blank(),
                                                    axis.line.x   = element_line(color="black"),
                                                    axis.line.y   = element_blank(),
                                                    axis.ticks.x = element_line(color="black"),
                                                    axis.ticks.y = element_blank(),
                                                    axis.text.x = element_text(size=12,angle=75,vjust=-.005),
                                                    axis.text.y = element_blank(),
                                                    strip.text = element_blank(),
                                                    strip.background = element_blank(),
                                                    panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                    panel.background = element_rect(fill="transparent"),
                                                    plot.background = element_rect(fill="transparent"),
                                                    plot.margin = margin(20,10,20,20)),
                            plot1 = ggplot(., aes(x=blup, y=1, group=genotype_q1, fill=genotype_q1)) +
                                                geom_blank(aes(x=min_blup)) +
                                                geom_blank(aes(x=max_blup)) +
                                                geom_boxplot(alpha=0.5) +
                                                scale_y_discrete(drop=FALSE) +
                                                geom_jitter(width=0.1, alpha=0.3) +
                                                geom_vline(mapping=aes(xintercept=0),color="grey30",linetype="dashed") +
                                                guides(fill = 'none') +
                                                facet_grid(cols=vars(genotype_q1)) +
                                                theme(axis.title.x = element_blank(),
                                                    axis.title.y = element_blank(),
                                                    axis.line.x   = element_line(color="black"),
                                                    axis.line.y   = element_blank(),
                                                    axis.text.x = element_text(size=12,angle=75,vjust=-.005),
                                                    axis.ticks.x = element_line(color="black"),
                                                    axis.ticks.y = element_blank(),
                                                    strip.text = element_text(face="bold",size=18),
                                                    strip.background = element_rect(fill="transparent"),
                                                    panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                    panel.background = element_rect(fill="transparent"),
                                                    plot.background = element_rect(fill="transparent"),
                                                    plot.margin = margin(10,10,20,24)),
                            plot2 = ggplot(., aes(y=genotype_q2, x=blup, fill=genotype_q2)) +
                                                geom_blank(aes(x=min_blup)) +
                                                geom_blank(aes(x=max_blup)) +
                                                geom_boxplot(alpha=0.5) +
                                                scale_y_discrete(drop=FALSE) +
                                                geom_jitter(width=0.1, alpha=0.3) +
                                                geom_vline(mapping=aes(xintercept=0),color="grey30",linetype="dashed") +
                                                guides(fill = 'none') +
                                                theme(axis.title.x = element_blank(),
                                                    axis.title.y = element_blank(),
                                                    axis.text.x = element_text(size=12,angle=75,vjust=-.005),
                                                    axis.text.y = element_text(face="bold",size=18),
                                                    axis.line   = element_line(color="black"),
                                                    axis.ticks = element_line(color="black"),
                                                    panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                    panel.background = element_rect(fill="transparent"),
                                                    plot.background = element_rect(fill="transparent"),
                                                    plot.margin = margin(20,22,24,20))) %>%
                            mutate(mad=unlist(mad),
                                nmad=unlist(nmad),
                                plot_filename = paste0(normalizePath(paste0(workflow,'/traits/',model,'--',trait,'/',trait),mustWork=TRUE), '/effects_interactions_plot.blups.chr',chr,'_',round.digits(position,2),'cM--chr',chr2,'_',round.digits(position2,2),'cM.png')) %>%
                            group_walk(~ {
                                            plot = plot_grid(.x$plot0[[1]],.x$plot2[[1]],.x$plot1[[1]],labels="auto",label_size=14, nrow=2,ncol=2,rel_widths=c(1,.4),rel_heights=c(1,.85))
                                            print(paste0("Saving ",.x$plot_filename))
                                            #png(filename=.x$plot_filename, width=640, height=320, bg="white")
                                            ggsave(filename=.x$plot_filename, plot = plot, device="png", bg="transparent", dpi=300, width=25, height=10, units="cm")
                                        })

    effs.ints.all.tb <- effs.ints.plots.tb %>%
                            left_join(qtl.ints.tb, by=c("method","trait","model","chr","position","chr2","position2")) %>%
                            left_join(qtl.nints.tb, by=c("method","trait","model","chr","position"), suffix=c("","1")) %>%
                            left_join(qtl.nints.tb %>% rename(chr2=chr,position2=position), by=c("method","trait","model","chr2","position2"), suffix=c("","2")) %>%
                            mutate(model_name = model_to_name(trait.cfg.tb,model,trait), 
                                trait_name = trait_to_name(trait.cfg.tb,model,trait),
                                position = round.digits(position,2),
                                position2 = round.digits(position2,2),
                                qtl_lod = paste0(round.digits(qtl_lod,2),"$^{",unlist(map(qtl_pvalue,siginfo)),"}$"))
} else {
    #Only reloading functions.  Load previous state.
    load(paste0(workflow,"/.RData.13_geninteffects.",num_top_qtls))
}

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

generateReducedTable <- function(tbl, caption=NULL, tfont_size=10) {
    etbl1 <- tbl %>%
        generateTableRemoveRepeats() %>%
        select(trait_name,chr,position,chr2,position2,qtl_lod,mad,nmad,marker_variance,marker_variance1,marker_variance2,plot_filename) %>%
        mutate(mad=round.digits(mad,3),
               nmad=percent(round.digits(nmad,4)),
               marker_variance=percent(marker_variance/100),
               marker_variance1=percent(marker_variance1/100),
               marker_variance2=percent(marker_variance2/100))
    if( is_html_output() ) {
        etbl1 <- etbl1 %>%
            mutate(nmad=color_bar(spec_color(nmad,option="E",alpha=0.3,scale_from=c(0,20)))(nmad),
                   marker_variance=color_bar(spec_color(marker_variance,option="E",alpha=0.3,scale_from=c(0,20)))(marker_variance),
                   marker_variance1=color_bar(spec_color(marker_variance1,option="E",alpha=0.3,scale_from=c(0,20)))(marker_variance1),
                   marker_variance2=color_bar(spec_color(marker_variance2,option="E",alpha=0.3,scale_from=c(0,20)))(marker_variance2))
    } else if ( is_pdf_output() ) {
        etbl1 <- etbl1 %>%
            mutate(nmad = knitr:::escape_latex(nmad),
                   marker_variance = knitr:::escape_latex(marker_variance),
                   marker_variance1 = knitr:::escape_latex(marker_variance1),
                   marker_variance2 = knitr:::escape_latex(marker_variance2))
    }
    etbl2 <- etbl1 %>%
        select(!plot_filename) %>%
        rename("Trait"=trait_name,
               "LG1"=chr,
               "Position1 (cM)"=position,
               "LG2"=chr2,
               "Position2 (cM)"=position2,
               "pLOD[note]"=qtl_lod,
               "MAD[note]"=mad,
               "NMAD[note]"=nmad,
               "Int[note] Variance"=marker_variance,
               "QTL1 Variance"=marker_variance1,
               "QTL2 Variance"=marker_variance2) %>%
        mutate("Interaction Effect Boxplots[note]"="") %>%
        kable(align='lrrrrrrrrrrc', caption=caption, escape = FALSE, longtable=TRUE, booktabs=TRUE)
    etbl3  <- etbl2 %>%
        kable_paper("striped", full_width=FALSE) %>%
        column_spec(1, width = "1cm", bold=TRUE) %>%
        column_spec(2, width = "0.5cm") %>%
        column_spec(3, width = "0.75cm") %>%
        column_spec(4, width = "0.5cm") %>%
        column_spec(5, width = "0.75cm") %>%
        column_spec(6, width = "0.75cm") %>%
        column_spec(7, width = ".5cm") %>%
        column_spec(8, width = ".5cm") %>%
        column_spec(9, width = ".75cm") %>%
        column_spec(10, width = ".75cm") %>%
        column_spec(11, width = ".75cm") %>%
        column_spec(12, width = "6cm", image=spec_image(tbl$plot_filename,width=708,height=283,res=300),
                        popover=paste0("<img src='",tbl$plot_filename,"' width='1024' height='410'>")) %>%
        add_footnote(c("QTL Penalized LOD Score w/ significance codes:\n*** pvalue $≥$ 0 and pvalue<0.001\n**  pvalue $≥$ 0.001 and pvalue<0.01\n*   pvalue $≥$ 0.01 and pvalue<0.05\n.  pvalue $≥$ 0.05 and pvalue<0.01\nNS  Not Significant\n",
                       "Max Average Difference - Maximum of all differences between mean interaction effects and and expected mean effect under additive models only.",
                       "Normalized Max Average Difference - Normalized value of Max Average Difference, using BLUP range as a normalization factor.",
                       "Percent of model variance explained by interaction effect.",
                       "Interaction effects are shown in subplot *a*, with boxplots showing the distribution of trait BLUPs organized by QTL 1 genotype (x-axis) and 
                       QTL 2 genotype (y-axis).  Subplot *c* is the marginal BLUP distributions factored by QTL genotype 1 (x-axis).
                       Subplot *b* is the marginal BLUP distributions factored by QTL genotype 2 (y-axis)."),
						escape=FALSE)
    if( is_html_output() ) {
        #This is necessary b/c some component of knitr kable is escaping some html that I don't want escaped, despite using the flag 'escape=FALSE'
        etbl3 <- gsub("&lt;","<",etbl3,fixed=T)
        etbl3 <- gsub("&gt;",">",etbl3,fixed=T)
        etbl3 <- gsub("&quot;","\"",etbl3,fixed=T)
        etbl3 <- gsub("data-placement=\"right\"","data-placement=\"auto\"",etbl3,fixed=T)
    } else if( is_pdf_output() ) {
        etbl3 <- etbl3 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }
    return(etbl3)
}
#generateReducedTable(effs.ints.all.tb %>% filter((model=='all-years') & (trait=="berry_length")))

save.image(paste0(workflow,"/.RData.13_geninteffects.",num_top_qtls))

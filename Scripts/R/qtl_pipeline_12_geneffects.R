library(formattable)
library(ggfittext)
library(jsonlite)
library(kableExtra)
library(knitr)
library(qtl)
library(tidyverse)

#Defaults
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top QTLs to show per trait

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
num_top_qtls <- as.numeric(num_top_qtls) #Convert to numeric

source('./usefulFunctions.R')
source(paste0(workflow,"/configs/model.cfg"))
trait.cfg.tb    <- read_csv(file=paste0(workflow,'/configs/model-traits.cfg.csv'), col_names=TRUE)

extract_effects <- function(method, model, trait, chr, position) {
    #print(paste0("method = ",method, ", model = ", model, "trait = ", trait))
    qtl.mname <- paste0(chr,'@',round.digits(position,1))
    cross <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/cross.rds'))
    qtl   <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/',trait,'/', ifelse(method == 'scanone', 'scanone.qtl', 'scansw'), '.rds'))
    effects <- effectplot(cross, pheno.col=trait, mname1=qtl.mname, draw=FALSE)
    re <- "(.+)\\.((AC)|(BC)|(AD)|(BD))$"
    marker <- unique(gsub(re, "\\1", names(effects$Means), fixed=FALSE))
    cm <- gsub(re, "\\2", names(effects$Means), fixed=FALSE)
    names(effects$Means) <- cm
    names(effects$SEs) <- cm
    AvB <- as.numeric((effects$Means["AC"] + effects$Means["AD"]) - (effects$Means["BC"] + effects$Means["BD"]))
    CvD <- as.numeric((effects$Means["AC"] + effects$Means["BC"]) - (effects$Means["AD"] + effects$Means["BD"]))
    Int <- as.numeric((effects$Means["AC"] + effects$Means["BD"]) - (effects$Means["AD"] + effects$Means["BC"]))
    effects.means.tb <- tibble(genotype=names(effects$Means), effect_mean=as.numeric(effects$Means))
    effects.ses.tb <- tibble(genotype=names(effects$SEs), effect_se=as.numeric(effects$SEs))
    #The qtl$prob contains the list of significant QTLs and the probability of a given genotype at the QTL.  I would like to show a boxplot of
    #blup values at the different genotypes for each QTL, but since the genotype is a mixed distribution at each QTL, I will only include genotypes with a higher
    #than, say, 95% probability of being a given genotype, and choose that as the representative genotype (assigning the entire BLUP and/or trait to that genotype for organization)
    #Determine which index corresponds to QTL
    qtl.i <- which(qtl$name == qtl.mname)
    #print(paste0("qtl.mname = ", qtl.mname, ", qtl.i = ",qtl.i))
    #print(paste0("qtl$name = ", qtl$name))
    qtl.p <- data.frame(qtl$prob[[qtl.i]])
    qtl.p$id_i <- rownames(qtl.p)
    qtl.p.nest <- qtl.p %>% 
                    mutate(blup=cross$pheno[,trait]) %>% 
                    pivot_longer(!c(id_i,blup),names_to='genotype', values_to='probability') %>% 
                    filter(probability > 0.95) %>%
                    group_by(genotype) %>%
                    nest() %>%
                    left_join(effects.means.tb, by="genotype") %>%
                    left_join(effects.ses.tb, by="genotype") %>%
                    mutate(method=method, model=model, trait=trait, chr=chr, position=position, AvB=AvB, CvD=CvD, Int=Int, blups=data) %>%
                    select(-data) %>%
                    group_by(method,model,trait,chr,position) %>%
                    nest()
    return(qtl.p.nest)
}

#For each qtl in the collated file, use it's position and consensus position to calculate the effects.  Store this information in the collated file?
generate_collated_effects <- function(qtl.collated.tb,num_top_qtls) {
    qtl.collated.filtered.tb <- qtl.collated.tb %>%
                                    arrange(method,trait,model,desc(marker_variance)) %>%
                                    group_by(method,trait,model)
                                    
    effects.tb <- NULL
    for( i in 1:nrow(qtl.collated.filtered.tb) ) {
        e.tb <- qtl.collated.filtered.tb[i,]
        effs.tb <- extract_effects(e.tb$method, e.tb$model, e.tb$trait, e.tb$chr, e.tb$position)
        if( is.null(effects.tb) ) {
            effects.tb = effs.tb
        } else {
            effects.tb <- rbind(effects.tb, effs.tb)
        }
    }
    effects.collated.tb  <- qtl.collated.filtered.tb %>%
                            left_join(effects.tb, by=c("method","trait","model","chr","position"))
    return(effects.collated.tb)
}


qtl.tb  <- read_csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), col_names=TRUE) %>%
                    filter(is.na(chr2) & is.na(position2)) %>% #Remove interaction effects for now
                    group_by(trait) %>%
                    mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>%
                    mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                    ungroup()
colnames(qtl.tb) <- gsub('.','_',colnames(qtl.tb),fixed=TRUE)

effs.tb <- generate_collated_effects(qtl.tb, num_top_qtls)
colnames(effs.tb) <- gsub('.','_',colnames(effs.tb),fixed=TRUE)
saveRDS(effs.tb,file=paste0(workflow,'/traits/effects_collated.rds'), compress=TRUE)
write_json(effs.tb, paste0(workflow, "/traits/effects_collated.json"), auto_unbox=T, pretty=T)
#Write a csv file without the grouped blups
write.csv(effs.tb %>% unnest(data) %>% select(!blups), file=paste0(workflow,'/traits/effects_collated.csv'), row.names=FALSE)
effs.tb <- readRDS(paste0(workflow,'/traits/effects_collated.rds')) #We can start here to load older state

#Plot generation
effs.1.tb <- effs.tb %>% 
                group_by(method,trait,model) %>%
                mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker_variance))),rep(FALSE,length(marker_variance)-min(num_top_qtls,length(marker_variance))))) %>%
                filter(top_qtls == TRUE) %>%
                select(-top_qtls) %>%
                unnest(data) %>%
                unnest(blups)

effs.filtered1.tb <- effs.1.tb %>% 
                        group_by(method,trait) %>%
                            mutate(min_blup = min(blup), max_blup=max(blup)) %>%
                            ungroup() %>%
                            group_by(method, trait, model, chr, position) %>%
                            mutate(genotype = factor(levels=rev(c("AC","AD","BC","BD")),genotype)) %>%
                            do(plot  = ggplot(., aes(x=factor(genotype), y=blup, fill=factor(genotype))) +
                                                geom_blank(aes(y=min_blup)) +
                                                geom_blank(aes(y=max_blup)) +
                                                geom_boxplot(alpha=0.5) +
                                                geom_jitter(width=0.1, alpha=0.3) +
                                                geom_hline(mapping=aes(yintercept=0),color="grey30",linetype="dashed") +
                                                coord_flip() +
                                                guides(fill = 'none') +
                                                theme(axis.text.x = element_text(face="bold", size=12, angle = 60, hjust = 1),
                                                      axis.text.y = element_text(face="bold", size=12),
                                                      axis.title.x = element_blank(),
                                                      axis.title.y = element_blank(),
                                                      axis.line   = element_line(color="black"),
                                                      axis.ticks = element_line(color="black"),
                                                      panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                      panel.background = element_rect(fill="transparent"),
                                                      plot.background = element_rect(fill="transparent"))) %>%
                            mutate(plot_filename = paste0(normalizePath(paste0(workflow,'/traits/',model,'--',trait,'/',trait),mustWork=TRUE), '/effects_plot.blups.chr',chr,'_',round.digits(position,2),'cm.',method,'.png')) %>%
                            group_walk(~ {
                                            print(paste0("Saving ",.x$plot_filename))
                                            #png(filename=.x$plot_filename, width=640, height=320, bg="white")
                                            ggsave(filename=.x$plot_filename, plot = .x$plot[[1]], device="png", bg="transparent", dpi=300, width=20, height=5, units="cm")
                                        })
effs.2.tb <- effs.tb %>% 
                group_by(method,trait,model) %>%
                mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker_variance))),rep(FALSE,length(marker_variance)-min(num_top_qtls,length(marker_variance))))) %>%
                filter(top_qtls == TRUE) %>%
                select(-top_qtls) %>%
                unnest(data) %>%
                pivot_wider(names_from=genotype,values_from=c(effect_mean,effect_se,blups)) %>%
                pivot_longer(c(AvB,CvD,Int), names_to="effect_type",values_to="effect_value")

min_effect_value <- floor(min(effs.2.tb$effect_value))
max_effect_value <- ceiling(max(effs.2.tb$effect_value))

effs.filtered2.tb <- effs.2.tb %>%
                            group_by(method,trait) %>%
                            mutate(min_effect_value = min(effect_value), max_effect_value = max(effect_value)) %>%
                            ungroup() %>%
                            mutate(effect.label = gsub("AvB","A.-B.",effect_type)) %>%
                            mutate(effect.label = gsub("CvD",".C-.D",effect.label)) %>%
                            mutate(effect.label = factor(levels=rev(c("A.-B.",".C-.D","Int")), effect.label)) %>%
                            group_by(method, trait, model, chr, position) %>%
                            do(plot  = ggplot(.,aes(x=effect.label, y=effect_value, fill=effect.label)) +
                                                geom_blank(aes(y=min_effect_value)) +
                                                geom_blank(aes(y=max_effect_value)) +
                                                geom_col(alpha=0.5) +
                                                geom_bar_text(aes(label=round.digits(effect_value,2)),angle=25) +
                                                geom_hline(mapping=aes(yintercept=0),color="grey30",linetype="dashed") +
                                                coord_flip() +
                                                guides(fill = 'none') +
                                                theme(axis.text.x = element_text(face="bold", size=12, angle = 60, hjust = 1),
                                                      axis.text.y = element_text(face="bold", size=12),
                                                      axis.title.x = element_blank(),
                                                      axis.title.y = element_blank(),
                                                      axis.line   = element_line(color="black"),
                                                      axis.ticks = element_line(color="black"),
                                                      panel.grid = element_line(color="darkgray",size=0.25,linetype=3),
                                                      panel.background = element_rect(fill="transparent"),
                                                      plot.background = element_rect(fill="transparent"))) %>%
                            mutate(plot_mpieffects_filename = paste0(normalizePath(paste0(workflow,'/traits/',model,'--',trait,'/',trait), mustWork=TRUE),'/effects_plot.mpieffects.chr',chr,'_',round.digits(position,2),'cm.',method,'.png')) %>%
                            group_walk(~ {
                                            print(paste0("Saving ",.x$plot_mpieffects_filename))
                                            ggsave(filename=.x$plot_mpieffects_filename, plot = .x$plot[[1]], device="png", bg="transparent", dpi=300, width=10, height=5, units="cm")
                                        })

#Replace GxY interaction significance values with those for a trait's 'all-year' model.
qtl.tb  <- qtl.tb %>%
                        select(-chr2,-position2) %>%
                        group_by(trait) %>%
                        mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>%
                        mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                        ungroup()

#Table generation
effs.complete.tb <-  qtl.tb %>% 
                        arrange(method,trait,model,desc(marker_variance)) %>%
                        inner_join(effs.filtered1.tb, by=c("method","trait","model","chr","position"), suffix=c("qtl","effs")) %>%
                        inner_join(effs.filtered2.tb, by=c("method","trait","model","chr","position"), suffix=c("qtl","effs")) %>%
                        mutate(marker = gsub(".+cM_(.+)", "\\1",nearest_marker), 
                               trait_name = trait_to_name(trait.cfg.tb,model,trait), 
                               model_name = model_to_name(trait.cfg.tb,model,trait),
                               marker_variance=marker_variance/100,
                               model_variance=model_variance/100)
                            

generateTableSignifSymbols <- function(tbl) {
    return( tbl %>%
                mutate(trait_name=paste0(trait_name,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$"),
                       model_name=paste0(model_name,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) )
}

#arrange(trait_name,model_name,desc(marker.variance)) %>%
generateTableRemoveRepeats <- function(tbl) {
    return( tbl %>%
        mutate(trait_repeat=(trait_name == c("",trait_name[-length(trait_name)]))) %>%
        mutate(model_repeat=(trait_repeat & (model_name == c("",model_name[-length(model_name)])))) %>%
        mutate(trait_name = ifelse(trait_repeat, "",trait_name),
               model_name = ifelse(model_repeat, "",model_name)) )
}

generateReducedTable <- function(tbl, caption=NULL) {
    etbl1 <- tbl %>%
            generateTableRemoveRepeats() %>%
            mutate(model_variance = ifelse(model_repeat, "", color_bar("lightgreen")(percent(model_variance))),
                   marker = ifelse(is.na(marker)," ",marker),
                   position = paste0(round.digits(position,2),"±",round.digits(interval/2,2)),
                   marker_variance = color_bar("lightblue")(percent(marker_variance)))
    etbl2 <- etbl1 %>%
             select(trait_name,chr,position,qtl_lod,marker_variance,model_variance) %>%
             rename("Trait"=trait_name, 
                    "LG"=chr, 
                    "Marker Location[note]"=position, 
                    "Variance Explained by QTL"=marker_variance,
                    "Model Variance[note]"=model_variance,
                    "pLOD[note]"=qtl_lod) %>% 
             mutate('Effect Size Boxplots[note]'="") %>%
             mutate('Effect Difference Plots[note]'="")
    if( is_html_output() ) {
        etbl3 <- etbl2 %>% kable("html", caption=caption, table.attr="id=\"kableTable\"", align='lllrrrllcc', escape = FALSE)
    } else {
        etbl3 <- etbl2 %>% kable(caption=caption, align='lllrrrllcc', escape = FALSE)
    }
    etbl4 <- etbl3 %>%
                kable_paper("striped", full_width=TRUE) %>%
                column_spec(1, bold=TRUE, width = "1.5cm") %>%
                column_spec(2, width = "1cm") %>%
                column_spec(3, width = "2.5cm") %>%
                column_spec(4, width = "1.5cm") %>%
                column_spec(5, width = "2.5cm") %>%
                column_spec(6, width = "2.5cm") %>%
                column_spec(7, width = "5cm", image=spec_image(etbl1$plot_filename,1280,320)) %>%
                column_spec(8, width = "4cm", image=spec_image(etbl1$plot_mpieffects_filename,640,320)) %>%
                add_footnote(c("QTL location ± 1.5pLOD interval (cM)",
                        "Variance of model with all significant QTLs fitted.",
                        "QTL Penalized LOD Score w/ significance codes:\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n",
                        "Boxplots of nearest marker BLUPs grouped by genotypes.  Haplotypes A and B are from maternal parent P1, and haplotypes C and D are from paternal parent P2.",
                        "Effect differences for mean QTL effect size estimates for each progeny genotype.  A.-B. is the maternal effect, calculated as (AC+AD)-(BC+BD).  .C-.D is the paternal effect, calculated as (AC + BC) – (AD + BD).  Int is the interaction effect, calculated as (AC + BD)-(AD+BC) (Sewell et al., 2002)."))
    return(etbl4)
}

generateTable <- function(tbl, meths=c("scanone","stepwiseqtl"), caption=NULL) {
    etbl1 <- tbl %>%
            filter(method %in% meths) %>%
            generateTableSignifSymbols() %>%
            generateTableRemoveRepeats() %>%
            mutate(model_variance = ifelse(model_repeat, "", color_bar("lightgreen")(percent(model_variance))),
                   marker = ifelse(is.na(marker)," ",marker),
                   position = paste0(round.digits(position,2),"±",round.digits(interval/2,2)),
                   marker_variance = color_bar("lightblue")(percent(marker_variance)))
    etbl2 <- etbl1 %>%
             select(trait_name,model_name,model_variance,chr,position,qtl_lod,marker_variance) %>%
             rename("Trait[note]"=trait_name, 
                    "Model[note]"=model_name, 
                    "Model Variance[note]"=model_variance, 
                    "LG"=chr, 
                    "Marker Location[note]"=position, 
                    "pLOD[note]"=qtl_lod, 
                    "Variance Explained by QTL"=marker_variance) %>%
             mutate('Effect Size Boxplots[note]'="") %>%
             mutate('Effect Difference Plots[note]'="")
    if( is_html_output() ) {
        etbl3 <- etbl2 %>% kable("html", caption=caption, table.attr="id=\"kableTable\"", align='lllrrrllcc', escape = FALSE)
    } else {
        etbl3 <- etbl2 %>% kable(caption=caption, align='lllrrrllcc', escape = FALSE)
    }
    etbl4 <- etbl3 %>%
                row_spec(row=(which(etbl1$trait_name != "")-1)[-1], hline_after = TRUE) %>%
                kable_paper("striped", full_width=TRUE) %>%
                column_spec(1, width = "1.5cm") %>%
                column_spec(4, width = "1cm") %>%
                column_spec(5, width = "2.5cm") %>%
                column_spec(8, width = "5cm", image=spec_image(paste0('file://',etbl1$plot_filename),1280,320)) %>%
                column_spec(9, width = "4cm", image=spec_image(paste0('file://',etbl1$plot_mpieffects_filename),640,320)) %>%
                add_footnote(c(paste0("All QTLs in table derived from running R/qtl package function ",meths,"().\n","Significance codes for model genotype$*$year effects appended to trait:\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n"), 
                        "Significance codes for model genotype effects appended to model",
                        "Variance of model with all significant QTLs fitted.",
                        "QTL location ± 1.5pLOD interval (cM)",
                        "Penalized LOD Score w/ significance codes for QTL appended",
                        "Boxplots of nearest marker BLUPs grouped by genotypes.  Haplotypes A and B are from maternal parent P1, and haplotypes C and D are from paternal parent P2.",
                        "Effect differences for mean QTL effect size estimates for each progeny genotype.  A.-B. is the maternal effect, calculated as (AC+AD)-(BC+BD).  .C-.D is the paternal effect, calculated as (AC + BC) – (AD + BD).  Int is the interaction effect, calculated as (AC + BD)-(AD+BC) (Sewell et al., 2002)."))
    return(etbl4)
}

#Table for all traits
#ktable.one <- generateTable(effs.complete.tb, "scanone", "QTLs (scanone) with Effect Plots for All Traits")
#ktable.one
#cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/effectplots.scanone.png'),mustWork=TRUE),"'"))

#ktable.sw <- generateTable(effs.complete.tb, "stepwiseqtl", "QTLs (stepwiseqtl) with Effect Plots for All Traits")
#ktable.sw
#cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/effectplots.stepwiseqtl.png'),mustWork=TRUE),"'"))

#Break out into more meaninful nuggets
#Table per trait
#effs.complete.tb %>%
#    group_by(method,trait) %>%
#    group_walk( ~ {
#                    ktbl <- generateTable(cbind(.y,.x), .y$method, paste0("QTLs (", .y$method, ") with Effect Plots for trait ", .y$trait))
#                    print(ktbl)
#                  } )


save.image(paste0(workflow,"/.RData.12_geneffects.",num_top_qtls))
load(paste0(workflow,"/.RData.12_geneffects.",num_top_qtls))

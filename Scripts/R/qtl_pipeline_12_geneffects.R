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
collate_effects <- get0("collate_effects", ifnotfound=TRUE)
reload_table_functions <- get0("reload_table_functions", ifnotfound=FALSE)

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
num_top_qtls <- as.numeric(num_top_qtls) #Convert to numeric
reload_table_functions=as.logical(reload_table_functions)

if( !reload_table_functions ) {
    source('./usefulFunctions.R')
    source(paste0(workflow,"/configs/model.cfg"))
    trait.cfg.tb    <- read_csv(file=paste0(workflow,'/configs/model-traits.cfg.csv'), col_names=TRUE)

    extract_effects <- function(method, model, trait, chr, position) {
        elist <- extract_effects_primitive(workflow, method, model, trait, chr, position) 
        effects.means.tb <- elist$means
        effects.ses.tb <- elist$ses
        qtl.mname <- elist$qtl.mname
        cross <- elist$cross
        AvB <- elist$AvB #Maternal Effect
        CvD <- elist$CvD #Paternal Effect
        Int <- elist$Int #Interaction Effect
        #The qtl$prob contains the list of significant QTLs and the probability of a given genotype at the QTL.  I would like to show a boxplot of
        #blup values at the different genotypes for each QTL, but since the genotype is a mixed distribution at each QTL, I will only include genotypes with a higher
        #than, say, 95% probability of being a given genotype, and choose that as the representative genotype (assigning the entire BLUP and/or trait to that genotype for organization)
        #Determine which index corresponds to QTL
        qtl   <- readRDS(file=paste0(workflow,'/traits/',model,'--',trait,'/',trait,'/', ifelse(method == 'scanone', 'scanone.qtl', 'scansw'), '.rds'))
        qtl.i <- which(qtl$name == qtl.mname)
        print(paste0("qtl.mname = ", qtl.mname, ", qtl.i = ",qtl.i))
        print(paste0("qtl$name = ", qtl$name))
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
                        dplyr::select(-data) %>%
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
                        group_by(trait) %>%
                        mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>%
                        mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                        ungroup()
    colnames(qtl.tb) <- gsub('.','_',colnames(qtl.tb),fixed=TRUE)

    if( collate_effects ) {
        effs.tb <- generate_collated_effects(qtl.tb %>% filter(is.na(chr2) & is.na(position2)), num_top_qtls)
        colnames(effs.tb) <- gsub('.','_',colnames(effs.tb),fixed=TRUE)
        saveRDS(effs.tb,file=paste0(workflow,'/traits/effects_collated.rds'), compress=TRUE)
    } else {
        effs.tb <- readRDS(paste0(workflow,'/traits/effects_collated.rds')) #We can start here to load older state
    }

    effs.epistatics.tb <- qtl.tb %>% 
            filter(!is.na(chr2) & !is.na(position2)) %>% #Only include interaction effects
            mutate(data = NA) # Insert an empty data field so it can be merged wtih effs.tb -- interactions for now will not have effects made available

    effs.both.tb <- bind_rows(effs.tb, effs.epistatics.tb)
    write_json(effs.both.tb, paste0(workflow, "/traits/effects_collated.json"), auto_unbox=T, pretty=T)
    #Write a csv file without the grouped blups
    write.csv(effs.tb %>% unnest(data) %>% dplyr::select(!blups), file=paste0(workflow,'/traits/effects_collated.csv'), row.names=FALSE)

    #Plot generation
    effs.1.tb <- effs.tb %>% 
                    group_by(method,trait,model) %>%
                    mutate(top_qtls = c(rep(TRUE,min(num_top_qtls,length(marker_variance))),rep(FALSE,length(marker_variance)-min(num_top_qtls,length(marker_variance))))) %>%
                    filter(top_qtls == TRUE) %>%
                    dplyr::select(-top_qtls) %>%
                    unnest(data) %>%
                    unnest(blups)

    effs.1b.tb <- effs.1.tb %>% 
                    group_by(method,trait) %>%
                    mutate(min_blup = min(blup), max_blup=max(blup)) %>%
                    ungroup() %>%
                    group_by(method, trait, model, chr, position) %>%
                    mutate(genotype = factor(levels=rev(c("AC","AD","BC","BD")),genotype))

    effs.filtered1.tb <- effs.1b.tb %>%
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
                    dplyr::select(-top_qtls) %>%
                    unnest(data) %>%
                    pivot_wider(names_from=genotype,values_from=c(effect_mean,effect_se,blups)) %>%
                    pivot_longer(c(AvB,CvD,Int), names_to="effect_type",values_to="effect_value")

    min_effect_value <- floor(min(effs.2.tb$effect_value))
    max_effect_value <- ceiling(max(effs.2.tb$effect_value))

    #Filter out only columns desired for final effs.complete.tb
    effs.2a.tb <- effs.2.tb %>% 
                    pivot_wider(names_from=effect_type, values_from=effect_value) %>%
                    dplyr::select(method,trait,model,chr,position,effect_mean_AC,effect_mean_AD,effect_mean_BC,effect_mean_BD,effect_se_AC,effect_se_AD,effect_se_BC,effect_se_BD, AvB, CvD, Int)


    effs.2b.tb <- effs.2.tb %>%
                                group_by(method,trait) %>%
                                mutate(min_effect_value = min(effect_value), max_effect_value = max(effect_value)) %>%
                                ungroup() %>%
                                mutate(effect.label = gsub("AvB","A.-B.",effect_type)) %>%
                                mutate(effect.label = gsub("CvD",".C-.D",effect.label)) %>%
                                mutate(effect.label = factor(levels=rev(c("A.-B.",".C-.D","Int")), effect.label))


    effs.filtered2.tb <- effs.2b.tb %>%
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
                                mutate(plot_mpieffects_filename = paste0(normalizePath(paste0(workflow,'/traits/',model,'--',trait,'/',trait), mustWork=TRUE),'/effects_plot.mpieffects.chr',chr,'_',round.digits(position,2),'cm.',method,'.apng')) %>%
                                group_walk(~ {
                                                print(paste0("Saving ",.x$plot_mpieffects_filename))
                                                ggsave(filename=.x$plot_mpieffects_filename, plot = .x$plot[[1]], device="png", bg="transparent", dpi=300, width=10, height=5, units="cm")
                                            })

    #Replace GxY interaction significance values with those for a trait's 'all-year' model.
    qtl.tb  <- qtl.tb %>%
                            dplyr::select(-chr2,-position2) %>%
                            group_by(trait) %>%
                            mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>%
                            mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                            ungroup()

    #Table generation
    effs.complete.tb <-  qtl.tb %>% 
                            arrange(method,trait,model,desc(marker_variance)) %>%
                            inner_join(effs.filtered1.tb, by=c("method","trait","model","chr","position"), suffix=c("qtl","effs")) %>%
                            inner_join(effs.filtered2.tb, by=c("method","trait","model","chr","position"), suffix=c("qtl","effs")) %>%
                            inner_join(effs.2a.tb, by=c("method","trait","model","chr","position")) %>%
                            mutate(marker = gsub(".+cM_(.+)", "\\1",nearest_marker), 
                                trait_name = trait_to_name(trait.cfg.tb,model,trait), 
                                model_name = model_to_name(trait.cfg.tb,model,trait),
                                marker_variance=marker_variance/100,
                                model_variance=model_variance/100)
                                
} else {
    #Only reloading functions.  Load previous state.
    load(paste0(workflow,"/.RData.12_geneffects.",num_top_qtls))
}

generateTableSignifSymbols <- function(tbl) {
    return( tbl %>%
                mutate(trait_name=paste0(trait_name,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$"),
                    model_name=paste0(model_name,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) )
}

#arrange(trait_name,model_name,desc(marker.vaariance)) %>%
generateTableRemoveRepeats <- function(tbl) {
    return( tbl %>%
        mutate(trait_repeat=(trait_name == c("",trait_name[-length(trait_name)]))) %>%
        mutate(model_repeat=(trait_repeat & (model_name == c("",model_name[-length(model_name)])))) %>%
        mutate(trait_name = ifelse(trait_repeat, "",trait_name),
            model_name = ifelse(model_repeat, "",model_name)) )
}

generateReducedTable <- function(tbl, caption=NULL, tfont_size=10) {
    etbl1 <- tbl %>%
            generateTableRemoveRepeats() %>%
            mutate(model_variance = percent(model_variance),
                marker = ifelse(is.na(marker)," ",marker),
                position = round.digits(position,2),
                position_left = round.digits(interval_left,2),
                position_right = round.digits(interval_right,2),
                qtl_lod=round.digits(qtl_lod,2),
                marker_variance = percent(marker_variance))
    if( is_html_output() ) {
        #colorbars only render correctly under html format
        etbl1 <- etbl1 %>%
            mutate(model_variance = ifelse(model_repeat, "", color_bar("lightgreen")(model_variance)),
                marker_variance = color_bar("lightblue")(marker_variance))
    } else if (is_pdf_output() ) {
        etbl1 <- etbl1 %>%
            mutate(model_variance = knitr:::escape_latex(model_variance),
                   marker_variance = knitr:::escape_latex(marker_variance))
    }
    etbl2 <- etbl1 %>%
            dplyr::select(trait_name,chr,position,position_left,position_right,qtl_lod,marker_variance,model_variance) %>%
            rename("Trait"=trait_name, 
                    "LG"=chr, 
                    "Position (cM)"=position, 
                    "1.5-LOD Min (cM)"=position_left,
                    "1.5-LOD Max (cM)"=position_right,
                    "Variance Explained by QTL"=marker_variance,
                    "Model Variance[note]"=model_variance,
                    "pLOD[note]"=qtl_lod) %>% 
            mutate('Effect Size Boxplots[note]'="") %>%
            mutate('Effect Difference Plots[note]'="")
    if( is_word_output() ) {
        etbl4 <- etbl2 %>% knitr::kable(caption = caption, format = 'pipe', align = 'llrrrrrrrcc', escape = FALSE, longtable = TRUE, booktabs = TRUE)
    } else {
        etbl3 <- etbl2 %>% kable(caption=caption, align='lllrrrllcc', booktabs=TRUE, escape = FALSE, longtable = TRUE)
        etbl4 <- etbl3 %>%
                    kable_paper("striped", full_width=FALSE) %>%
                    column_spec(1, bold=TRUE, width = "1.0cm") %>%
                    column_spec(2, width = "0.5cm") %>%
                    column_spec(3, width = "0.75cm") %>%
                    column_spec(4, width = "0.75cm") %>%
                    column_spec(5, width = "0.75cm") %>%
                    column_spec(6, width = "0.75cm") %>%
                    column_spec(7, width = "1.0cm") %>%
                    column_spec(8, width = "1.0cm") %>%
                    column_spec(9, width = "6.0cm", image=spec_image(etbl1$plot_filename,width=702,height=175,res=300)) %>% #702x175 #To calculate width, take the width defined for column, convert to inches (/2.54), and multiply by res (300 dpi).  Then to find height, simply use aspect ratio to find height in pixels
                    column_spec(10, width = "3.0cm", image=spec_image(etbl1$plot_mpieffects_filename,width=350,height=175,res=300)) %>% #350x175
                    add_footnote(c("Variance of model with all significant QTLs fitted.",
                            "QTL Penalized LOD Score w/ significance codes:\n*** pvalue $≥$ 0 and pvalue<0.001\n**  pvalue $≥$ 0.001 and pvalue<0.01\n*   pvalue $≥$ 0.01 and pvalue<0.05\n.  pvalue $≥$ 0.05 and pvalue<0.01\nNS  Not Significant\n",
                            "Boxplots of nearest marker BLUPs grouped by genotypes.  Haplotypes A and B are from maternal parent P1, and haplotypes C and D are from paternal parent P2.",
                            "Effect differences for mean QTL effect size estimates for each progeny genotype.  A.-B. is the maternal effect, calculated as (AC+AD)-(BC+BD).  .C-.D is the paternal effect, calculated as (AC + BC) – (AD + BD).  Int is the interaction effect, calculated as (AC + BD)-(AD+BC) (Sewell et al., 2002)."),
                            escape=FALSE)
    }
    if( is_pdf_output() ) {
        etbl4 <- etbl4 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }
    return(etbl4)
}

effs.table.tbl <- effs.complete.tb %>% 
                    filter(model=='all-years' &
                            method=='stepwiseqtl' &
                            !is.na(nearest_marker)) %>%
                    group_by(trait,model) %>%
                    mutate(top_qtls = c(rep(TRUE,min(1,length(marker_variance))),rep(FALSE,length(marker_variance)-min(1,length(marker_variance))))) %>%
                    filter(top_qtls == TRUE) %>%
                    dplyr::select(-top_qtls) %>%
                    ungroup() %>%
                    arrange(trait,desc(marker_variance))

generateReducedTableNoPlots <- function(tbl, caption=NULL, tfont_size=10) {
    scipen_save = getOption('scipen') #Save current scipen value
    options(scipen=-3) #Set to negative value to render small doubles in scientific notation.
    #View(etbl1)
    #etbl1 <- effs.table.tbl %>%
    etbl1 <- tbl %>%
            generateTableRemoveRepeats() %>%
            mutate(model_variance = percent(round.digits(model_variance,3)),
                marker = ifelse(is.na(marker)," ",marker),
                position = round.digits(position,1),
                position_left = round.digits(interval_left,1),
                position_right = round.digits(interval_right,1),
                qtl_lod=round.digits(qtl_lod,1),
                marker_variance = percent(round.digits(marker_variance,3))) %>%
            mutate(AC_stat = paste0(signif.digits.char(effect_mean_AC,2), "±", signif.digits.char(effect_se_AC,2)),
                   AD_stat = paste0(signif.digits.char(effect_mean_AD,2), "±", signif.digits.char(effect_se_AD,2)),
                   BC_stat = paste0(signif.digits.char(effect_mean_BC,2), "±", signif.digits.char(effect_se_BC,2)),
                   BD_stat = paste0(signif.digits.char(effect_mean_BD,2), "±", signif.digits.char(effect_se_BD,2))) %>%
            mutate(AvB = signif.digits.char(AvB, 2),
                   CvD = signif.digits.char(CvD, 2),
                   Int = signif.digits.char(Int, 2))
    options(scipen=scipen_save)
    if( is_html_output() ) {
        #colorbars only render correctly under html format
        etbl1 <- etbl1 %>%
            mutate(model_variance = ifelse(model_repeat, "", color_bar("lightgreen")(model_variance)),
                marker_variance = color_bar("lightblue")(marker_variance))
    } else if (is_pdf_output() ) {
        etbl1 <- etbl1 %>%
            mutate(model_variance = knitr:::escape_latex(model_variance),
                   marker_variance = knitr:::escape_latex(marker_variance))
    }
    if( is_word_output() ) {
        etbl2 <- etbl1 %>%
                dplyr::select(trait_name,chr,position_left,position,position_right,model_variance,marker_variance,AvB,CvD,Int,AC_stat,AD_stat,BC_stat,BD_stat) %>%
                rename("Trait"=trait_name, 
                       "LG"=chr, 
                       "1.5-LOD min"=position_left,
                       "Position"=position, 
                       "1.5-LOD max"=position_right,
                       "R~m~^2^"=model_variance,
                       "R~q~^2^"=marker_variance,
                       "A.-B."=AvB,
                       ".C-.D"=CvD,
                       "Int"=Int,
                       "eff~AC~"=AC_stat,
                       "eff~AD~"=AD_stat,
                       "eff~BC~"=BC_stat,
                       "eff~BD~"=BD_stat)
        etbl4 <- flextable(etbl2) %>%
            ftExtra::colformat_md(part="header") %>%
            align(align="center", part="header") %>%
            align(align="right", part="body") %>%
            set_caption(caption) %>%
            flextable::footnote(i = 1, 
                    j = c("1.5-LOD min","Position","1.5-LOD max","R~m~^2^","R~q~^2^","A.-B.",".C-.D","Int","eff~AC~","eff~AD~","eff~BC~","eff~BD~"),
                    value=as_paragraph(c("1.5 below peak LOD left position (cM).",
                          "QTL position of peak LOD (cM).",
                          "1.5 below peak LOD right position (cM).",
                          "Percent of phenotypic variance explained by all signficant additive effect QTL fit by model.",
                          "Percent of additive genetic variance explained by QTL.",
                          "Maternal effect size: (AC+AD)-(BC+BD)",
                          "Paternal effect size: (AC+BC)–(AD+BD)",
                          "Interaction effect size: (AC+BD)-(AD+BC)",
                          "Effect size of AC genotype.",
                          "Effect size of AD genotype.",
                          "Effect size of BC genotype.",
                          "Effect size of BD genotype.")),
                    ref_symbols = c("a","b","c","d","e","f","g","h","i","j","k","l"),
                    part = "header") %>%
            theme_booktabs(bold_header=TRUE) %>%
            set_table_properties(layout="autofit")
    } else {
        etbl2 <- etbl1 %>%
                dplyr::select(trait_name,chr,position,position_left,position_right,qtl_lod,marker_variance,model_variance,AvB,CvD,Int,AC_stat,AD_stat,BC_stat,BD_stat) %>%
                rename("Trait"=trait_name, 
                        "LG"=chr, 
                        "Position (cM)"=position, 
                        "1.5-LOD Min (cM)"=position_left,
                        "1.5-LOD Max (cM)"=position_right,
                        "Variance Explained by QTL"=marker_variance,
                        "Model Variance[note]"=model_variance,
                        "pLOD[note]"=qtl_lod,
                        "A.-B.[note]"=AvB,
                        ".C-.D[note]"=CvD,
                        "Int[note]"=Int,
                        "$eff_{AC}$[note]"=AC_stat,
                        "$eff_{AD}$[note]"=AD_stat,
                        "$eff_{BC}$[note]"=BC_stat,
                        "$eff_{BD}$[note]"=BD_stat)
        etbl3 <- etbl2 %>% kable(caption=caption, align='llrrrrrrrrrrrr', booktabs=TRUE, escape = FALSE, longtable = TRUE)
        etbl4 <- etbl3 %>%
                    kable_paper("striped", full_width=FALSE) %>%
                    column_spec(1, bold=TRUE, width = "1.0cm") %>%
                    column_spec(2, width = "0.5cm") %>%
                    column_spec(3, width = "0.75cm") %>%
                    column_spec(4, width = "0.75cm") %>%
                    column_spec(5, width = "0.75cm") %>%
                    column_spec(6, width = "0.75cm") %>%
                    column_spec(7, width = "1.0cm") %>%
                    column_spec(8, width = "1.0cm") %>%
                    column_spec(9, width = "0.75cm") %>%
                    column_spec(10, width = "0.75cm") %>%
                    column_spec(11, width = "0.75cm") %>%
                    column_spec(12, width = "1.0cm") %>%
                    column_spec(13, width = "1.0cm") %>%
                    column_spec(14, width = "1.0cm") %>%
                    column_spec(15, width = "1.0cm") %>%
                    add_footnote(c("Variance of model with all significant QTLs fitted.",
                            "QTL Penalized LOD Score w/ significance codes:\n*** pvalue $≥$ 0 and pvalue<0.001\n**  pvalue $≥$ 0.001 and pvalue<0.01\n*   pvalue $≥$ 0.01 and pvalue<0.05\n.  pvalue $≥$ 0.05 and pvalue<0.01\nNS  Not Significant\n",
                            "Maternal effect size - (AC+AD)-(BC+BD)",
                            "Paternal effect size - (AC+BC)–(AD+BD)",
                            "Interaction effect size - (AC+BD)-(AD+BC)",
                            "AC effect size",
                            "AD effect size",
                            "BC effect size",
                            "BD effect size"),
                            escape=FALSE)
    }
    if( is_pdf_output() ) {
        etbl4 <- etbl4 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }
    return(etbl4)
}

generateTable <- function(tbl, meths=c("scanone","stepwiseqtl"), caption=NULL, tfont_size=10) {
    etbl1 <- tbl %>%
            filter(method %in% meths) %>%
            generateTableSignifSymbols() %>%
            generateTableRemoveRepeats() %>%
            mutate(model_variance = percent(model_variance),
                marker = ifelse(is.na(marker)," ",marker),
                position = round.digits(position,2),
                position_left = round.digits(interval_left,2),
                position_right = round.digits(interval_right,2),
                qtl_lod=round.digits(qtl_lod,2),
                marker_variance = percent(marker_variance))
    if( is_html_output() ) {
        #colorbars only render correctly under html format
        etbl1 <- etbl1 %>%
            mutate(model_variance = ifelse(model_repeat, "", color_bar("lightgreen")(model_variance)),
                marker_variance = color_bar("lightblue")(marker_variance))
    } else if (is_pdf_output() ) {
        etbl1 <- etbl1 %>%
            mutate(model_variance = knitr:::escape_latex(model_variance),
                   marker_variance = knitr:::escape_latex(marker_variance))
    }

    etbl2 <- etbl1 %>%
            dplyr::select(trait_name,model_name,chr,position,position_left,position_right,qtl_lod,marker_variance,model_variance) %>%
            rename("Trait[note]"=trait_name, 
                    "Model[note]"=model_name, 
                    "LG"=chr, 
                    "Position (cM)"=position, 
                    "1.5-LOD Min (cM)"=position_left,
                    "1.5-LOD Max (cM)"=position_right,
                    "Variance Explained by QTL"=marker_variance,
                    "Model Variance[note]"=model_variance, 
                    "pLOD[note]"=qtl_lod) %>% 
            mutate('Effect Size Boxplots[note]'="") %>%
            mutate('Effect Difference Plots[note]'="")
        etbl3 <- etbl2 %>% kable(caption = caption, align = 'llrrrrrrrcc', escape = FALSE, longtable = TRUE, booktabs = TRUE)
        etbl4 <- etbl3 %>%
                    row_spec(row=(which(etbl1$trait_name != "")-1)[-1], hline_after = TRUE) %>%
                    kable_paper("striped", full_width=FALSE) %>%
                    column_spec(1, bold=TRUE, width = "1.0cm") %>%
                    column_spec(2, bold=TRUE, width = "1.0cm") %>%
                    column_spec(3, width = "0.5cm") %>%
                    column_spec(4, width = "0.75cm") %>%
                    column_spec(5, width = "0.75cm") %>%
                    column_spec(6, width = "0.75cm") %>%
                    column_spec(7, width = "0.75cm") %>%
                    column_spec(8, width = "1.0cm") %>%
                    column_spec(9, width = "1.0cm") %>%
                    column_spec(10, width = "6.0cm", image=spec_image(etbl1$plot_filename,width=702,height=175,res=300)) %>% #702x175 #To calculate width, take the width defined for column, convert to inches (/2.54), and multiply by res (300 dpi).  Then to find height, simply use aspect ratio to find height in pixels
                    column_spec(11, width = "3.0cm", image=spec_image(etbl1$plot_mpieffects_filename,width=350,height=175,res=300)) %>% #350x175
                    add_footnote(c(paste0("All QTLs in table derived from running R/qtl package function ",meths,"().\n","Significance codes for model genotype$*$year effects appended to trait:\n*** pvalue $≥$ 0 and pvalue<0.001\n**  pvalue $≥$ 0.001 and pvalue<0.01\n*   pvalue $≥$ 0.01 and pvalue<0.05\n.  pvalue $≥$ 0.05 and pvalue<0.01\nNS  Not Significant\n"), 
                            "Significance codes for model genotype effects appended to model",
                            "Variance of model with all significant QTLs fitted.",
                            "Penalized LOD Score w/ significance codes for QTL appended",
                            "Boxplots of nearest marker BLUPs grouped by genotypes.  Haplotypes A and B are from maternal parent P1, and haplotypes C and D are from paternal parent P2.",
                            "Effect differences for mean QTL effect size estimates for each progeny genotype.  A.-B. is the maternal effect, calculated as (AC+AD)-(BC+BD).  .C-.D is the paternal effect, calculated as (AC + BC) – (AD + BD).  Int is the interaction effect, calculated as (AC + BD)-(AD+BC) (Sewell et al., 2002)."),
                            escape=FALSE)
    if( is_pdf_output() ) {
        etbl4 <- etbl4 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }
    return(etbl4)
}

save.image(paste0(workflow,"/.RData.12_geneffects.",num_top_qtls))

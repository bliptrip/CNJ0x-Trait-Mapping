#!/usr/bin/env RScript

#NOTE: This particular script generates informative data plots for the BLUPs: Things like pairwise correlation plots, boxplots, etc.
# Think of combining BLUP data into a table format similar to the input phenotype table, and then can display data in a way similar
# to how this was done for the original phenotypes.
#
# GGally is our friend here!

# loading libraries
library(cowplot)
library(equatags)
library(flextable)
library(formattable)
library(GGally)
library(ggplot2)
library(ggthemes)
library(jsonlite)
library(kableExtra)
library(knitr)
library(RColorBrewer)
library(tidyverse)
options(knitr.kable.NA='') #So we don't display NA's in table renderings.

workflow <- get0("workflow", ifnotfound="../../Workflows/9")
P1_Name <- get0("P1_Name", ifnotfound="Mullica_Queen")
P2_Name <- get0("P2_Name", ifnotfound="Crimson_Queen")
POP_Name <- get0("POP_Name", ifnotfound="CNJ02")
table.width <- get0("table.width", ifnotfound=0.85)
reload_table_functions <- get0("reload_table_functions", ifnotfound=FALSE)

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

reload_table_functions=as.logical(reload_table_functions)

if( !reload_table_functions ) {
    source('./usefulFunctions.R')
    source(paste0(workflow,"/configs/model.cfg"))

    cat(paste0("Population name: ",POP_Name))

    randomEffects2String <- function(randoms) {
        randoms <- gsub("vs(id, Gu = A)", "Z_{g}g", randoms, fixed=TRUE)
        randoms <- gsub("id:year", "Z_{ge}ge", randoms, fixed=TRUE)
        randoms <- gsub("vs(spl2D(row, column))", "Z_{s}s", randoms, fixed=TRUE)
        randoms <- gsub("vs(columnf)", "Z_{c}c", randoms, fixed=TRUE)
        randoms <- gsub("vs(rowf)", "Z_{r}r", randoms, fixed=TRUE)
        return(randoms)
    }

    collateBLUPs <- function(trait.cfg, trait.path, args.l) {
        models.tb.p         <- args.l[[2]]
        trait_name          <- trait.cfg$trait
        model_id         <- trait.cfg$model
        cat(paste0("Processing trait: ", trait_name, ", Model: ", model_id,"\n"))
        model <- readRDS(paste0(trait.path,"/mmer.rds"))
        varcomp      <- model$sigma
        g.idx     <- which(grepl("^u:id",names(varcomp)))
        ge.idx    <- which(grepl("^id:year",names(varcomp)))
        e.idx      <- which(grepl("^units",names(varcomp)))
        vg        <- varcomp[[g.idx]][trait_name,trait_name]
        vge       <- ifelse(length(ge.idx) > 0,varcomp[[ge.idx]][trait_name,trait_name],NA)
        ve        <- varcomp[[e.idx]][trait_name,trait_name]
        if( model_id == "all-years" ) {
            num_years  <- length(model$Beta$Estimate) #Number of years based on model estimate
            if(is.na(vge)) {
                h2 <- vg/(vg+(ve/num_years))
            } else {
                h2 <- vg/(vg+(vge/num_years)+(ve/num_years))
            }
        } else {
            h2 <- vg/(vg+ve)
        }
        u                    <- model$U[["u:id"]]
        blup_adjust  <- sum(model$Beta$Estimate)
        blups                <- u[[trait_name]] + blup_adjust
        ids                 <- names(blups)
        anova.file            <- paste0(trait.path,'/anova.csv')
        GLRpvalue            <- NA
        GZRpvalue            <- NA
        GxYLRpvalue            <- NA
        GxYZRpvalue            <- NA
        model_formula          <- randomEffects2String(as.character(model$call$random)[[2]])
        if( file.exists(anova.file) ) {
            anova.df=read.csv(file=anova.file, header=TRUE, row.names=1)
            GLRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% dplyr::select(PrChisq))
            GZRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% dplyr::select(PrNorm))
            if( model_id == "all-years" ) {
                GxYLRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% dplyr::select(PrChisq))
                GxYZRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% dplyr::select(PrNorm))
            }
        }
        if( model_id == "all-years" ) {
            pheno.pre <- model$dataOriginal[,c("id",trait_name)] %>% group_by(id) %>% dplyr::summarize(mean=mean(.data[[trait_name]],na.rm=TRUE))
            pheno <- pheno.pre$mean
            names(pheno) <- pheno.pre$id
        } else {
            pheno <- model$dataOriginal[[trait_name]]
            names(pheno) <- model$dataOriginal$id
        }
        pheno <- pheno[ids] #Make sure that the raw phenotype data is in the same order as BLUPs
        pheno.sd   <- sd(pheno, na.rm=TRUE)
        pheno.u    <- mean(pheno, na.rm=TRUE)
        phenoZ     <- (pheno - pheno.u)/pheno.sd
        blupsZ     <- (blups - pheno.u)/pheno.sd
        models.tb.p$value <- models.tb.p$value %>% 
                                add_row(id=factor(ids), model=factor(trait.cfg$model), model_label=factor(trait.cfg$model_label), trait=factor(trait.cfg$trait), model_formula=model_formula, value=pheno, valueZ=phenoZ, type=factor("Raw"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2) %>%
                                add_row(id=factor(ids), model=factor(trait.cfg$model), model_label=factor(trait.cfg$model_label), trait=factor(trait.cfg$trait), model_formula=model_formula, value=blups, valueZ=blupsZ, type=factor("BLUPs"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2) 
    }

    pheno.means.df <-read.csv(file=pheno_dpath2fpath(pheno_file))
    pheno.means.df$year <- as.factor(pheno.means.df$year) #Needed for modeling column effects
    traits.df   <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=T)
#Just build a 'long' version of collated BLUPs as tibble table, and then use tidyverse 'pivot_wider()' to flatten
    model.collated.long.tb       <- tibble(id=factor(), model=factor(), model_label=factor(), trait=factor(), model_formula=character(), value=numeric(), valueZ=numeric(), type=factor(), label=factor(), label_short=factor(), GLRpvalue=numeric(), GZRpvalue=numeric(), GxYLRpvalue=numeric(), GxYZRpvalue=numeric(), vg=numeric(), vge=numeric(), ve=numeric(), h2=numeric())
    model.collated.long.tb.p   <- newPointer(model.collated.long.tb)
    selectedmodels.df          <- read.csv(file=paste0(workflow,"/traits/selectedModels.csv"))
    loopThruTraits(workflow, collateBLUPs, loopArgs=list(selectedmodels.df,model.collated.long.tb.p))
    model.collated.long.tb     <- model.collated.long.tb.p$value %>%
                                    group_by(trait) %>%
                                    mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>% #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                                    mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                                    ungroup()
#Save long version, as it can be manipulated in many ways to produce meaningful plots/graphs
    saveRDS(model.collated.long.tb, file=paste0(workflow,'/traits/blups_collated.long.rds'), compress=TRUE)
    model.collated.raw.tb <- model.collated.long.tb %>% 
                                filter(type == "Raw") %>% 
                                dplyr::select(c('id','model','trait','value','valueZ')) %>%
                                rename(raw=value) %>%
                                rename(rawZ=valueZ)
    model.collated.blup.tb <- model.collated.long.tb %>% 
                                filter(type == "BLUPs") %>%
                                dplyr::select(-type) %>%
                                rename(blup=value) %>%
                                rename(blupZ=valueZ)
    model.collated.wide.tb <- model.collated.blup.tb %>% 
                                left_join(model.collated.raw.tb, by=c('id','model','trait') ) %>%
                                rename(genotype=id) %>%
                                mutate(id=seq(1,n()))

    saveRDS(model.collated.wide.tb, file=paste0(workflow,'/traits/blups_collated.wide.rds'), compress=TRUE)
    write_json(model.collated.wide.tb, paste0(workflow, '/traits/blups_collated.wide.json'), auto_unbox=T, pretty=T) #Also write to JSON file for easy import by other programming languages/tools
    write_csv(model.collated.wide.tb, file=paste0(workflow,'/traits/blups_collated.wide.csv'))

#Spread out
#model.collated.tb <- spread(model.collated.long.tb, trait, blup)
#write.csv(model.collated.tb, file=paste0(workflow,'/traits/blups_collated.csv'), row.names=F)

#Generate one plot per unique model

#Generate plot across only models with valid 'year' -- not 'all-years'
#model.collated.years.tb <- model.collated.tb %>% filter(model != 'all-years')
#full_traits <- c()
#for (trait in unique(model.collated.long.tb$trait)) {
#    if (!any(is.na(model.collated.years.tb[,trait]))) {
#        full_traits <- c(full_traits, trait)
#    }
#}
#g <- ggscatmat(model.collated.years.tb[,c('id','model',full_traits)], color="model", alpha=0.8) +
#         theme_minimal() +
#         theme(axis.title = element_blank(),
#               axis.text.x  = element_text(size=14,angle=60),
#               axis.text.y  = element_text(size=14),
#               strip.text = element_text(face="bold",size=18),
#               legend.text=element_text(face="bold",size=18),
#               legend.title=element_text(face="bold",size=20))
#png(filename=paste0(workflow,'/traits/plots/blups_collated.traits.years.png'), width=1280, height=960, bg="white")
#g
#dev.off()
#
##Generate a plot across all models.  NOTE: Need to remove any traits with any NA's as these will result in ggscatmat returning an error (unbalanced)
#model.collated.allmodels.tb <- model.collated.tb
#full_traits <- c()
#for (trait in unique(model.collated.long.tb$trait)) {
#    if (!any(is.na(model.collated.allmodels.tb[,trait]))) {
#        full_traits <- c(full_traits, trait)
#    }
#}
#g <- ggscatmat(model.collated.allmodels.tb[,c('id','model',full_traits)], color="model", alpha=0.8) +
#         theme_minimal() +
#         theme(axis.title = element_blank(),
#               axis.text.x  = element_text(size=8,angle=60),
#               axis.text.y  = element_text(size=8),
#               strip.text = element_text(face="bold",size=10),
#               legend.text=element_text(face="bold",size=8),
#               legend.title=element_text(face="bold",size=10))
#png(filename=paste0(workflow,'/traits/plots/blups_collated.traits.allmodels.png'), width=1280, height=960, bg="white")
#g
#dev.off()

    model.collated.long.tb <- model.collated.long.tb %>% 
                                drop_na(value) %>%
                                mutate(model_trait=as.factor(paste0(model,' ',trait)), model_type=as.factor(paste0(model,' ',type)))
    all_model_types <- model.collated.long.tb %>% dplyr::select(model_trait) %>% unique()

#Inject empty entries

    model.collated.long.stats.tb <- model.collated.long.tb %>%
                                        filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                        group_by(model,trait,type,model_type) %>%
                                        arrange(desc(value)) %>%
                                        dplyr::summarize(min=min(value), 
                                                        minZ=min(valueZ), 
                                                        min_geno = gsub(paste0(POP_Name,"_([0-9]+_[0-9]+)"),"g\\1",id[which.min(value)]),
                                                        max=max(value), 
                                                        maxZ=max(valueZ), 
                                                        max_geno = gsub(paste0(POP_Name,"_([0-9]+_[0-9]+)"),"g\\1",id[which.max(value)]),
                                                        top_genos = paste0(gsub(paste0(POP_Name,"_([0-9]+_[0-9]+)"),"g\\1",head(id,5)),collapse=","),
                                                        bottom_genos = paste0(gsub(paste0(POP_Name,"_([0-9]+_[0-9]+)"),"g\\1",tail(id,5)),collapse=","),
                                                .groups="keep") %>%
                                        filter(min != max) %>%
                                        ungroup() %>%
                                        group_nest(trait, .key="stats")

    model.collated.long.ylims.tb <- model.collated.long.tb %>%
                                        filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                        group_by(trait,type) %>%
                                        filter(min(value)!= max(value)) %>%
                                        dplyr::summarize(min=min(value, na.rm=TRUE), minZ=min(valueZ, na.rm=TRUE),
                                                max=max(value, na.rm=TRUE), maxZ=max(valueZ, na.rm=TRUE),
                                                .groups="keep") %>%
                                        dplyr::mutate(range=max-min, 
                                                    rangeZ = maxZ-minZ) %>%
                                        dplyr::mutate(min = min - range/8,
                                                    max = max + range/8,
                                                    minZ = minZ - rangeZ/8,
                                                    maxZ = maxZ + rangeZ/8) %>%
                                        gather("limit","value",-c("trait","type"))



    model.collated.long.ylims.tb$model <- model.collated.long.tb$model[1]
    model.collated.long.ylims.tb$model_type <- paste0(model.collated.long.tb$model[1], " ", model.collated.long.ylims.tb$type)#Just a 'placeholder' to allow ggplot to work with using this dataset
    model.collated.long.ylims.tb <- model.collated.long.ylims.tb %>% 
                                        ungroup() %>% 
                                        group_nest(trait, .key="ylims")

#Display boxplots of Raw Values + BLUPs without Z-score scaling
    gs <- model.collated.long.tb %>% 
            group_by(model,trait,type) %>%
            filter(min(value) != max(value)) %>%
            ungroup() %>%
            group_nest(trait) %>%
            left_join(model.collated.long.stats.tb, by="trait") %>%
            left_join(model.collated.long.ylims.tb, by="trait") %>%
            group_by(trait) %>%
            do(
            plot = ggplot(.$data[[1]], aes(x=factor(model), y=value, group=factor(model_type), fill=factor(type))) +
                        geom_boxplot(alpha=0.4, position=position_dodge(1)) +
                        geom_jitter(alpha=0.1, position=position_jitterdodge(jitter.width=0.5, dodge.width=1)) +
                        geom_text(data=.$stats[[1]], position=position_dodge(1), aes(y=min, label=min_geno), size=8, fontface='bold', vjust=1.5, hjust=1.5, angle=45) +
                        geom_text(data=.$stats[[1]], position=position_dodge(1), aes(y=max, label=max_geno), size=8, fontface='bold', vjust=-1.5, hjust=0, angle=45) +
                        geom_blank(data=.$ylims[[1]] %>% filter((limit == "min") | (limit == "max"))) +
                        ylab("") +
                        xlab("Year") + 
                        guides(fill='none') +
                        ggtitle(label=.$data[[1]]$label_short) +
                        theme(  axis.text.x = element_text(face="bold", size=32, angle = 60, hjust = 1),
                                axis.text.y = element_text(face="bold", size=32),
                                legend.title = element_text(fac="bold", size=32),
                                legend.text = element_text(fac="bold", size=24),
                                legend.key.size = ggplot2::unit(5,"cm"),
                                axis.title  = element_text(face="bold",size=48),
                                plot.title   = element_text(face="bold",size=52, hjust=0.5),
                                plot.subtitle = element_text(size=48, hjust = 0.5),
                                plot.margin = ggplot2::unit(c(1,1,1,1),"cm")))

#Yield-related traits
##Upright Yield
    p1 <- (gs %>% filter(trait == "berry_length"))$plot[[1]] + ylab("Value")
    p2 <- (gs %>% filter(trait == "berry_weight"))$plot[[1]]
    p3 <- (gs %>% filter(trait == "berry_width"))$plot[[1]]
    p4 <- (gs %>% filter(trait == "total_berry_weight"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
    p5 <- (gs %>% filter(trait == "num_peds"))$plot[[1]] + ylab("Value")
    p6 <- (gs %>% filter(trait == "num_seeds"))$plot[[1]]
    p7 <- (gs %>% filter(trait == "UMFM"))$plot[[1]]
    p8 <- (gs %>% filter(trait == "ULvW"))$plot[[1]]
    pg <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4,rel_widths=c(1,1,1,1.25))
    png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.upright_yield.png'), width=3840, height=2560, bg="white")
    pg
    dev.off()
##Upright Yield-Secondary Traits
    p1 <- (gs %>% filter(trait == "upright_length"))$plot[[1]] + ylab("Value")
    p2 <- (gs %>% filter(trait == "secondary_growth"))$plot[[1]]
    p3 <- (gs %>% filter(trait == "dry_wt_leaves"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
    pg <- plot_grid(p1,p2,p3,nrow=1,ncol=3,rel_widths=c(1,1,1.25))
    png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.upright_yield2.png'), width=2880, height=1280, bg="white")
    pg
    dev.off()

#Biennial-bearing traits
    p1 <- (gs %>% filter(trait == "BBI_UTBM"))$plot[[1]] + ylab("Value")
    p2 <- (gs %>% filter(trait == "BBI_TY"))$plot[[1]]
    p3 <- (gs %>% filter(trait == "BBI_SFY"))$plot[[1]]
    p4 <- (gs %>% filter(trait == "rebud"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
    pg <- plot_grid(p1,p2,p3,p4,nrow=1,ncol=4,rel_widths=c(1,1,1,1.25))
    png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.bbi.png'), width=3840, height=1280, bg="white")
    pg
    dev.off()

#Whole Plot Traits as collected by Vorsa's Group
    p1 <- (gs %>% filter(trait == "TY"))$plot[[1]] + ylab("Value")
    p2 <- (gs %>% filter(trait == "SFY"))$plot[[1]]
    p3 <- (gs %>% filter(trait == "MFM"))$plot[[1]]
    p4 <- (gs %>% filter(trait == "PFR"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
    p5 <- (gs %>% filter(trait == "Tacy"))$plot[[1]] + ylab("Value")
    p6 <- (gs %>% filter(trait == "Brix"))$plot[[1]]
    p7 <- (gs %>% filter(trait == "TA"))$plot[[1]]
    p8 <- (gs %>% filter(trait == "PAC"))$plot[[1]]
    pg <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4,rel_widths=c(1,1,1,1.25))
    png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.plot_traits.png'), width=3840, height=2560, bg="white")
    pg
    dev.off()



    g <- model.collated.long.tb %>% 
            group_by(model,trait,type) %>%
            filter(min(value) != max(value)) %>%
            ungroup() %>%
            arrange(model,trait,type) %>%
            ggplot(aes(x=factor(model), y=value, fill=factor(trait))) +
            geom_boxplot(alpha=0.4, position=position_dodge(1)) +
            geom_jitter(alpha=0.1, position=position_jitterdodge(jitter.width=0.5, dodge.width=1)) +
            geom_text(data=model.collated.long.stats.tb %>% unnest(stats), position=position_dodge(1), aes(y=min, label=min_geno), size=8, fontface='bold', vjust=1.5, hjust=1.5, angle=45) +
            geom_text(data=model.collated.long.stats.tb %>% unnest(stats), position=position_dodge(1), aes(y=max, label=max_geno), size=8, fontface='bold', vjust=-1.5, hjust=0, angle=45) +
            geom_blank(data=model.collated.long.ylims.tb %>% unnest(ylims) %>% filter((limit == "min") | (limit == "max"))) +
            facet_grid(rows=vars(type), cols=vars(trait), scales="free") +
            guides(fill = 'none') +
            xlab("Year") + 
            ylab("Trait Values") +
            ggtitle(label=paste0(POP_Name," Population Raw Values and BLUPs"), subtitle="Yield-Related Traits") +
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
            group_by(model,trait,type) %>%
            filter(min(value) != max(value)) %>%
            ungroup() %>%
            arrange(model,trait,type) %>%
            ggplot(aes(x=factor(model), y=valueZ, fill=factor(trait))) +
            geom_blank(data=model.collated.long.ylims.tb %>% unnest(ylims) %>% filter((limit == "minZ") | (limit == "maxZ")), aes(x=factor(model), y=value)) +
            geom_boxplot(alpha=0.4) +
            geom_jitter(width=0.1, alpha=0.4) +
            geom_text(data=model.collated.long.stats.tb %>% unnest(stats), aes(y=minZ, label=min_geno), size=8, fontface='bold', vjust=1.5, hjust=1.5, angle=45) +
            geom_text(data=model.collated.long.stats.tb %>% unnest(stats), aes(y=maxZ, label=max_geno), size=8, fontface='bold', vjust=-1.5, hjust=0, angle=45) +
            facet_grid(rows=vars(type), cols=vars(trait), scales="free_y") +
            guides(fill = 'none') +
            xlab("Year") + 
            ylab("Trait Values") +
            ggtitle(label=paste0(POP_Name," Population Raw Values and BLUPs Z-Score Scaled"), subtitle="Yield-Related Traits") +
            theme( axis.text.x = element_text(face="bold", size=32, angle = 60, hjust = 1),
                axis.text.y = element_text(face="bold", size=32),
                strip.text  = element_text(face="bold", size=36, angle = 30),
                axis.title  = element_text(face="bold",size=48),
                plot.title   = element_text(face="bold",size=52, hjust=0.5),
                plot.subtitle = element_text(size=48, hjust = 0.5))
    png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.allmodels.zscore.png'), width=2560, height=1920, bg="white")
    g
    dev.off()

#Generate a table of BLUP summaries by model and trait, with min value, min genotype, max value, max genotype, P1 blup value, P2 blup value

    p1.blup_summary.tb <- model.collated.long.tb %>%
                            filter(id == P1_Name) %>%
                            filter(model == "all-years")
                            
    p1.blup_summary.tb <- model.collated.long.tb %>%
                            filter(id == P1_Name) %>%
                            group_by(model,trait,type) %>%
                            dplyr::summarize(P1 = mean(value), P1z = mean(valueZ))

    p2.blup_summary.tb <- model.collated.long.tb %>%
                            filter(id == P2_Name) %>%
                            group_by(model,trait,type) %>%
                            dplyr::summarize(P2 = mean(value), P2z = mean(valueZ))

    model.collated.summary.tb <- model.collated.long.tb %>%
                                        filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                        group_by(model,model_label,trait,label,label_short,model_formula,type) %>%
                                        arrange(desc(value)) %>%
                                        dplyr::summarize(min=min(value), 
                                                        mean=mean(value), 
                                                        sd=sd(value), 
                                                        range=abs(min(value)-max(value)),
                                                        minZ=min(valueZ), 
                                                        min_geno = gsub(paste0(POP_Name,"_([0-9]+)_([0-9]+))"),"$g\\1_{\\2}$",id[which.min(value)]),
                                                        bottom_genos = paste0(gsub(paste0(POP_Name,"_([0-9]+)_([0-9]+)"),"$g\\1_{\\2}$",tail(id,5)),collapse=", "),
                                                        max=max(value), 
                                                        maxZ=max(valueZ), 
                                                        max_geno = gsub(paste0(POP_Name,"_([0-9]+)_([0-9]+)"),"$g\\1_{\\2}$",id[which.max(value)]),
                                                        top_genos = paste0(gsub(paste0(POP_Name,"_([0-9]+)_([0-9]+)"),"$g\\1_{\\2}$",head(id,5)),collapse=", "),
                                                        GLRpvalue=mean(GLRpvalue),
                                                        GZRpvalue=mean(GZRpvalue),
                                                        GxYLRpvalue=mean(GxYLRpvalue),
                                                        GxYZRpvalue=mean(GxYZRpvalue),
                                                        vg=mean(vg), 
                                                        vge=mean(vge), 
                                                        ve=mean(ve), 
                                                        h2=mean(h2)) %>%
                                        filter(min != max) %>%
                                        left_join(p1.blup_summary.tb, by=c("model","trait","type")) %>%
                                        left_join(p2.blup_summary.tb, by=c("model","trait","type")) %>%
                                        mutate('Parent Range' = signif.digits(abs(P2-P1),2)) %>%
                                        ungroup() #This is necessary for trait_repeat and model_repeat calculations to work below

    model.collated.summary.sub.tb <- model.collated.summary.tb %>% select(!c(type,min,mean,sd,range,minZ,min_geno,max,maxZ,max_geno,bottom_genos,top_genos,P1,P2,'Parent Range'))
    model.collated.summary.wide.tb <- model.collated.summary.tb %>% pivot_wider(id_cols=c(model,model_label,trait,label,label_short),names_from=c(type),values_from=c(min,mean,sd,range,minZ,min_geno,max,maxZ,max_geno,top_genos,bottom_genos,P1,P2,'Parent Range')) %>% inner_join(model.collated.summary.sub.tb, by=c("model","model_label","trait","label","label_short")) %>% distinct(across(c(model,trait)),.keep_all=TRUE) 
    write_csv(model.collated.summary.wide.tb,file=paste0(workflow,"/traits/blups_collated.summary.wide.csv"))

    model.collated.blup_summary.tb <- model.collated.summary.tb %>%
            filter(type == "BLUPs")

    model.collated.raw_summary.tb <- model.collated.summary.tb %>%
            filter(type == "Raw")
} else {
    load(paste0(workflow,"/.RData.01_genBLUP"))
}

generateTableSignifSymbols <- function(values.collated.tb) {
    return( values.collated.tb %>%
                mutate(label=paste0(label,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$")) %>%
                mutate(model_label=paste0(model_label,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) )
}

generateTableRemoveRepeats <- function(values.collated.tb) {
    return( values.collated.tb %>%
        mutate( label = ifelse(trait_repeat, " ", label),
                model_label = ifelse(model_repeat, " ", model_label)) )
}


generateFullTableHelper <- function(values.collated.tb) {
    tbl1 <- values.collated.tb %>%
        arrange(trait,model) %>%
        mutate(label = as.character(label),
               model_label = as.character(model_label)) %>% #Necessary for trait_repeat and model_repeat to be calculated correctly
        mutate(trait_repeat=(label == c("",label[-length(label)])),
                model_repeat=(model_label == c("",model_label[-length(model_label)]))) %>%
        mutate(min_Raw=signif.digits(min_Raw,2)) %>%
        mutate(min_BLUPs=signif.digits(min_BLUPs,2)) %>%
        mutate(max_Raw=signif.digits(max_Raw,2)) %>%
        mutate(max_BLUPs=signif.digits(max_BLUPs,2)) %>%
        mutate(mean_Raw=ifelse(!is.na(mean_Raw),paste0(signif.digits(mean_Raw,2),"±",signif.digits(sd_Raw,2))," ")) %>% 
        mutate(mean_BLUPs=ifelse(!is.na(mean_BLUPs),paste0(signif.digits(mean_BLUPs,2),"±",signif.digits(sd_BLUPs,2))," ")) %>% 
        mutate(range_Raw=signif.digits(range_Raw,2)) %>%
        mutate(range_BLUPs=signif.digits(range_BLUPs,2)) %>%
        mutate(P1_Raw=signif.digits(P1_Raw,2)) %>%
        mutate(P1_BLUPs=signif.digits(P1_BLUPs,2)) %>%
        mutate(P2_Raw=signif.digits(P2_Raw,2)) %>%
        mutate(P2_BLUPs=signif.digits(P2_BLUPs,2)) %>%
        mutate(vg=signif.digits(vg,2)) %>%
        mutate(vge=signif.digits(vge,2)) %>%
        mutate(ve=signif.digits(ve,2)) %>%
        mutate(h2=signif.digits(h2,2))
    return(tbl1)
}

generateFullTable <- function(values.collated.tb, caption=NULL, tfont_size=10) {
    tbl1 <- generateFullTableHelper(values.collated.tb) %>%
                generateTableSignifSymbols() %>%
                generateTableRemoveRepeats() %>%
                dplyr::select(label,model_label,vg,vge,ve,h2,min_Raw,min_BLUPs,mean_Raw,mean_BLUPs,max_Raw,max_BLUPs,range_Raw,range_BLUPs,min_geno_Raw,min_geno_BLUPs,max_geno_Raw,max_geno_BLUPs,P1_Raw,P1_BLUPs,P2_Raw,P2_BLUPs,'Parent Range_Raw','Parent Range_BLUPs')
    tbl2 <- tbl1 %>%
        rename("Trait[note]"=label,"Model[note]"=model_label,'$σ_{g}^{2}$[note]'=vg,'$σ_{g ε}^{2}$[note]'=vge,'$σ_{ε}^{2}$[note]'=ve,'$h^{2}$[note]'=h2,'$Min_r$'=min_Raw,'$Min_b$'=min_BLUPs,'$Max_r$'=max_Raw,'$Max_b$'=max_BLUPs,"$μ_r$±$SE$"=mean_Raw,"$μ_b$±$SE$"=mean_BLUPs,'$R_r$[note]'=range_Raw,'$R_b$[note]'=range_BLUPs,'$MinG_r$[note]'=min_geno_Raw,'$MinG_b$'=min_geno_BLUPs,'$MaxG_r$[note]'=max_geno_Raw,'$MaxG_b$'=max_geno_BLUPs,'$P_{1r}$[note]'=P1_Raw,'$P_{1b}$'=P1_BLUPs,'$P_{2r}$[note]'=P2_Raw,'$P_{2b}$'=P2_BLUPs,'$PR_r$[note]'='Parent Range_Raw','$PR_b$[note]'='Parent Range_BLUPs') %>%
        kable(align='llrrrrrrrrrrrrccccrrrrrr', booktabs=TRUE, escape = FALSE, longtable = TRUE, caption=caption) %>%
        row_spec(row=which(tbl1$label != " ")-1, hline_after = TRUE) %>%
        column_spec(column=c(5,7,9,11,13,15,17,19,21), italic=TRUE, color="gray") %>%
        column_spec(column=c(1), bold=TRUE, width="1.0cm") %>%
        column_spec(column=c(2,13,15), width="0.65cm") %>%
        column_spec(column=c(12,14), width="0.75cm") %>%
        column_spec(column=c(16), border_left=TRUE) %>%
        column_spec(column=c(8,10,12,14,16,18,20,22,24), italic=TRUE, color="gray") %>%
        add_header_above(c(" ", " ", "F1 Progeny"=16, "Parents"=6)) %>%
        kable_paper("striped", full_width=FALSE) %>%
        add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                        "Signficance codes for Genotype Effects",
                        "Additive genomic variance of model.",
                        "Additive genomic by year interaction effect variance of model.",
                        "Residual variance of model.",
                        "Narrow-sense genomic heritability for trait/model.",
                        "Range of trait values in progeny.",
                        "Range of BLUP values in progeny.",
                        paste0("F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => ",POP_Name,"_<num>.  eg: g1_77 => ",POP_Name,"_1_77"),
                        "F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for maximum genotype.",
                        "Maternal trait value",
                        "Paternal trait value",
                        "Range (difference) between two parental trait values.",
                        "Range (difference) between two parental BLUP values."
                        ),
                    escape=TRUE)
    if( is_pdf_output() ) {
        tbl2 <- tbl2 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }

    return(tbl2)
}

generateReducedTableFlex <- function(values.collated.tb, caption=NULL) {
    tbl1 <- generateFullTableHelper(values.collated.tb) %>%
                dplyr::select(label,model_label,model_formula,vg,vge,ve,h2,min_Raw,min_BLUPs,mean_Raw,mean_BLUPs,max_Raw,max_BLUPs,range_Raw,range_BLUPs,min_geno_Raw,min_geno_BLUPs,max_geno_Raw,max_geno_BLUPs,P1_Raw,P1_BLUPs,P2_Raw,P2_BLUPs,'Parent Range_Raw','Parent Range_BLUPs') %>%
                select(!model_label)
    tbl2 <- flextable(tbl1) %>%
        align(align="center", part="header") %>%
        align(align="right", part="body") %>%
        align(j =  ~label + model_formula + min_geno_Raw + max_geno_Raw + min_geno_BLUPs + max_geno_BLUPs, align="left", part="body") %>%
        italic(j = grep("_BLUPs", colnames(tbl1), ignore.case=TRUE), part="body") %>%
        bold(i = ~ label != " ", j = 1) %>%
        color(j = grep("_BLUPs", colnames(tbl1), ignore.case=TRUE), color="darkgray", part="body") %>%
        add_header_row(top=FALSE, values=c(" ", " ", "F1 Progeny", "Parents"), colwidths=c(1,1,16,6)) %>%
        flextable::footnote(i = 1, 
                 j = ~vg + vge + ve + h2 + min_geno_Raw + max_geno_Raw + P1_Raw + P2_Raw,
                 value = as_paragraph(c("Additive genomic variance of model.",
                        "Additive genomic by year interaction effect variance of model.",
                        "Residual variance of model.",
                        "Narrow-sense genomic heritability for trait/model.",
                        paste0("F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => ",POP_Name,"_1_<num>.  eg: g1_77 => ",POP_Name,"_1_77"),
                        "F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for maximum genotype.",
                        "Maternal trait value",
                        "Paternal trait value")),
                ref_symbols = c("a","b","c","d","e","f","g","h"),
                part = "header") %>%
        set_caption(caption) %>%
        mk_par(i = 1, j = ~vg + vge + ve + h2 + min_geno_Raw + max_geno_Raw + P1_Raw + P2_Raw, 
               c(as_paragraph("σ",as_sub("g"),as_sup("a")),
                 as_paragraph("σ",as_sub("gε"),as_sup("b")),
                 as_paragraph("σ",as_sub("ε"),as_sup("c")),
                 as_paragraph("h",as_sup("2 d")),
                 as_paragraph("Min Geno",as_sub("r"),as_sup("e")),
                 as_paragraph("Max Geno",as_sub("r"),as_sup("f")),
                 as_paragraph("P",as_sub("1r"),as_sup("g")),
                 as_paragraph("P",as_sub("2r"),as_sup("h")) ),
               part = "header") %>%
        mk_par(i  =  1, 
               j  = c('label','model_formula', 'min_Raw','min_BLUPs','mean_Raw','mean_BLUPs','max_Raw','max_BLUPs','range_Raw','range_BLUPs','min_geno_BLUPs','max_geno_BLUPs','P1_BLUPs','P2_BLUPs','Parent Range_Raw','Parent Range_BLUPs'),
               c(as_paragraph("Trait"),
                 as_paragraph("Model Formula"),
                 as_paragraph("Min",as_sub("r")),
                 as_paragraph("Min",as_sub("b")), 
                 as_paragraph("μ",as_sub("r")),
                 as_paragraph("μ",as_sub("b")),
                 as_paragraph("Max",as_sub("r")),
                 as_paragraph("Max",as_sub("b")),
                 as_paragraph("Range",as_sub("r")),
                 as_paragraph("Range",as_sub("b")),
                 as_paragraph("Min Geno",as_sub("b")),
                 as_paragraph("Max Geno",as_sub("b")),
                 as_paragraph("P",as_sub("1b")),
                 as_paragraph("P",as_sub("2b")),
                 as_paragraph("PRange",as_sub("r")),
                 as_paragraph("PRange",as_sub("b"))),
               part="header") %>%
        mk_par(j =~ model_formula,
               value = as_paragraph(as_equation(., width=0.1, height=0.5, unit="cm")),
               part="body",
               use_dot=TRUE) %>%
        theme_zebra() %>%
        vline(j=18, part="body") %>%
        set_table_properties(width=table.width)
    return(tbl2)
}

generateFullTableFlex <- function(values.collated.tb, caption=NULL) {
    tbl1 <- generateFullTableHelper(values.collated.tb) %>%
                generateTableSignifSymbols() %>%
                generateTableRemoveRepeats() %>%
                mutate(label=gsub("$","",label,fixed=TRUE)) %>%
                mutate(label=gsub(" ",'\\ ',label,fixed=TRUE)) %>%
                mutate(model_label=gsub("$","",model_label,fixed=TRUE)) %>%
                mutate(model_label=gsub(" ",'\\ ',model_label,fixed=TRUE)) %>%
                dplyr::select(label,model_label,model_formula, vg,vge,ve,h2,min_Raw,min_BLUPs,mean_Raw,mean_BLUPs,max_Raw,max_BLUPs,range_Raw,range_BLUPs,min_geno_Raw,min_geno_BLUPs,max_geno_Raw,max_geno_BLUPs,P1_Raw,P1_BLUPs,P2_Raw,P2_BLUPs,'Parent Range_Raw','Parent Range_BLUPs') %>%
                rename(`Trait`=label,`Model`=model_label,`Model Formula`=model_formula, `σ_g^2`=vg,`σ_{gε}^2`=vge,`σ_{ε}^2`=ve,`h^2`=h2,`Min_r`=min_Raw,`Min_b`=min_BLUPs,`Max_r`=max_Raw,`Max_b`=max_BLUPs,`μ_r±SE`=mean_Raw,`μ_b±SE`=mean_BLUPs,`R_r`=range_Raw,`R_b`=range_BLUPs,`MinG_r`=min_geno_Raw,`MinG_b`=min_geno_BLUPs,`MaxG_r`=max_geno_Raw,`MaxG_b`=max_geno_BLUPs,`P_{1r}`=P1_Raw,`P_{1b}`=P1_BLUPs,`P_{2r}`=P2_Raw,`P_{2b}`=P2_BLUPs,`PR_r`=`Parent Range_Raw`,`PR_b`=`Parent Range_BLUPs`)
    tbl2 <- flextable(tbl1) %>%
        set_caption(caption) %>%
        align(align="center", part="header") %>%
        align(align="right", part="body") %>%
        align(j =  colnames(tbl1) %in% c('Trait','Model', 'Model\\ Formula', 'MinG_r','MaxG_r','MinG_b','MaxG_b'), align="left", part="body") %>%
        italic(j = grep("_b", colnames(tbl1), ignore.case=TRUE), part="body") %>%
        bold(i = ~ Trait != " ", j = 1) %>%
        color(j = grep("_b", colnames(tbl1), ignore.case=TRUE), color="gray", part="body") %>%
        mk_par(j = "Trait", part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        mk_par(j = "Model", part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        mk_par(j = 'Model Formula', part="body", value=as_paragraph(as_equation(., width=1, height=0.5)), use_dot=TRUE) %>%
        add_header_row(top=FALSE, values=c(" ", " ", " ", "F1\\ Progeny", "Parents"), colwidths=c(1,1,1,16,6)) %>%
        mk_par(part = "header", value = as_paragraph(as_equation(.,width = .1, height = .2)), use_dot = TRUE) %>%
        flextable::footnote(i = 1, 
                 j = colnames(tbl1) %in% c('Trait','Model','σ_g^2','σ_{gε}^2','σ_{ε}^2','h^2','MinG_r','MaxG_r','P_{1r}','P_{2r}'),
                 value = as_paragraph(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                        "Signficance codes for Genotype Effects",
                        "Additive genomic variance of model.",
                        "Additive genomic by year interaction effect variance of model.",
                        "Residual variance of model.",
                        "Narrow-sense genomic heritability for trait/model.",
                        paste0("F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => ",POP_Name,"_1_<num>.  eg: g1_77 => ",POP_Name,"_1_77"),
                        "F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for maximum genotype.",
                        "Maternal trait value",
                        "Paternal trait value")),
                ref_symbols = c("a","b","c","d","e","f","g","h","i","j"),
                part = "header") %>%
        theme_zebra() %>%
        hline( i = (which(tbl1$Trait != " ")-1)[-1], part="body") %>%
        set_table_properties(width=table.width)
    return(tbl2)
}

generateReducedTable <- function(values.collated.tb, caption=NULL, tfont_size=10) {
    tbl1 <- generateFullTableHelper(values.collated.tb) %>%
                mutate(model_formula = gsub("^","$",model_formula,fixed=FALSE)) %>%
                mutate(model_formula = gsub("$","$",model_formula,fixed=FALSE)) %>%
                arrange(desc(h2)) %>% 
                dplyr::select(label,model_formula,h2,min_Raw,min_BLUPs,mean_Raw,mean_BLUPs,max_Raw,max_BLUPs,range_Raw,range_BLUPs,bottom_genos_Raw,bottom_genos_BLUPs,top_genos_Raw,top_genos_BLUPs,P1_Raw,P1_BLUPs,P2_Raw,P2_BLUPs,'Parent Range_Raw','Parent Range_BLUPs')
    tbl2 <- tbl1 %>%
        rename("Trait"=label,'Model Terms[note]'=model_formula,'$h^{2}$[note]'=h2,'$Min_r$'=min_Raw,'$Min_b$'=min_BLUPs,'$Max_r$'=max_Raw,'$Max_b$'=max_BLUPs,"$μ_r$±$SE$"=mean_Raw,"$μ_b$±$SE$"=mean_BLUPs,'$R_r$[note]'=range_Raw,'$R_b$[note]'=range_BLUPs,'$MinG_r$[note]'=bottom_genos_Raw,'$MinG_b$'=bottom_genos_BLUPs,'$MaxG_r$[note]'=top_genos_Raw,'$MaxG_b$'=top_genos_BLUPs,'$P_{1r}$[note]'=P1_Raw,'$P_{1b}$'=P1_BLUPs,'$P_{2r}$[note]'=P2_Raw,'$P_{2b}$'=P2_BLUPs,'$PR_r$[note]'='Parent Range_Raw','$PR_b$[note]'='Parent Range_BLUPs') %>%
        kable(align='lcrrrrrrrrrccccrrrrrr', booktabs=TRUE, escape = FALSE, longtable = TRUE, caption=caption) %>%
        column_spec(column=c(5,7,9,11,13,15,17,19,21), italic=TRUE, color="gray") %>%
        column_spec(column=c(1), bold=TRUE, width="1.0cm") %>%
        column_spec(column=c(2,13,15), width="0.65cm") %>%
        column_spec(column=c(12,14), width="0.75cm") %>%
        column_spec(column=c(16), border_left=TRUE) %>%
        add_header_above(c(" ", " ", "F1 Progeny"=13, "Parents"=6)) %>%
        kable_paper("striped", full_width=FALSE) %>%
        add_footnote(c("Optimum model selected using AIC criterium.",
                        "Narrow-sense genomic heritability for trait/model.",
                        "Range of trait values in progeny.",
                        "Range of BLUP values in progeny.",
                        paste0("F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => ",POP_Name,"_1_<num>.  eg: g1_77 => ",POP_Name,"_1_77"),
                        "F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for minimum genotype.",
                        "Maternal trait value",
                        "Paternal trait value",
                        "Range (difference) between two parental trait values.",
                        "Range (difference) between two parental BLUP values."
                        ),
                     escape=TRUE)
    if( is_pdf_output() ) {
        tbl2 <- tbl2 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  font_size=tfont_size,
                                  repeat_header_continued = TRUE)
    }
    return(tbl2)
}

generateTable <- function(values.collated.tb, type, caption=NULL, tfont_size=10) {
    print(tbl2)
    tbl1 <- values.collated.tb %>%
        arrange(type,trait,model) %>%
        mutate(label=paste0(label,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$")) %>%
        mutate(model_label=paste0(model_label,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) %>%
        mutate(trait_repeat=(label == c("",label[-length(label)])),
                model_repeat=(model_label == c("",model_label[-length(model_label)]))) %>%
        mutate(label = ifelse(trait_repeat, " ", cell_spec(label, bold=TRUE)),
               model_label = ifelse(model_repeat, " ", cell_spec(model_label, bold=TRUE))) %>%
        mutate(min=signif.digits(min,3)) %>%
        mutate(max=signif.digits(max,3)) %>%
        mutate(mean=ifelse(!is.na(mean),paste0(signif.digits(mean,3),"±",signif.digits(sd,3))," ")) %>% 
        mutate(vg=signif.digits(vg,3)) %>%
        mutate(vge=signif.digits(vge,3)) %>%
        mutate(ve=signif.digits(ve,3)) %>%
        mutate(h2=signif.digits(h2,3)) %>%
        mutate(range=signif.digits(range,3)) %>%
        mutate(P1=signif.digits(P1,3)) %>%
        mutate(P2=signif.digits(P2,3)) %>%
        dplyr::select(label,model_label,vg,vge,ve,h2,min,max,mean,range,min_geno,max_geno,P1,P2,'Parent Range')
    tbl2 <- tbl1 %>%
        rename("Trait[note]"=label,"Model[note]"=model_label,'$σ_{g}^{2}$[note]'=vg,'$σ_{g ε}^{2}$[note]'=vge,'$σ_{ε}^{2}$[note]'=ve,'$h^{2}$[note]'=h2,Min=min,Max=max,"μ±SE"=mean,Range=range,'MinG[note]'=min_geno,'Max\ Geno[note]'=max_geno,'P1[note]'=P1,'P2[note]'=P2) %>%
        kable(align='llrrrrrrrrccrrr', booktabs=TRUE, escape = FALSE, longtable = TRUE, caption=caption) %>%
        row_spec(row=which(values.collated.tb$trait != " ")-1, hline_after = TRUE) %>%
        add_header_above(c(" ", " ", "F1 Progeny"=10, "Parents"=3)) %>%
        kable_paper("striped", full_width=FALSE) %>%
        add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                        "Signficance codes for Genotype Effects",
                        "Additive genomic variance of model.",
                        "Additive genomic by year interaction effect variance of model.",
                        "Residual variance of model.",
                        "Narrow-sense genomic heritability for trait/model.",
                        paste0("F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => ",POP_Name,"_1_<num>.  eg: g1_77 => ",POP_Name,"_1_77"),
                        "F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for maximum genotype.",
                        "Maternal trait value",
                        "Paternal trait value"),
                     escape=TRUE)
    if( is_pdf_output() ) {
        tbl2 <- tbl2 %>%
                    kable_styling(latex_options=c("repeat_header"),
                                  repeat_header_continued = TRUE,
                                  font_size=tfont_size)
    }
    summary.path <- normalizePath(paste0(workflow,'/traits/'), mustWork = TRUE);
    summary.file <- paste0(summary.path,'/',type,'.summary.table.png')
    print(tbl2)
    return(tbl2)
}

save.image(paste0(workflow,"/.RData.01_genBLUP"))

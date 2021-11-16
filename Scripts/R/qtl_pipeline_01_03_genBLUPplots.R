#!/usr/bin/env RScript

#NOTE: This particular script generates informative data plots for the BLUPs: Things like pairwise correlation plots, boxplots, etc.
# Think of combining BLUP data into a table format similar to the input phenotype table, and then can display data in a way similar
# to how this was done for the original phenotypes.
#
# GGally is our friend here!

# loading libraries
library(cowplot)
library(formattable)
library(GGally)
library(ggplot2)
library(ggthemes)
library(jsonlite)
library(kableExtra)
library(knitr)
library(RColorBrewer)
library(tidyverse)

source('./usefulFunctions.R')

workflow <- get0("workflow", ifnotfound="../../Workflows/1")
P1_Name <- get0("P1_Name", ifnotfound="Mullica_Queen")
P2_Name <- get0("P2_Name", ifnotfound="Crimson_Queen")

source(paste0(workflow,"/configs/model.cfg"))

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

collateBLUPs <- function(trait.cfg, trait.path, args.l) {
    models.tb.p         <- args.l[[2]]
	num_years			<- args.l[[3]]
    trait_name          <- trait.cfg$trait
    model_label         <- trait.cfg$model
    cat(paste0("Processing trait: ", trait_name, ", Model: ", model_label,"\n"))
    model <- readRDS(paste0(trait.path,"/mmer.rds"))
    varcomp	  <- model$sigma
    g.idx     <- which(grepl("^u:id",names(varcomp)))
    ge.idx    <- which(grepl("^id:year",names(varcomp)))
    e.idx	  <- which(grepl("^units",names(varcomp)))
    vg        <- varcomp[[g.idx]][trait_name,trait_name]
    vge       <- ifelse(length(ge.idx) > 0,varcomp[[ge.idx]][trait_name,trait_name],NA)
    ve        <- varcomp[[e.idx]][trait_name,trait_name]
    if( model_label == "all-years" ) {
		if(is.na(vge)) {
			h2 <- vg/(vg+(ve/num_years))
		} else {
			h2 <- vg/(vg+(vge/num_years)+(ve/num_years))
		}
    } else {
		h2 <- vg/(vg+ve)
    }
    u                    <- model$U[["u:id"]]
    intercept            <- (model$Beta %>% filter(Effect == "(Intercept)"))$Estimate
    blups                <- u[[trait_name]] + intercept
    blupsZ                <- scale(blups)[,1]
    anova.file            <- paste0(trait.path,'/anova.csv')
    GLRpvalue            <- NA
    GZRpvalue            <- NA
    GxYLRpvalue            <- NA
    GxYZRpvalue            <- NA
    if( file.exists(anova.file) ) {
        anova.df=read.csv(file=anova.file, header=TRUE, row.names=1)
        GLRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% dplyr::select(PrChisq))
        GZRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% dplyr::select(PrNorm))
        if( model_label == "all-years" ) {
            GxYLRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% dplyr::select(PrChisq))
            GxYZRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% dplyr::select(PrNorm))
        }
    }
    if( model_label == "all-years" ) {
        pheno <- model$dataOriginal[,c("id",trait_name)] %>% group_by(id) %>% dplyr::summarize(mean=mean(.data[[trait_name]]))
        pheno <- pheno[["mean"]]
    } else {
        pheno <- model$dataOriginal[[trait_name]]
    }
    phenoZ        <- scale(pheno)[,1]
    ids           <- names(blups)
    models.tb.p$value <- models.tb.p$value %>% 
                            add_row(id=factor(ids), model=factor(trait.cfg$model), model_label=factor(trait.cfg$model_label), trait=factor(trait.cfg$trait), value=pheno, valueZ=phenoZ, type=factor("Raw"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2) %>%
                            add_row(id=factor(ids), model=factor(trait.cfg$model), model_label=factor(trait.cfg$model_label), trait=factor(trait.cfg$trait), value=blups, valueZ=blupsZ, type=factor("BLUPs"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2) 
}

pheno.means.df <-read.csv(file=pheno_dpath2fpath(pheno_file))
pheno.means.df$year <- as.factor(pheno.means.df$year) #Needed for modeling column effects
n.years <- nlevels(pheno.means.df$year)
traits.df   <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=T)
#Just build a 'long' version of collated BLUPs as tibble table, and then use tidyverse 'pivot_wider()' to flatten
model.collated.long.tb	   <- tibble(id=factor(), model=factor(), model_label=factor(), trait=factor(), value=numeric(), valueZ=numeric(), type=factor(), label=factor(), label_short=factor(), GLRpvalue=numeric(), GZRpvalue=numeric(), GxYLRpvalue=numeric(), GxYZRpvalue=numeric(), vg=numeric(), vge=numeric(), ve=numeric(), h2=numeric())
model.collated.long.tb.p   <- newPointer(model.collated.long.tb)
selectedmodels.df          <- read.csv(file=paste0(workflow,"/traits/selectedModels.csv"))
loopThruTraits(workflow, collateBLUPs, loopArgs=list(selectedmodels.df,model.collated.long.tb.p,n.years))
model.collated.long.tb     <- model.collated.long.tb.p$value %>%
                                group_by(trait) %>%
                                mutate(GxYLRpvalue = rep(min(GxYLRpvalue,na.rm=TRUE),length(GxYLRpvalue))) %>% #Modify to have the interaction effect modeled in 'all-years' for all model,trait groupings
                                mutate(GxYZRpvalue = rep(min(GxYZRpvalue,na.rm=TRUE),length(GxYZRpvalue))) %>%
                                ungroup()
#Save long version, as it can be manipulated in many ways to produce meaningful plots/graphs
saveRDS(model.collated.long.tb, file=paste0(workflow,'/traits/blups_collated.long.rds'), compress=TRUE)
model.collated.long.tb <- readRDS(file=paste0(workflow,'/traits/blups_collated.long.rds'))
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

write_json(model.collated.wide.tb, paste0(workflow, '/traits/blups_collated.long.json'), auto_unbox=T, pretty=T) #Also write to JSON file for easy import by other programming languages/tools

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
                                    dplyr::summarize(min=min(value), minZ=min(valueZ), min_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.min(value)]),
                                              max=max(value), maxZ=max(valueZ), max_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.max(value)]),
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
           fake = for r in .$data[[1]] {
           }
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

#Generate a table of BLUP summaries by model and trait, with min value, min genotype, max value, max genotype, MQ blup value, CQ blup value

p1.blup_summary.tb <- model.collated.long.tb %>%
                        filter(id == P1_Name) %>%
                        group_by(model,trait,type) %>%
                        dplyr::summarize(MQ = mean(value), MQz = mean(value))

p2.blup_summary.tb <- model.collated.long.tb %>%
                        filter(id == P2_Name) %>%
                        group_by(model,trait,type) %>%
                        dplyr::summarize(CQ = mean(value), CQz = mean(value))
model.collated.summary.tb <- model.collated.long.tb %>%
                                    filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                    group_by(model,trait,type) %>%
                                    dplyr::summarize(min=min(value), mean=mean(value), sd=sd(value), range=abs(min(value)-max(value)),
                                              minZ=min(valueZ), min_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.min(value)]),
                                              max=max(value), maxZ=max(valueZ), max_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.max(value)]),
                                              GLRpvalue=mean(GLRpvalue),GZRpvalue=mean(GZRpvalue),GxYLRpvalue=mean(GxYLRpvalue),GxYZRpvalue=mean(GxYZRpvalue),
											  vg=mean(vg), vge=mean(vge), ve=mean(ve), h2=mean(h2)) %>%
                                    filter(min != max) %>%
                                    left_join(p1.blup_summary.tb, by=c("model","trait","type")) %>%
                                    left_join(p2.blup_summary.tb, by=c("model","trait","type")) %>%
                                    mutate('Parent Range' = round.digits(abs(CQ-MQ),2)) %>%
                                    ungroup() #This is necessary for trait_repeat and model_repeat calculations to work below

model.collated.blup_summary.tb <- model.collated.summary.tb %>%
        filter(type == "BLUPs")

model.collated.raw_summary.tb <- model.collated.summary.tb %>%
        filter(type == "Raw")


generateTable <- function(values.collated.tb, type) {
	tbl <- values.collated.tb %>%
			arrange(type,trait,model) %>%
			mutate(trait=paste0(trait,"$^{",unlist(map(GxYLRpvalue,siginfo)),"}$")) %>%
			mutate(model=paste0(model,"$^{",unlist(map(GLRpvalue,siginfo)),"}$")) %>%
			mutate(trait_repeat=(trait == c("",trait[-length(trait)])),
					model_repeat=(model == c("",model[-length(model)]))) %>%
			mutate(trait = ifelse(trait_repeat, "", cell_spec(trait, "html", bold=TRUE)),
					model = ifelse(model_repeat, "", cell_spec(model, "html", bold=TRUE))) %>%
			mutate(min=round.digits(min,2)) %>%
			mutate(max=round.digits(max,2)) %>%
			mutate(mean=paste0(round.digits(mean,2),"\u00b1",round.digits(sd,2))) %>%
			mutate(vg=round.digits(vg,2)) %>%
			mutate(vge=round.digits(vge,2)) %>%
			mutate(ve=round.digits(ve,2)) %>%
			mutate(h2=round.digits(h2,3)) %>%
			mutate(range=round.digits(range,2)) %>%
			mutate(MQ=round.digits(MQ,2)) %>%
			mutate(CQ=round.digits(CQ,2)) %>%
			dplyr::select(trait,model,vg,vge,ve,h2,min,max,mean,range,min_geno,max_geno,MQ,CQ,'Parent Range') %>%
			rename("Trait[note]"=trait,"Model[note]"=model,'$\\sigma_{g}^{2}$[note]'=vg,'$\\sigma_{g \\epsilon}^{2}$[note]'=vge,'$\\sigma_{\\epsilon}^{2}$[note]'=ve,'$h^{2}$[note]'=h2,Min=min,Max=max,"Mean \u00b1 SE"=mean,Range=range,'Min Geno[note]'=min_geno,'Max Geno[note]'=max_geno,'MQ[note]'=MQ,'CQ[note]'=CQ) %>%
			kable("html", align='llrrrrrrrrccrrr', booktabs=TRUE, escape = FALSE, table.attr="id=\"kableTable\"", caption=ifelse(type=="raw","Raw Phenotype Statistics for Trait Models", "BLUP Statistics for Trait Models")) %>%
			row_spec(row=which(values.collated.tb$trait != ""), extra_css = "border-top: 1px solid #ddd") %>%
			add_header_above(c("", "", "F1 Progeny"=10, "Parents"=3)) %>%
			kable_paper("striped", full_width=FALSE) %>%
			add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
						"Signficance codes for Genotype Effects",
						"Additive genomic variance of model.",
						"Additive genomic by year interaction effect variance of model.",
						"Residual variance of model.",
						"Narrow-sense genomic heritability for trait/model.",
						"F1 progeny genotype with minimum trait value.  Genotype identifier is shortened for visibility.  Translation is g<num> => CNJ02_1_<num>.  eg: g77 => CNJ02_1_77",
						"F1 progeny genotype with maximum trait value.  Genotype identifier is same format as for maximum genotype.",
						"Maternal Mullica Queen trait value",
						"Paternal Crimson Queen trait value"))
  summary.path <- normalizePath(paste0(workflow,'/traits/'), mustWork = TRUE);
  summary.file <- paste0(summary.path,'/',type,'.summary.table.png')
	cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",summary.file,"'"))
	print(tbl)
}

generateTable(model.collated.raw_summary.tb, "raw")
generateTable(model.collated.blup_summary.tb, "blups")
generateTable(model.collated.blup_summary.tb %>% filter(trait %in% c("Berry Length","Berry Width", "Berry Weight")), "blups")
generateTable(model.collated.blup_summary.tb %>% filter(trait %in% c("Number of Berries", "Total Berry Weight")), "blups")

#summary.path <- normalizePath(paste0(workflow,'/traits/'), mustWork = TRUE);
#summary.file <- paste0(summary.path,'/blups.summary.table.png');
#cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",summary.file,"'"))

save.image(paste0(workflow,"/.RData.01_genBLUP"))

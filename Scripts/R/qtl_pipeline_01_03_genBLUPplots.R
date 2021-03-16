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
library(kableExtra)
library(knitr)
library(RColorBrewer)
library(tidyverse)

source(paste0(workflow,"/configs/model.cfg"))
source('./usefulFunctions.R')

workflow <- "../../Workflows/1"
P1_Name <- "Mullica_Queen"
P2_Name <- "Crimson_Queen"

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

collateBLUPs <- function(trait.cfg, trait.path, args.l) {
    selectedmodels.df   <- args.l[[1]]
    models.tb.p         <- args.l[[2]]
	num_years			<- args.l[[3]]
    trait_name          <- trait.cfg$trait
    model_label         <- trait.cfg$model
    selectedmodel.label <- (selectedmodels.df %>% filter((trait == trait_name) & (model == model_label)) %>% select(selected_model))[[1]] #Choose the correct model
    anova.drop1.selectedmodel.file <- paste0(trait.path,"/anova.models.drop1.rds")
    if( file.exists(anova.drop1.selectedmodel.file) ) {
		anova.drop1.selectedmodel.l <- readRDS(anova.drop1.selectedmodel.file)
		model						<- anova.drop1.selectedmodel.l[[selectedmodel.label]]
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
            GLRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% select(PrChisq))
            GZRpvalue    <- as.numeric(anova.df %>% filter(dropterm == 'vs(id, Gu = A)') %>% select(PrNorm))
            if( model_label == "all-years" ) {
                GxYLRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% select(PrChisq))
                GxYZRpvalue <- as.numeric(anova.df %>% filter(dropterm == 'id:year') %>% select(PrNorm))
            }
        }
        if( model_label == "all-years" ) {
            pheno <- model$data[,c("id",trait_name)] %>% group_by(id) %>% summarize(mean=mean(.data[[trait_name]]))
            pheno <- pheno[["mean"]]
        } else {
            pheno <- model$data[[trait_name]]
        }
        phenoZ        <- scale(pheno)[,1]
        ids           <- names(blups)
        model_label   <- ifelse(trait.cfg$model == 'all-years', 'All Years', trait.cfg$model)
        models.tb.p$value <- models.tb.p$value %>% 
                                add_row(id=factor(ids), model=factor(trait.cfg$model, labels=model_label), trait=factor(trait_name, labels=factor(trait.cfg$label)), value=blups, valueZ=blupsZ, type=factor("BLUPs"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2)
        models.tb.p$value <- models.tb.p$value %>% 
                                add_row(id=factor(ids), model=factor(trait.cfg$model, labels=model_label), trait=factor(trait_name, labels=factor(trait.cfg$label)), value=pheno, valueZ=phenoZ, type=factor("Raw"), label=factor(trait.cfg$label), label_short=factor(trait.cfg$label_short), GLRpvalue=GLRpvalue, GZRpvalue=GZRpvalue, GxYLRpvalue=GxYLRpvalue, GxYZRpvalue=GxYZRpvalue, vg=vg, vge=vge, ve=ve, h2=h2) 
    }
}

pheno.means.df <-read.csv(file=pheno_dpath2fpath(pheno_file))
pheno.means.df$year <- as.factor(pheno.means.df$year) #Needed for modeling column effects
n.years <- nlevels(pheno.means.df$year)
traits.df   <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=T)
#Just build a 'long' version of collated BLUPs as tibble table, and then use tidyverse 'pivot_wider()' to flatten
model.collated.long.tb	   <- tibble(id=factor(), model=factor(), trait=factor(), value=numeric(), valueZ=numeric(), type=factor(), label=factor(), label_short=factor(), GLRpvalue=numeric(), GZRpvalue=numeric(), GxYLRpvalue=numeric(), GxYZRpvalue=numeric(), vg=numeric(), vge=numeric(), ve=numeric(), h2=numeric())
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
                            mutate(model_trait=as.factor(paste0(model,' ',trait)), model_type=as.factor(paste0(model,' ',type)))

model.collated.long.stats.tb <- model.collated.long.tb %>%
                                    filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                    group_by(model,trait,type,model_type) %>%
                                    summarize(min=min(value), minZ=min(valueZ), min_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.min(value)]),
                                              max=max(value), maxZ=max(valueZ), max_geno = gsub("CNJ02_1_([0-9]+)","g\\1",id[which.max(value)]),
                                              .groups="keep") %>%
                                    filter(min != max) %>%
                                    ungroup() %>%
                                    group_nest(trait, .key="stats")

model.collated.long.ylims.tb <- model.collated.long.tb %>%
                                    filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                    group_by(trait,type) %>%
                                    filter(min(value)!= max(value)) %>%
                                    summarize(min=min(value, na.rm=TRUE) - 2.5, minZ=min(valueZ, na.rm=TRUE)-1,
                                              max=max(value, na.rm=TRUE) + 2.5, maxZ=max(valueZ, na.rm=TRUE)+1,
                                              .groups="keep") %>%
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
                    ggtitle(label=.$trait) +
                    theme(    axis.text.x = element_text(face="bold", size=32, angle = 60, hjust = 1),
                            axis.text.y = element_text(face="bold", size=32),
                            legend.title = element_text(fac="bold", size=32),
                            legend.text = element_text(fac="bold", size=24),
                            legend.key.size = ggplot2::unit(5,"cm"),
                            axis.title  = element_text(face="bold",size=48),
                            plot.title   = element_text(face="bold",size=52, hjust=0.5),
                            plot.subtitle = element_text(size=48, hjust = 0.5),
                            plot.margin = ggplot2::unit(c(1,1,1,1),"cm")))
                        
p1 <- (gs %>% filter(trait == "Berry Length"))$plot[[1]] + ylab("Value")
p2 <- (gs %>% filter(trait == "Berry Weight"))$plot[[1]]
p3 <- (gs %>% filter(trait == "Berry Width"))$plot[[1]]
p4 <- (gs %>% filter(trait == "Total Berry Weight"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
#p3 <- (gs %>% filter(trait == "Berry Width"))$plot[[1]] + guides(fill=guide_legend(title="Type of Value"))
#p4 <- (gs %>% filter(trait == "Total Berry Weight"))$plot[[1]]
p5 <- (gs %>% filter(trait == "Number of Pedicels"))$plot[[1]] + ylab("Value")
p6 <- (gs %>% filter(trait == "Number of Berries"))$plot[[1]]
p7 <- (gs %>% filter(trait == "Number of Seeds"))$plot[[1]]

#pg <- plot_grid(p1,p2,p3,NULL,p4,NULL,p5,p6,p7,nrow=3,ncol=3,rel_widths=c(1,1,1.25))
pg <- plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4,rel_widths=c(1,1,1,1.25))

#png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.plotgrid.png'), width=2560, height=3840, bg="white")
png(filename=paste0(workflow,'/traits/plots/blups_collated.boxplot.plotgrid.wide.png'), width=3840, height=2560, bg="white")
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
                        summarize(MQ = mean(value), MQz = mean(value))

p2.blup_summary.tb <- model.collated.long.tb %>%
                        filter(id == P2_Name) %>%
                        group_by(model,trait,type) %>%
                        summarize(CQ = mean(value), CQz = mean(value))
model.collated.summary.tb <- model.collated.long.tb %>%
                                    filter(!(id %in% c(P1_Name,P2_Name))) %>%
                                    group_by(model,trait,type) %>%
                                    summarize(min=min(value), mean=mean(value), sd=sd(value), range=abs(min(value)-max(value)),
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
			select(trait,model,vg,vge,ve,h2,min,max,mean,range,min_geno,max_geno,MQ,CQ,'Parent Range') %>%
			rename("Trait[note]"=trait,"Model[note]"=model,'$\\sigma_{g}^{2}$[note]'=vg,'$\\sigma_{g \\epsilon}^{2}$[note]'=vge,'$\\sigma_{\\epsilon}^{2}$[note]'=ve,'$h^{2}$[note]'=h2,Min=min,Max=max,"Mean \u00b1 SE"=mean,Range=range,'Min Geno[note]'=min_geno,'Max Geno[note]'=max_geno,'MQ[note]'=MQ,'CQ[note]'=CQ) %>%
			kable("html", align='llrrrrrrrrccrrr', booktabs=TRUE, escape = FALSE, table.attr="id=\"kableTable\"") %>%
			row_spec(row=which(values.collated.tb$trait != ""), extra_css = "border-top: 1px solid #ddd") %>%
			add_header_above(c("", "", "F1 Progeny"=10, "Parents"=3)) %>%
			kable_paper("striped", full_width=FALSE) %>%
			add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
						"Signficance codes for Genotype Effects",
						"Additive genomic variance of model.",
						"Additive genomic by year interaction effect variance of model.",
						"Residual variance of model.",
						"Narrow-sense genomic heritability for trait/model.",
						"F1 progeny genotype with minimum trait BLUP value.  Genotype identifier is shortened for visibility.  Translation is g<num> => CNJ02_1_<num>.  eg: g77 => CNJ02_1_77",
						"F1 progeny genotype with minimum trait BLUP value.  Genotype identifier is same format as for minimum genotype.",
						"Maternal Mullica Queen trait BLUP value - Data only available for years 2012, 2013",
						"Paternal Crimson Queen trait BLUP value - Data only available for years 2012, 2013"))

	cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/',type,'.summary.table.png'),mustWork=TRUE),"'"))
	print(tbl)
}

generateTable(model.collated.raw_summary.tb, "raw")
generateTable("blups")
generateTable(model.collated.blup_summary.tb %>% filter(trait %in% c("Berry Length","Berry Width", "Berry Weight")), "blups")
generateTable(model.collated.blup_summary.tb %>% filter(trait %in% c("Number of Berries", "Total Berry Weight")), "blups")

model.collated.blup_summary.tb %>%
    select(trait,model,vg,vge,ve,h2,min,max,mean,range,min_geno,max_geno,MQ,CQ,'Parent Range') %>%
    rename("Trait[note]"=trait,"Model[note]"=model,'$\\sigma_{g}^{2}$'=vg,'$\\sigma_{g \\epsilon}^{2}$'=vge,'$\\sigma_{\\epsilon}^{2}$'=ve,'$h^{2}$'=h2,Min=min,Max=max,"Mean \u00b1 SE"=mean,Range=range,'Min Geno[note]'=min_geno,'Max Geno[note]'=max_geno,'MQ[note]'=MQ,'CQ[note]'=CQ) %>%
    kable("html", align='llrrrrrrrrccrrr', booktabs=TRUE, escape = FALSE, table.attr="id=\"kableTable\"", caption=ifelse(type=="raw","Raw Phenotype Statistics for Trait Models", "BLUP Statistics for Trait Models")) %>%
    row_spec(row=which(model.collated.blup_summary.tb$trait != ""), extra_css = "border-top: 1px solid #ddd") %>%
    add_header_above(c("", "", "F1 Progeny"=10, "Parents"=3)) %>%
    kable_paper("striped", full_width=FALSE) %>%
    add_footnote(c("Significance codes for Genotype:Year Effects\n*** pvalue≥0 and pvalue<0.001\n**  pvalue≥0.001 and pvalue<0.01\n*   pvalue≥0.01 and pvalue<0.05\n.  pvalue≥0.05 and pvalue<0.01\nNS  Not Significant\n", 
                   "Signficance codes for Genotype Effects",
                   "F1 progeny genotype with minimum trait BLUP value.  Genotype identifier is shortened for visibility.  Translation is g<num> => CNJ02_1_<num>.  eg: g77 => CNJ02_1_77",
                   "F1 progeny genotype with minimum trait BLUP value.  Genotype identifier is same format as for minimum genotype.",
                   "Maternal Mullica Queen trait BLUP value - Data only available for years 2012, 2013",
                   "Paternal Crimson Queen trait BLUP value - Data only available for years 2012, 2013"))


cat(paste0("In Firefox javascript console, type: ':screenshot --dpi 8 --file --selector #kableTable --filename ",normalizePath(paste0(workflow,'/traits/blups.summary.table.png'),mustWork=TRUE),"'"))

save.image(".RData.01_genBLUP")

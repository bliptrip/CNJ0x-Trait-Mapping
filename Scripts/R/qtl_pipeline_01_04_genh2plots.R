#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script generates the correlation plots derived from means of two Vorsa populations along with the correlations of 
# both populations combined..

# loading libraries
source('./usefulFunctions.R')

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

workflow <- get0("workflow", ifnotfound="../../Workflows/1")
POP_Name <- get0("POP_Name", ifnotfound="cnj02")
show_groups <- get0("show_groups", ifnotfound=T)
show_legend <- get0("show_legend", ifnotfound=T)

show_groups <- as.logical(show_groups)
show_legend <- as.logical(show_legend)

library(RColorBrewer)
library(corrplot)
library(GGally)
library(ggthemes)
library(lattice)
library(plotly)
library(tidyverse)
library(ggh4x)

blups.collated.tb <- readRDS(paste0(workflow,"/traits/blups_collated.wide.rds"))

h2.tb <- blups.collated.tb %>%
            group_by(model,trait,label,label_short) %>%
            summarise(h2=mean(h2)) %>%
            mutate(se=0) #h2 SE is not valid in current paradigm
h2.collated.path <- paste0(workflow, "/traits/h2.csv")
write_csv(h2.tb, file=h2.collated.path)

decorate_h2_label <- function(model, h2) {
    model_h2            <- vector("character", length(model))
    nul.tf              <- (h2 == 0) | is.na(h2)
    model_h2[nul.tf]   <- as.character(model[nul.tf])
    model_h2[!nul.tf]  <- paste0(model[!nul.tf],": ",h2[!nul.tf])
    return(model_h2)
}

generateH2Plot <- function(h2.data, traits, subtitle, display_groups=T,display_legend=T,model_levels=c("2011","2012","2013","2014","all-years"),model_labels=c("2011","2012","2013","2014","All Years")) {
    palette <- brewer.pal(n=length(model_labels),"RdYlBu")
    names(palette) <- model_labels
    h2.data <- h2.data %>% filter(label_short %in% traits)
    #h2.data$model <- as.factor(h2.data$model) #Cast to factor to have plots correctly treat as discrete
    #h2.data$label <- as.factor(h2.data$label)
    #h2.data$label_short <- as.factor(h2.data$label_short)
    h2.nul.idx <- which((h2.data$h2
                         == 0) | (h2.data$se == 0))
    h2.data$se[h2.nul.idx] <- NA
        
    h2.d2 <- h2.data %>% 
             mutate(trait=factor(trait,levels=unique(trait),ordered=T),
                    group=factor(label_short,levels=unique(label_short),ordered=T),
                    model=factor(model,levels=model_levels,labels=model_labels,ordered=T),
                    value=h2) %>%
             select(trait,group,model,value,se)
    
    h2.complete <- data.frame( trait=rep(levels(h2.d2$trait),each=nlevels(h2.d2$model)),
                               group=rep(levels(h2.d2$group),each=nlevels(h2.d2$model)),
                               model=rep(levels(h2.d2$model),times=nlevels(h2.d2$group)),
                               value=NA,
                               se=NA) %>%
                    mutate( group_model = paste0(group,"-",model) ) %>%
                    mutate( group_model = factor(group_model,levels=unique(group_model),ordered=T),
                            group = factor(group,levels=levels(h2.d2$group),ordered=T),
                            model = factor(model,levels=model_levels,labels=model_labels,ordered=T) )
    h2.d2 <- h2.d2 %>% 
             mutate( group_model = paste0(group,"-",model) )
    h2.amend <- h2.complete %>% 
                    filter( !(group_model %in% unique(h2.d2$group_model)) )
    data <- h2.d2 %>% bind_rows(h2.amend) 
    data$value[data$value < 0.1] <- NA
    data$h2 <- as.character(round(data$value,2))
    data <- data %>% 
                mutate(model_h2=decorate_h2_label(model,h2)) %>%
                mutate(model_trait=paste0(model,"-",group),
                       model_h2=factor(model_h2,levels=unique(model_h2),ordered=T)) %>%
                group_by(group) %>%
                arrange(desc(model),.by_group=T)
    #data$value[is.na(data$value)] <- 0 #Forces empty datasets to render on barplots
    
    
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar=3
    to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
    colnames(to_add) = colnames(data)
    to_add$group=rep(levels(data$group), each=empty_bar)
    data=rbind(data, to_add)
    data=data %>% arrange(group)
    data$id=seq(1, nrow(data))
    
    # Get the name and the y position of each label
    label_data=data
    number_of_bar=nrow(label_data)
    angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data=data %>% 
        group_by(group) %>% 
        summarize(start=min(id), end=max(id) - empty_bar) %>% 
        rowwise() %>% 
        mutate(title=mean(c(start, end)))
    
    # prepare a data frame for grid (scales)
    grid_data = base_data
    grid_data$start = base_data$end + 1
    grid_data$end = grid_data$start + (empty_bar-1)
    
    # Make the plot
    p = ggplot(data, aes(x=as.factor(id), y=value, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
            geom_bar(stat="identity", color="black", alpha=0.9) +
            scale_fill_manual(values=palette) +
            #geom_errorbar(aes(x=as.factor(id), ymin=value-(se/2), ymax=value+(se/2))) +
            #geom_pointrange(aes(x=as.factor(id), y=value, ymin=value-(se/2), ymax=value+(se/2)), na.rm=TRUE) +
            # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
            geom_segment(data=grid_data, aes(x = start, y = 1.0, xend = end, yend = 1.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = start, y = 0.75, xend = end, yend = 0.75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = start, y = 0.5, xend = end, yend = 0.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = start, y = 0.25, xend = end, yend = 0.25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            # Add text showing the value of each 100/75/50/25 lines
            annotate("text", x = rep(max(data$id),4), y = c(0.25, 0.5, 0.75, 1.0), label = c("0.25", "0.50", "0.75", "1.00") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
            ylim(-0.5,1.25) +
            theme_minimal() +
            #theme_solarized(light=F) +
            theme(
                    legend.position = "none",
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.margin = unit(rep(-1,4), "cm") 
                ) +
            coord_polar() + 
            geom_text(data=label_data, aes(x=id, y=value+0.1, label=model_h2, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
            # Add base line information
            geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), colour = "black", alpha=0.9, size=0.6 , inherit.aes = FALSE )  +
            geom_text(data=base_data, aes(x = title, y = -.2, label=group), hjust=0.5, colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    
    ggsave(filename=paste0(workflow,"/traits/plots/h2_barplot.polar.",subtitle,".png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

    #Make a non-polar version of plot
    #scale_x_discrete(drop=F) is needed to avoid dropping model factors that have an NA for a particular trait (group)
    p = ggplot(data %>% mutate(id=-1*id), aes(x=id, y=value, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    		geom_bar(stat="identity", color='black', width=0.65, alpha=0.9, show.legend=display_legend) +
            scale_fill_manual(values=palette, drop=FALSE) +
            scale_x_discrete(drop=FALSE) +
            geom_segment(data=grid_data, aes(x = -start, y = 1.0, xend = -end, yend = 1.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = -start, y = 0.75, xend = -end, yend = 0.75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = -start, y = 0.5, xend = -end, yend = 0.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = -start, y = 0.25, xend = -end, yend = 0.25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            annotate("text", x = -1*rep(max(grid_data$end),4), y = c(0.25, 0.5, 0.75, 1.0), label = c("0.25", "0.50", "0.75", "1.00") , color="grey", angle=15, fontface="bold",vjust=2) +
            geom_text(aes(x=id, y=value, label=as.character(h2)), color="black", fontface="bold",size=3.5,alpha=0.8, hjust=-0.2, vjust=0.5, inherit.aes = FALSE )
            # Add base line information
    if (display_groups == T) {
        p <- p +
            geom_segment(data=base_data, aes(x = -end, y = 1.3, xend = -start, yend = 1.3), colour = "black", alpha=0.9, size=0.6 , inherit.aes = FALSE)  +
            geom_text(data=base_data, aes(x = -title, y = 1.35, label=group), hjust=0, colour = "black", alpha=0.8, fontface="bold", inherit.aes = FALSE)
    }
    p <- p +
            coord_flip(ylim=c(0.0,1.5),xlim=c(-1*(max(grid_data$end)+10),0)) +
            labs(title=POP_Name) +
            theme_minimal() +
            theme(
                    plot.title = element_text(face="bold",size=16,hjust=0.5),
                    legend.title = element_blank(),
                    legend.text = element_text(face="bold",size=12),
                    legend.position = "bottom",
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.margin = unit(c(0,.25,.25,0),"in"),
                    text = element_text(size=11),
                    legend.box.spacing = unit(-0.36,"in")) +
                force_panelsizes(rows = unit(7.5, "in"),
                                 cols = unit(5, "in"))

    ggsave(filename=paste0(workflow,"/traits/plots/h2_barplot.cartesian.",subtitle,".png"), plot=p, device="png", bg="white", width=5, height=9, units="in", dpi=300)
    
    p = ggplot(h2.tb,aes(x=model,y=h2,fill=model,label=sprintf("%0.2f", round(h2, digits = 2)))) +
    		geom_bar(stat="identity", color='black', alpha=0.9, show.legend=FALSE) +
            geom_text(size=5,vjust=-1) +
            scale_fill_manual(values=palette) +
            ylab("Heritability") +
            xlab("Year") +
    		theme_minimal() +
            theme(axis.title = element_text(face="bold",size=16,hjust=0.5),
                  axis.text.x  = element_text(size=14,angle=60),
                  axis.text.y  = element_text(size=14),
                  strip.text = element_text(face="bold",size=18)) +
    		facet_grid(. ~ label_short)
    
    ggsave(filename=paste0(workflow,"/traits/plots/h2_barplot.simple.",subtitle,".png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)
}

h2.tb <- read_csv(h2.collated.path)
upright_yield_traits = c("UBL","UBW","UBM","UTBM","UNS","UMFM","ULvW")
generateH2Plot(h2.tb, upright_yield_traits, paste0("upright_yield.",POP_Name),display_groups=show_groups,display_legend=show_legend)
bbi_traits = c("BBITBM","BBITY","BBISFY","URB")
generateH2Plot(h2.tb, bbi_traits, paste0("bbi.",POP_Name),display_groups=show_groups,display_legend=show_legend)
upright_yield2_traits = c("UL","USL","UDM","URB")
generateH2Plot(h2.tb, upright_yield2_traits, paste0("upright_yield2.",POP_Name),display_groups=show_groups,display_legend=show_legend)
upright_ped_traits = c("UNP","UNB","UNAB","UN0","UNAF")
generateH2Plot(h2.tb, upright_ped_traits, paste0("upright_ped.",POP_Name),display_groups=show_groups,display_legend=show_legend)
upright_other_traits = c("UCD","UCLP","UCLS","UBBL")
generateH2Plot(h2.tb, upright_other_traits, paste0("upright_other.",POP_Name),display_groups=show_groups,display_legend=show_legend)
upright_chimera_shape_traits = c("UKLvW","UKEC","UKSO","UKTO", "UKUX", "UKUY")
generateH2Plot(h2.tb, upright_chimera_shape_traits, paste0("upright_shape.",POP_Name),display_groups=show_groups,display_legend=show_legend)
plot_traits = c("TY","SFY","MFM","PFR","Tacy","Brix","TA","PAC")
generateH2Plot(h2.tb, plot_traits, paste0("plot_traits.",POP_Name),display_groups=show_groups,display_legend=show_legend)


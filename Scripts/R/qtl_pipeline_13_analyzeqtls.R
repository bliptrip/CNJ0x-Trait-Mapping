#!/usr/bin/env RScript

# loading libraries
source('./usefulFunctions.R')

workflow <- "../../Workflows/1"

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

library(Hmisc) #For rcorr -- does pairwise correlations with p-value stats included
library(RColorBrewer)
library(corrplot)
library(GGally)
library(ggthemes)
library(lattice)
library(plotly)
library(tidyverse)

decorate_h2_label <- function(model, h2) {
    model_h2            <- vector("character", length(model))
    nul.tf              <- (h2 == 0)
    model_h2[nul.tf]   <- model[nul.tf]
    model_h2[!nul.tf]  <- paste0(model[!nul.tf],": ",round(h2[!nul.tf],3))
    return(model_h2)
}

circostraits.cfg="model-traits.cfg.csv"
circostraits.cfg.df <- read.csv(file=paste0(workflow,'/configs/',circostraits.cfg),header=T) %>% filter(is.na(mask) | (mask != "TRUE"))

#Just summarize to traits and their labels, short labels 
circostraits.labels.df <- circostraits.cfg.df %>% group_by(trait) %>% summarize(label=label[1],label_short=label_short[1])
#Just summarize to models and their associated labels
circostraits.models.df <- circostraits.cfg.df %>% group_by(model) %>% summarize(model_label=model_label[1])

#filter out models of all-years
qtl.collated.df   <-  read.csv(file=paste0(workflow,"/traits/qtl_collated.consensus.csv"),header=T) %>%
                        left_join(circostraits.labels.df, by="trait") %>%
                        left_join(circostraits.models.df, by="model")

#qtl.ay.df         <- qtl.collated.df %>%
#                        filter((model=="all-years") & (method=="stepwiseqtl")) %>%
#                        group_by(label_short) %>%
#                        select(label_short, chr, marker.variance) %>%
#                        #Pad all chromosomes with 0's to make graph show all chromosomes
#                        rbind(list(label_short=unique(.$label_short),chr=rep(1,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(2,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(3,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(4,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(5,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(6,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(7,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(8,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(9,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(10,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(11,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short))))) %>%
#                        rbind(list(label_short=unique(.$label_short),chr=rep(12,length(unique(.$label_short))),marker.variance=rep(0,length(unique(.$label_short)))))
#

qtl.ay.summary.df <-  qtl.collated.df %>%
                        filter((model=="all-years") & (method=="stepwiseqtl")) %>%
                        group_by(label_short, chr) %>%
                        summarize(chr_qtl_var=sum(marker.variance)) %>%
                        mutate(trait_chr=paste0(label_short, chr))



#Derived from https://www.r-graph-gallery.com/297-circular-barplot-with-groups/
#Generate the data frame
data=data.frame(
            individual=qtl.ay.summary.df$trait_chr,
            group=factor(qtl.ay.summary.df$label_short),
            value=qtl.ay.summary.df$chr_qtl_var,
            chr=qtl.ay.summary.df$chr
        )

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
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
            geom_bar(stat="identity", alpha=0.5) +
            # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
            geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            # Add text showing the value of each 100/75/50/25 lines
            annotate("text", x = rep(max(data$id),4), y = c(10,20,30,40), label = c("10%", "20%", "30%", "40%") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
            ylim(-20,50) +
            coord_polar() + 
            theme_minimal() +
            #theme_solarized(light=F) +
            theme(
                    legend.position = "none",
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.margin = unit(rep(-1,4), "cm") 
                ) +
            geom_text(data=label_data, aes(x=id, y=value+0.1, label=chr, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
            # Add base line information
            geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), colour = "black", alpha=0.9, size=0.6 , inherit.aes = FALSE )  +
            geom_text(data=base_data, aes(x = title, y = -2, label=group), hjust=c(1,1,1,1,0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p

ggsave(filename=paste0(workflow,"/traits/plots/all-years_trait_marker_variance_barplot.png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

#Now try to summarize and display variance of stepwise-qtl on each year separately, and put in a bar plot.
generateTotalLGVariancePlots <- function(qtl.df, qtl.method, traits) {
    qtl.ey.sw.summary.df <-    qtl.df %>%
                                    filter(method==qtl.method) %>%
                                    filter(trait %in% traits) %>%
                                    group_by(label_short, model_label, chr) %>%
                                    summarize(chr_qtl_var=sum(marker.variance))
    qtl.ey.sw.summary.df$chr <-    as.factor(qtl.ey.sw.summary.df$chr)

    o = ggplot(qtl.ey.sw.summary.df) +
            geom_bar(aes(x=as.factor(chr),y=chr_qtl_var,fill=chr), stat="identity", color='black', alpha=1.0, show.legend=FALSE) +
            scale_fill_manual(values=brewer.pal(12,"Paired")) +
            ylab("Total QTL Variance") +
            xlab("Chromosome") +
            theme_minimal() +
            theme(axis.title = element_text(face="bold",size=16),
                axis.text.x  = element_text(size=14,angle=60),
                axis.text.y  = element_text(size=14),
                strip.text = element_text(face="bold",size=18)) +
            facet_grid(model_label ~ label_short)

    ggsave(filename=paste0(workflow,"/traits/plots/each-year_trait_marker_variance_barplot.",qtl.method,".png"), plot=o, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

    nchr = nlevels(qtl.ey.sw.summary.df$chr)
    p = ggplot(qtl.ey.sw.summary.df) +
            geom_bar(aes(x=as.factor(chr),y=chr_qtl_var,fill=chr), stat="identity", color='black', alpha=1.0, show.legend=FALSE) +
            scale_fill_manual(values=brewer.pal(12,"Paired")) +
            coord_polar(start=(pi/nchr),direction=1) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                axis.text.x  = element_text(size=14,angle=60),
                axis.text.y  = element_text(size=14),
                strip.text = element_text(face="bold",size=18)) +
            facet_grid(model_label ~ label_short)

    ggsave(filename=paste0(workflow,"/traits/plots/each-year_trait_marker_variance_barplot_polar.",qtl.method,".png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)
}

yield.traits <- c("berry_length","berry_width","berry_weight","total_berry_weight")
generateTotalLGVariancePlots(qtl.collated.df, "scanone", yield.traits)
generateTotalLGVariancePlots(qtl.collated.df, "stepwiseqtl", yield.traits)

save.image(".RData.13_analyze")

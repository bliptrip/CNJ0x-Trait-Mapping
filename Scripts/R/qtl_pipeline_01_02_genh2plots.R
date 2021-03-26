#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script generates the correlation plots derived from means of two Vorsa populations along with the correlations of 
# both populations combined..

# loading libraries
source('./usefulFunctions.R')

workflow <- get0("workflow", ifnotfound="../../Workflows/1")

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

collateH2sCB <- function(trait.cfg, trait.path, h2.collated.df.p) {
    h2.trait.df <- read.csv(file=paste0(trait.path,"/h2.csv"))
    append.pointer(h2.collated.df.p, c(trait.cfg$model, trait.cfg$trait, h2.trait.df$Estimate, h2.trait.$SE))
}

h2.df            <- data.frame(model=character(),trait=character(), h2=numeric(), h2_se=numeric(), stringsAsFactors=FALSE)
h2.collated.df.p <- newPointer(h2.df)
loopThruTraits(workflow, collateH2sCB, h2.collated.df.p)
h2.df            <- h2.collated.df.p$value
h2.collated.path <- paste0(workflow, "/traits/h2.csv")
write.csv(h2.df, file=h2.collated.path)


decorate_h2_label <- function(model, h2) {
    model_h2            <- vector("character", length(model))
    nul.tf              <- (h2 == 0)
    model_h2[nul.tf]   <- model[nul.tf]
    model_h2[!nul.tf]  <- paste0(model[!nul.tf],": ",round(h2[!nul.tf],3))
    return(model_h2)
}

h2.df <- read.csv(h2.collated.path, header=T, stringsAsFactors=F)
h2.df$model <- as.factor(h2.df$model) #Cast to factor to have plots correctly treat as discrete
h2.df$trait.abbrev <- rename_traits(h2.df$trait, trait.abbrev.map.df)
h2.df$trait.abbrev <- factor(h2.df$trait.abbrev, levels=c('BL','BW','BM','TBM','NB','NS'))
h2.nul.idx <- which(h2.df$h2 == 0)
h2.df$se[h2.nul.idx] <- NA
h2.df <- h2.df %>%
            mutate(model_trait=paste0(model,"-",trait.abbrev)) %>%
            mutate(model_h2=decorate_h2_label(model,h2))


#Derived from https://www.r-graph-gallery.com/297-circular-barplot-with-groups/
#Generate the data frame
data=data.frame(
        individual=h2.df$model_trait,
        group=h2.df$trait.abbrev,
        value=h2.df$h2,
        model=as.character(h2.df$model),
        model_h2=as.character(h2.df$model_h2),
        se=h2.df$se
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
p = ggplot(data, aes(x=as.factor(id), y=value, fill=as.factor(model))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        geom_bar(stat="identity", color="black", alpha=0.9) +
        scale_fill_manual(values=brewer.pal(3,"Accent")) +
        #geom_errorbar(aes(x=as.factor(id), ymin=value-(se/2), ymax=value+(se/2))) +
        #geom_pointrange(aes(x=as.factor(id), y=value, ymin=value-(se/2), ymax=value+(se/2)), na.rm=TRUE) +
            # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
            geom_segment(data=grid_data, aes(x = end, y = 1.0, xend = start, yend = 1.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
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
            geom_text(data=base_data, aes(x = title, y = -.2, label=group), hjust=c(0.5,0.5,0.5,0.35,0.35,0.35), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(filename=paste0(workflow,"/traits/plots/h2_polar_barplot.png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

#Make a non-polar version of plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
        #geom_errorbar(aes(x=as.factor(id), ymin=value-(se/2), ymax=value+(se/2))) +
        geom_pointrange(aes(x=as.factor(id), y=value, ymin=value-(se/2), ymax=value+(se/2)), na.rm=TRUE) +
            # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
            geom_segment(data=grid_data, aes(x = end, y = 1.0, xend = start, yend = 1.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
            # Add text showing the value of each 100/75/50/25 lines
            annotate("text", x = rep(max(data$id),4), y = c(0.25, 0.5, 0.75, 1.0), label = c("0.25", "0.50", "0.75", "1.00") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
            geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
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
            geom_text(data=label_data, aes(x=id, y=value+0.1, label=model_h2, hjust=0), color="black", fontface="bold",alpha=0.8, size=4, angle=0, inherit.aes = FALSE ) +
            # Add base line information
            geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), colour = "black", alpha=0.9, size=0.6 , inherit.aes = FALSE )  +
            geom_text(data=base_data, aes(x = title, y = -.2, label=group), hjust=c(0.5,0.5,0.5,0.35,0.35,0.35), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(filename=paste0(workflow,"/traits/plots/h2_barplot.png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

p = ggplot(h2.df,aes(x=model,y=h2,fill=model,label=sprintf("%0.2f", round(h2, digits = 2)))) +
		geom_bar(stat="identity", color='black', alpha=0.9, show.legend=FALSE) +
        geom_text(size=5,vjust=-1) +
        scale_fill_manual(values=brewer.pal(3,"Accent")) +
        ylab("Heritability") +
        xlab("Year") +
		theme_minimal() +
        theme(axis.title = element_text(face="bold",size=16),
              axis.text.x  = element_text(size=14,angle=60),
              axis.text.y  = element_text(size=14),
              strip.text = element_text(face="bold",size=18)) +
		facet_grid(. ~ trait.abbrev)

#Add for error bars on the heritabilities
#geom_errorbar(aes(x=model,ymin=h2-(se/2),ymax=h2+(se/2))) +

ggsave(filename=paste0(workflow,"/traits/plots/h2_barplot_simple.png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

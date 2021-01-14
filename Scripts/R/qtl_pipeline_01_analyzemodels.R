#year_var=numeric(), This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#!/usr/bin/env RScript

#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular run permutations script will be run on CHTC cluster - thus, it needs to be invoked as an RScripts


# loading libraries
source('./usefulFunctions.R')

workflow <- "../../Workflows/5"

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

library(ggplot2)
library(sommer)
library(tidyverse)

#Loop through previously generated models and gather summary stats on them, plot residuals, etc.  Consider comparing across years and also for all-years.
#Things to assess/plot:
# 1) Model significances.
# 2) Total model variances for BLUPs, residuals, and other covariates.
# 3) Plot of residuals versus index to see if residuals are homo- or heteroscedastic
# 4) Compare residuals across three years.
# 5) Look at year-to-year correlation of BLUPs.  Are they tight for 1st-order differences?  What about 2nd-order differences?
# 6) How about spearman-rank correlations.  What is the makeup of this compared to Pearson's correlation (should be similar).
# 7) Are the years significantly different for the three years?  What about pairwise comparisons of years?


readModelsCB  <- function(trait.cfg, trait.path, models.l) {
    model.map.l         <- models.l$model_map 
    model.collated.df   <- models.l$model_collated_table
    model               <- readRDS(file=paste0(trait.path,"/mmer.rds"))
    append.pointer(model.map.l, model)
    model_idx           <- length(model.map.l$value)
    year_var            <- NA
    if( "year" %in% names(model$var.comp) ) {
        year_var        <- model$var.comp$year
    }
    append.pointer(model.collated.df, c(trait.cfg$model, trait.cfg$trait, model_idx, model$var.comp$accession_name, year_var, model$var.comp$units))
    #Interesting components of model
    # $var.comp[$year, accession_name, ...]
    # $u.hat - Estimated BLUPs for RE covariates (years & genotypic values)
    # $Var.u.hat - Variance estimates for RE covariates (years & genotypic values)
    # $fitted.y
    # $fitted.u
    # $LL
    # $AIC
    # $BIC
    # $ZETA - the original Z-matrices specified for RE (year & genotypic values)
    #AMNOTE: Why is genomic relationship matrix of individual w/ itself not 1?  It looks like this might be intrinsic to the way additive relationship matrices
    # are calculated (XX'/c).
}

model.map.l.p       <- newPointer(list())
model.collated.df   <- data.frame(model=character(),traits=character(),model_idx=integer(), bv_var=numeric(), year_var=numeric(), res_var=numeric(), stringsAsFactors=FALSE)
model.collated.df.p <- newPointer(model.collated.df)
loopThruTraits(workflow, readModelsCB, list(model_map=model.map.l.p, model_collated_table=model.collated.df.p))

#Put the collated model data into a long-format for use by ggplot2
model.spread.df <- (model.collated.df.p$value) %>%
                        gather(type, var, -c(model, mtraits, model_idx))

model.abbrev.df <- data.frame(rows=c("2011","2012","2013","all-years"), abbrevs=c("11","12","13","AY"), stringsAsFactors=F, row.names="rows")
model.spread.df$model <- rename_traits(model.spread.df$model, model.abbrev.df)
model.spread.df$mtraits <- rename_traits(model.spread.df$mtraits, trait.abbrev.map.df)
                    
model.spread.df <- model.spread.df %>%
                    mutate(mtraits_model=paste0(mtraits,"-",model)) %>%
                    group_by(mtraits_model) %>%
                    arrange(type, .by_group=TRUE)

generate_h2 <- function(type, var) {
    i <- which(type == "bv_var")
    var = as.numeric(var)
    return((var[i])/sum(var, na.rm=TRUE))
}

model.h2.df <- model.spread.df %>% summarize(h2=generate_h2(type, var))
                    
                    

               

g <- ggplot(na.omit(model.spread.df), aes(x=as.factor(mtraits_model), y=as.numeric(var), fill=type)) +
        geom_bar(stat="identity")



#sort/group by trait and compare genotype variance, covariate variance, and total variance.  Consider doing a stacked bar plot of the different covariance components

#Derived from https://www.r-graph-gallery.com/297-circular-barplot-with-groups/
#Generate the data frame
model.spread.nona.df <- na.omit(model.spread.df)
data=data.frame(
        individual=model.spread.nona.df$mtraits_model,
        group=model.spread.nona.df$mtraits,
        value=as.numeric(model.spread.nona.df$var),
        type=model.spread.nona.df$type
        )

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar)
to_add$individual=paste0(rep(paste0(levels(data$group),"-NA"), each=empty_bar),as.character(seq(1,empty_bar)))
data=rbind(data, to_add)
data=data %>% arrange(group)
data$id=seq(1, nrow(data))

data$individual <- as.factor(as.character(data$individual))
idx <- 0
collapse_ids <- function(x) {
    .GlobalEnv$idx <- .GlobalEnv$idx + 1;
    return(rep(.GlobalEnv$idx, length(x)));
}
data <- data %>% group_by(individual) %>%
         mutate(id = collapse_ids(individual)) 

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
p = ggplot(data, aes(x=as.factor(id), y=value, fill=type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        geom_bar(aes(x=as.factor(id), y=value, fill=type), stat="identity", alpha=0.5) +
        # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
        geom_segment(data=grid_data, aes(x = end, y = 30.0, xend = start, yend = 30.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 22.5, xend = start, yend = 22.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 15.0, xend = start, yend = 15.0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 7.5, xend = start, yend = 7.5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        # Add text showing the value of each 100/75/50/25 lines
        annotate("text", x = rep(max(data$id),4), y = c(7.5,15.0,22.5,30.0), label = c("7.50","15.0","22.5","30.0") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
        ylim(-5,27.0) +
        theme_minimal() +
        theme(
                legend.position = "none",
                axis.text = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                plot.margin = unit(rep(-1,4), "cm") 
            ) +
        coord_polar() + 
        #geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
        # Add base line information
        geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), colour = "black", alpha=0.9, size=0.6 , inherit.aes = FALSE )  +
        geom_text(data=base_data, aes(x = title, y = -.2, label=group), hjust=c(0.5,0.5,0.5,0.35,0.35,0.35), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(filename=paste0(workflow,"/traits/plots/h2_polar_barplot.png"), plot=p, device="png", bg="white", width=33, height=25, units="cm", dpi=300)

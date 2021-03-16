#year_var=numeric(), This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#!/usr/bin/env RScript

#Analyzes previously fitted models by generating anova tables using a drop1 term likelihood ratio approach to determine significance of an effect,
#and based on whether g:y interaction effects are significant or not (where applicable), chooses a better model and generates relevant cross
#files for r/qtl using selected model BLUPs.
#


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

library(ggplot2)
library(qtl)
library(rlist)
library(sommer)
library(tidyverse)

source(paste0(workflow,'/configs/model.cfg'))

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
    trait               <- trait.cfg$trait
    model.map.l         <- models.l$model_map 
    model.collated.df   <- models.l$model_collated_table
    model               <- readRDS(file=paste0(trait.path,"/mmer.rds"))
    print(paste0("Calculating pvalues for Trait: ",trait.cfg$trait,"Model: ",trait.cfg$model))
    append.pointer(model.map.l, model)
    model_idx           <- length(model.map.l$value)
    var.comp <- summary(model)$varcomp
    var.i <- grep("u:id", rownames(var.comp)) #Return variance component with genotype identifier
    units.i <- grep("units", rownames(var.comp))
    append.pointer(model.collated.df, c(trait.cfg$model, trait.cfg$trait, model_idx, var.comp$VarComp[var.i], var.comp$VarComp[units.i]))
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
    
    
    #Use a 'drop1' type approach to determine significance of model terms
    #Get our A matrix, since it will be necessary.
    A <- model$A.mat
    #Gather fixed terms
    fixed <- as.formula(trait.cfg$fixed)
    yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
    response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    fixedtermss <- apply(data.frame(yuyuf),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fixedtermss <- fixedtermss[which(fixedtermss != "-1")]
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
    #Gather random terms
    random <- as.formula(trait.cfg$random)
    yuyur <- strsplit(as.character(random[2]), split = "[+-]")[[1]]
    randomtermss <- apply(data.frame(yuyur),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    randomtermss    <- randomtermss[which(randomtermss != "-1")]
    randomtermss    <- randomtermss[which(randomtermss != "1")]
    rcov <- as.formula(trait.cfg$rcov)
    yuyuu <- strsplit(as.character(rcov[2]), split = "[+-]")[[1]]
    rcovtermss <- apply(data.frame(yuyuu),1,function(x) {
        strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    allterms          <- c(fixedtermss,randomtermss,rcovtermss)
    allterms.size     <- length(allterms)
    anova.combined    <- tibble(dropterm=character(allterms.size), #Term that is dropped from the model
								Chisq=character(allterms.size),
								ChiDf=character(allterms.size),
								PrChisq=numeric(allterms.size),
								PrChisqInfo=character(allterms.size))
	anova.models	  <- list('full'=model)
    #rownames(anova.combined) <- as.character(allterms)
    #Now for each fixed term, drop it from model and calculate prob under assumption of chisqr distro for likelihood ratio
    i = 1
    for (j in seq_along(fixedtermss)) {
        usef <- setdiff(fixedtermss,fixedtermss[i])
        if( length(usef) > 0 ) {
            fixedf <- paste(response,"~ ",paste(usef,collapse = " + "))
        } else {
            fixedf <- paste(response,"~ 1")
        }
        mmer.expr <- paste0(c("mmer(fixed=", fixedf, ", random=", trait.cfg$random, ", rcov=", trait.cfg$rcov, ", data=model$dataOriginal"), collapse="")
        if( !is.empty(trait.cfg["mmer_args"]) ) {
                mmer.expr <- paste0(mmer.expr,",",trait.cfg["mmer_args"])
        }
        mmer.expr <- paste0(mmer.expr,")")
        print(mmer.expr)
        model.drop1 <- try(eval(parse(text=mmer.expr)))
        if( length(model.drop1) > 1 ) {
            anova.cols    <- c("Chisq", "ChiDf", "PrChisq", "PrChisqInfo")
            anova.drop1 <- anova(model, model.drop1)[2,]
            prchisq <- strsplit(anova.drop1$PrChisq, split = " ")[[1]]
            anova.drop1$PrChisq <- as.numeric(prchisq[1])
            anova.drop1$PrChisqInfo <- prchisq[2]
            anova.combined[i,anova.cols] = anova.drop1[,anova.cols]
			anova.combined[i,"dropterm"] = fixedtermss[i]
			anova.models <- eval(parse(text=paste0("list.append(anova.models,'",fixedtermss[i],"'=model.drop1)")))
            i = i + 1
        }
    }    
    #Now for each random term, drop it from model and calculate prob under assumption of chisqr distro for likelihood ratio
    k = 1
    for (j in seq_along(randomtermss)) {
        usef <- setdiff(randomtermss,randomtermss[k])
        randomf <- paste("~",paste(usef,collapse = " + "))
        mmer.expr <- paste0(c("mmer(fixed=", trait.cfg$fixed, ", random=", randomf, ", rcov=", trait.cfg$rcov, ", data=model$dataOriginal"), collapse="")
        if( !is.empty(trait.cfg["mmer_args"]) ) {
                mmer.expr <- paste0(mmer.expr,",",trait.cfg["mmer_args"])
        }
        mmer.expr <- paste0(mmer.expr,")")
        print(mmer.expr)
        model.drop1 <- try(eval(parse(text=mmer.expr)))
        if( length(model.drop1) > 1 ) {
            anova.cols    <- c("Chisq", "ChiDf", "PrChisq", "PrChisqInfo")
            anova.drop1 <- anova(model, model.drop1)[2,]
            prchisq <- strsplit(anova.drop1$PrChisq, split = " ")[[1]]
            anova.drop1$PrChisq <- as.numeric(prchisq[1])
            anova.drop1$PrChisqInfo <- prchisq[2]
            anova.combined[k+(i-1),anova.cols] = anova.drop1[,anova.cols]
			anova.combined[k+(i-1),"dropterm"] = randomtermss[k]
			anova.models <- eval(parse(text=paste0("list.append(anova.models,'",randomtermss[k],"'=model.drop1)")))
            k = k + 1
        }
    }    
    vcov <- read.csv(file=paste0(trait.path,'/vcov.csv'), row.names=1)
    vcov$PrNorm <- pnorm(abs(vcov$Zratio),lower.tail=FALSE)
    vcov$PrNormInfo <- unlist(map(vcov$PrNorm,siginfo))
	vcov$dropterm <- rownames(vcov)
	vcov <- as_tibble(vcov)
    vcov.cols <- c("dropterm","Zratio","PrNorm","PrNormInfo")
	anova.combined <- anova.combined %>%
						filter(dropterm != "") %>%
						mutate(PrChisqInfo = ifelse(is.na(PrChisqInfo),"NS",PrChisqInfo)) %>%
						left_join(vcov %>% select(vcov.cols), by="dropterm")
    #Where applicable, add Zratio values for random effects (as found in vcov matrix) and pnorm values for these
    saveRDS(anova.combined, file=paste0(trait.path,"/anova.rds"), compress=TRUE)
	saveRDS(anova.models, file=paste0(trait.path,"/anova.models.drop1.rds"), compress=TRUE)
    write.csv(anova.combined, file=paste0(trait.path,"/anova.csv"), row.names=TRUE)
}

model.map.l.p       <- newPointer(list())
model.collated.df   <- data.frame(model=character(),traits=character(),model_idx=integer(), bv_var=numeric(), res_var=numeric(), stringsAsFactors=FALSE)
model.collated.df.p <- newPointer(model.collated.df)
loopThruTraits(workflow, readModelsCB, list(model_map=model.map.l.p, model_collated_table=model.collated.df.p))

#Now that we've generated anova tables, choose the optimal model for all-years based on significance values of any g:y effects -- if not significant, then choose
#the model without this effect, as I've found that sommer's mmer() function will give invalid g:y variance values (negative values) under some situations, and if I remove
#this from the modeling term, the variance terms look more reasonable
#Use the selected model to generate relevant cross file for QTL analysis and BLUP plot generation/etc.
selectModelCB <- function(trait.cfg, trait.path, selectedModels) {
	trait		  <- trait.cfg$trait
	model_label	  <- trait.cfg$model
	anova.drop1.selectedmodel.label <- NA
	if( model_label == 'all-years' ) { #Interaction effects only apply for 'all-years' model - year-by-year model doesn't have an interaction effect
		anova.file <- paste0(trait.path,"/anova.csv")
		if( file.exists(anova.file) ) {
			anova.int	<- read_csv(anova.file) %>% filter(dropterm == "id:year") #Filter out row with g:y interaction effects
			if( anova.int$PrChisq > int_effects_alpha ) { #g:y interaction effect is not significant -- choose model that drops this term
				anova.drop1.selectedmodel.label <- "id:year"
			} else {
				anova.drop1.selectedmodel.label <- "full"
			}
		}
	} else {
		anova.drop1.selectedmodel.label <- "full" #Always assess the full model when looking within a given year
	}
	append.pointer(selectedModels, c(model_label,trait,anova.drop1.selectedmodel.label))
}
models.selected.df   <- data.frame(model=character(),trait=character(),selected_model=character())
models.selected.df.p <- newPointer(models.selected.df)
loopThruTraits(workflow, selectModelCB, models.selected.df.p)
models.selected.df <- models.selected.df.p$value
write.csv(models.selected.df, file=paste0(workflow,"/traits/selectedModels.csv"), row.names=F)

generateCrossFilesCB <- function(trait.cfg,trait.path,args.l) {
	gData			 <- args.l[[1]]
	map.df			 <- args.l[[2]]
	selectedModels	 <- args.l[[3]]
	trait_name		 <- trait.cfg$trait
	model_label		 <- trait.cfg$model
	anova.drop1.selectedmodel.file <- paste0(trait.path,"/anova.models.drop1.rds")
	selectedmodel.label <- (selectedModels$value %>% filter((trait == trait_name) & (model == model_label)) %>% select(selected_model))[[1]]
	if( file.exists(anova.drop1.selectedmodel.file) ) {
		anova.drop1.selectedmodel.l <- readRDS(anova.drop1.selectedmodel.file)
		model						<- anova.drop1.selectedmodel.l[[selectedmodel.label[[1]]]]
		blups						<- model$U[["u:id"]][[trait_name]]
		generate_cross_file(trait.cfg, blups, gData, map.df, trait.path)
	}
}

gData.rqtl	      <- readRDS(file=paste0(workflow,"/traits/rqtl.gdata.rds"))
superMap.df       <- read.csv(geno_rpath2fpath(geno_consensus_file),header=T)[,c("marker","LG","consensus")] %>%
                        mutate(binID = generate_bin_id(LG,consensus))
rownames(superMap.df) <- superMap.df$marker
loopThruTraits(workflow,generateCrossFilesCB,list(gData.rqtl,superMap.df,models.selected.df.p))

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

save.image(".RData.01_analyzemodels")

# This file contains code to analyze and remove statistical outliers in the Vorsa upright datasets

#Only install the following packages when running for the first time
install.packages(c("car"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("pryr"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("openxlsx","RColorBrewer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("lm"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

source('./usefulFunctions.R')

#File for calculating basic statistics from the NJ0[24] population datasets.

require(openxlsx)
#Companion to applied regression - contains outlierTest() for regressions.
library(car)
library(lme4)

#If using RScript command-line, can pass filenames as arguments on command-line

#Change these variables to the appropriate names 
dataset.filename  <- "../../Data/phenotypic data/DerivedData/Data-combined-collated.xlsx"
dataset.sheet     <- "Combined CNJ Upright Data"

#CNJ0[24] Dataframe
#NOTE: If you have problems reading the excel file, open in excel and resave.
cnjpop.df.raw  <- read.xlsx(dataset.filename, sheet=dataset.sheet, colNames=FALSE, startRow=3)
colnames(cnjpop.df.raw) <- c("population", "year", "accession", "accession_name", "upright", "num_peds", "num_no_fruit", "num_berries", "num_aborted_flower", "num_aborted_berries", "total_berry_weight", "upright_length", "secondary_growth", "dry_wt_leaves", "berry_length", "berry_width", "berry_weight", "calyx_diam", "calyx_lobe_form", "calyx_lobe_size", "shape_calyx_end", "shape_stem_end", "berry_skin", "berry_shape", "num_seeds", "rebud", "notes1")
saveRDSw(cnjpop.df.raw, "cnjpop.df.raw.rds", compress=TRUE)
cnjpop.df.raw <- readRDSw("cnjpop.df.raw.rds")

#Cleanup helper
cleanup_numeric_helper <- function(text) {
    if( grepl("(n/a)|(n\\a)|(no tissue)|(^-$)|(_)", text) ) {
        text <- NA
    }
    text <- gsub('*', '', text, fixed=T)
    text <- gsub('^', '', text, fixed=T)
    text <- gsub(',', '.', text)
    text <- gsub("^$", NA, text)
    #Remove duplicate decimal points.  Flag as an issue and have manually cross-checked
    #Convert to vector
    if( grepl(".", text, fixed=T) && !grepl("\\", text, fixed=T) ) {
        strcomps <- strsplit(text, ".",fixed=T)[[1]]
        if(length(strcomps) > 2) {
            strcomps_rest <- paste(strcomps[-(1)], collapse="")
            text <- paste(strcomps[1],strcomps_rest,sep=".")
        }
    }
    if(grepl("\\", text, fixed=T)) {
        #Simply take the means of multiple entries
        subtexts <- strsplit(text, "\\", fixed=T)[[1]]
        return(round(mean(as.numeric(subtexts),3)))
    } else {
        if( !is.na(text) ) {
            text <- round(as.numeric(text),3)
        }
        return(as.numeric(text))
    }
}

cnjpop.cols.include.prefix   <- c("population", "year", "accession", "accession_name", "upright")
cnjpop.cols.include   <- c(cnjpop.cols.include.prefix, "num_peds", "num_no_fruit", "num_berries", "num_aborted_flower", "num_aborted_berries", "total_berry_weight", "upright_length", "secondary_growth", "dry_wt_leaves", "berry_length", "berry_width", "berry_weight", "num_seeds")
cnjpop.df <- cnjpop.df.raw[,cnjpop.cols.include]

#Apply a function to cleanup datasets
#Berry weight
cnjpop.df$berry_weight = vapply(cnjpop.df$berry_weight, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
berry_weight_na.tf = is.na(cnjpop.df$berry_weight)
saveRDSw(berry_weight_na.tf, "berry_weight_na.tf.rds", compress=TRUE)
cnjpop.df.raw[berry_weight_na.tf, c(cnjpop.cols.include.prefix, "berry_weight")]   
#Use asterics in berry weight to guide which datasets to remove for berry_length and berry_width (since these are likely not legit for smashed/rotten berries)
#Berry length
cnjpop.df.berry_length_width.exclude.tf <- grepl("*^", cnjpop.df.raw$berry_weight, fixed=T)
cnjpop.df$berry_length[cnjpop.df.berry_length_width.exclude.tf] <- NA
cnjpop.df$berry_length <- vapply(cnjpop.df$berry_length, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
berry_length_na.tf = is.na(cnjpop.df$berry_length)
saveRDSw(berry_length_na.tf, "berry_length_na.tf.rds", compress=TRUE)
cnjpop.df.raw[berry_length_na.tf, c(cnjpop.cols.include.prefix, "berry_length")]   
#Berry width
cnjpop.df$berry_width[cnjpop.df.berry_length_width.exclude.tf] <- NA
cnjpop.df$berry_width <- vapply(cnjpop.df$berry_width, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
berry_width_na.tf = is.na(cnjpop.df$berry_width)
cnjpop.df.raw[berry_width_na.tf, c(cnjpop.cols.include.prefix, "berry_width")]   
saveRDSw(berry_width_na.tf, "berry_width_na.tf.rds", compress=TRUE)
#Berry length and width
#Combine the logic of which entries to remove as berry length and width should both be there and consistent.
berry_length_width_na.tf = berry_length_na.tf | berry_width_na.tf
saveRDSw(berry_length_width_na.tf, "berry_length_width_na.tf.rds", compress=TRUE)
#Number of seeds in largest berry
cnjpop.df$berry_width[cnjpop.df.berry_length_width.exclude.tf] <- NA
cnjpop.df$num_seeds = vapply(cnjpop.df$num_seeds, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_seeds_na.tf = is.na(cnjpop.df$num_seeds)
saveRDSw(num_seeds_na.tf, "num_seeds_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_seeds_na.tf, c(cnjpop.cols.include.prefix, "num_seeds")]   
#Num pedicels
cnjpop.df$num_peds = vapply(cnjpop.df$num_peds, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_peds_na.tf = is.na(cnjpop.df$num_peds)
saveRDSw(num_peds_na.tf, "num_peds_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_peds_na.tf, c(cnjpop.cols.include.prefix, "num_peds")]   
#Number of pedicels without fruit
cnjpop.df$num_no_fruit = vapply(cnjpop.df$num_no_fruit, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_no_fruit_na.tf = is.na(cnjpop.df$num_no_fruit)
saveRDSw(num_no_fruit_na.tf, "num_no_fruit_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_no_fruit_na.tf, c(cnjpop.cols.include.prefix, "num_no_fruit")]   
#Number of berries
cnjpop.df$num_berries = vapply(cnjpop.df$num_berries, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_berries_na.tf = is.na(cnjpop.df$num_berries)
saveRDSw(num_berries_na.tf, "num_berries_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_berries_na.tf, c(cnjpop.cols.include.prefix, "num_berries")]   
#Number of aborted flowers
cnjpop.df$num_aborted_flower = vapply(cnjpop.df$num_aborted_flower, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_aborted_flower_na.tf = is.na(cnjpop.df$num_aborted_flower)
saveRDSw(num_aborted_flower_na.tf, "num_aborted_flower_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_aborted_flower_na.tf, c(cnjpop.cols.include.prefix, "num_aborted_flower")]   
#Number of aborted berries
cnjpop.df$num_aborted_berries = vapply(cnjpop.df$num_aborted_berries, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
num_aborted_berries_na.tf = is.na(cnjpop.df$num_aborted_berries)
saveRDSw(num_aborted_berries_na.tf, "num_aborted_berries_na.tf.rds", compress=TRUE)
cnjpop.df.raw[num_aborted_berries_na.tf, c(cnjpop.cols.include.prefix, "num_aborted_berries")]   
#Total berry weight/upright
cnjpop.df$total_berry_weight = vapply(cnjpop.df$total_berry_weight, FUN=cleanup_numeric_helper, FUN.VALUE=as.numeric(NA))
total_berry_weight_na.tf = is.na(cnjpop.df$total_berry_weight)
saveRDSw(total_berry_weight_na.tf, "total_berry_weight_na.tf.rds", compress=TRUE)
cnjpop.df.raw[total_berry_weight_na.tf, c(cnjpop.cols.include.prefix, "total_berry_weight")]   

#Given interdependency of num peds, num pedicels without fruit, num berries, num aborted flowers, num aborted berries, remove union of all NA's of these
mvariable_union_na.tf <- (num_peds_na.tf | num_no_fruit_na.tf | num_berries_na.tf | num_aborted_flower_na.tf | num_aborted_berries_na.tf)
mvariables <- c("num_peds","num_no_fruit","num_berries","num_aborted_flower","num_aborted_berries")
cnjpop.df[which(mvariable_union_na.tf),mvariables] <- NA
cnjpop.df[which(mvariable_union_na.tf),mvariables]
mvariable_union_n_na.idx <- which(!mvariable_union_na.tf)

saveRDSw(cnjpop.df, "cnjpop.df.rds", compress=TRUE)
cnjpop.df <- readRDSw("cnjpop.df.rds")

#Do consistency checking on num peds, num pedicels without fruit, num berries, num aborted flowers, num aborted berries, total berry weight
#num_peds = num_no_fruit + num_berries + num_aborted_flower + num_aborted_berries
expected_num_peds = Reduce('+',cnjpop.df[,mvariables[-1]])
mvariable_inconsistent.tf <- (expected_num_peds != cnjpop.df$num_peds)
mvariable_inconsistent.idx <- which(mvariable_inconsistent.tf)
write.csv(file=pheno_dpath2fpath("upright_inconsistencies.csv"), cbind(cnjpop.df[mvariable_inconsistent.idx, c(cnjpop.cols.include.prefix, "num_peds")],expected_num_peds[mvariable_inconsistent.idx], cnjpop.df[mvariable_inconsistent.idx, mvariables[-1]]), row.names=T)
#find which inconsistencies are due to num_aborted_flowers + num_aborted_berries being encoded redundantly in num_no_fruit.
diff_num_peds  <- (expected_num_peds - cnjpop.df$num_peds)
num_aborted    <- Reduce('+',cnjpop.df[,c("num_aborted_flower","num_aborted_berries")])
equals.tf      <- ((diff_num_peds == num_aborted) & (cnjpop.df$num_no_fruit >= num_aborted)) & mvariable_inconsistent.tf #Check that num_no_fruit is more than num_aborted & mask out consistent entries
equals.idx     <- which(equals.tf)
#For thos that the difference in the number of expected pedicels from recorded value is equal to the number of aborted, we assume consistency in miscoding and remove those redundant entries from the total
cnjpop.df$num_no_fruit[equals.idx] <- cnjpop.df$num_no_fruit[equals.idx] - num_aborted[equals.idx]
expected_num_peds = Reduce('+',cnjpop.df[,mvariables[-1]])
write.csv(file=pheno_dpath2fpath("upright_inconsistencies_with_fixed.csv"), cbind(cnjpop.df[mvariable_inconsistent.idx, c(cnjpop.cols.include.prefix, "num_peds")],expected_num_peds[mvariable_inconsistent.idx], cnjpop.df[mvariable_inconsistent.idx, mvariables[-1]]), row.names=T)
#Remove the fixed inconsistencies
mvariable_inconsistent.tf <- (expected_num_peds != cnjpop.df$num_peds)
mvariable_inconsistent.idx <- which(mvariable_inconsistent.tf)
write.csv(file=pheno_dpath2fpath("upright_inconsistencies_without_fixed.csv"), cbind(cnjpop.df[mvariable_inconsistent.idx, c(cnjpop.cols.include.prefix, "num_peds")],expected_num_peds[mvariable_inconsistent.idx], cnjpop.df[mvariable_inconsistent.idx, mvariables[-1]]), row.names=T)
#Now fix entries where the number of aborted flowers is redundantly encoded in the number of pedicels without fruit.
diff_num_peds  <- (expected_num_peds - cnjpop.df$num_peds)
equals.tf      <- ((diff_num_peds == cnjpop.df$num_aborted_flower) & (cnjpop.df$num_no_fruit >= cnjpop.df$num_aborted_flower)) & mvariable_inconsistent.tf #Mask out consistent entries
equals.idx     <- which(equals.tf)
cnjpop.df$num_no_fruit[equals.idx] <- cnjpop.df$num_no_fruit[equals.idx] - cnjpop.df$num_aborted_flower[equals.idx]
expected_num_peds = Reduce('+',cnjpop.df[,mvariables[-1]])
write.csv(file=pheno_dpath2fpath("upright_inconsistencies_with_fixed2.csv"), cbind(cnjpop.df[mvariable_inconsistent.idx, c(cnjpop.cols.include.prefix, "num_peds")],expected_num_peds[mvariable_inconsistent.idx], cnjpop.df[mvariable_inconsistent.idx, mvariables[-1]]), row.names=T)
#Remove the fixed inconsistencies
mvariable_inconsistent.tf <- (expected_num_peds != cnjpop.df$num_peds)
mvariable_inconsistent.idx <- which(mvariable_inconsistent.tf)
write.csv(file=pheno_dpath2fpath("upright_inconsistencies_without_fixed2.csv"), cbind(cnjpop.df[mvariable_inconsistent.idx, c(cnjpop.cols.include.prefix, "num_peds")],expected_num_peds[mvariable_inconsistent.idx], cnjpop.df[mvariable_inconsistent.idx, mvariables[-1]]), row.names=T)
#About 1% are still inconsistent.  Just fill these fields with NAs and ignore.
cnjpop.df[mvariable_inconsistent.idx,mvariables] <- NA
which(cnjpop.df$num_no_fruit[mvariable_inconsistent.idx] < 0)

saveRDSw(mvariable_inconsistent.tf, "mvariable_inconsistent.tf.rds", compress=TRUE)
saveRDSw(cnjpop.df, "cnjpop.df.rds", compress=TRUE)


#Read the state of serialized objects here if don't want to run through the previous commands to cleanup datasets
berry_weight_na.tf <- readRDSw("berry_weight_na.tf.rds")
berry_length_na.tf <- readRDSw("berry_length_na.tf.rds")
berry_width_na.tf <- readRDSw("berry_width_na.tf.rds")
berry_length_width_na.tf <- readRDSw("berry_length_width_na.tf.rds")
num_seeds_na.tf <- readRDSw("num_seeds_na.tf.rds")
mvariable_inconsistent.tf <- readRDSw("mvariable_inconsistent.tf.rds")
cnjpop.df <- readRDSw("cnjpop.df.rds")

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartzw(file="p1_1_phenotype_scatterplots_outliers.pdf", title="Phenotype Scatterplots with Outliers")
pairs(formula=~num_peds+num_no_fruit+num_berries+num_aborted_flower+num_aborted_berries+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with Outliers")
dev.off()

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartzw(file="p1_2_phenotype_num_peds_scatterplots_outliers.pdf", title="Phenotype Scatterplots with Outliers")
pairs(formula=~num_peds+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with Outliers")
dev.off()

quartzw(file="p1_3_phenotype_num_no_fruit_scatterplots_outliers.pdf", title="Phenotype Scatterplots with Outliers")
pairs(formula=~num_no_fruit+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with Outliers")
dev.off()

quartzw(file="p1_4_phenotype_num_berries_scatterplots_outliers.pdf", title="Phenotype Scatterplots with Outliers")
pairs(formula=~num_berries+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with Outliers")
dev.off()

quartzw(file="p1_5_mvariables_scatterplots_outliers.pdf", title="Phenotype Scatterplots with Outliers")
pairs(formula=~num_peds+num_no_fruit+num_berries+num_aborted_flower+num_aborted_berries, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with Outliers")
dev.off()

#Boxplots!

quartzw(file="p2_1_largest_berry_phenotype_boxplots_outliers.pdf", title="Phenotype Boxplots with Outliers")
par(mfrow=c(2,2), las=1, cex=0.8)
cnjpop.berry_weight.year.boxplot       <- boxplot(formula=berry_weight~population, data=cnjpop.df, xlab="Population",ylab="Largest Berry Weight (g)",main="Largest Berry Weight")
cnjpop.berry_length.year.boxplot       <- boxplot(formula=berry_length~population, data=cnjpop.df, xlab="Population",ylab="Largest Berry Length (mm)",main="Largest Berry Length")
cnjpop.berry_width.year.boxplot        <- boxplot(formula=berry_width~population, data=cnjpop.df, xlab="Population",ylab="Largest Berry Width (mm)",main="Largest Berry Width")
cnjpop.berry_num_seeds.year.boxplot    <- boxplot(formula=num_seeds~population, data=cnjpop.df, xlab="Population",ylab="Largest Berry Number of Seeds",main="Largest Berry Number of Seeds")
dev.off()

quartzw(file="p2_2_mvariable_total_berry_weight_phenotype_boxplots_outliers.pdf", title="Phenotype Boxplots with Outliers")
par(mfrow=c(3,2), las=1, cex=0.8)
cnjpop.num_peds.year.boxplot       <- boxplot(formula=num_peds~population, data=cnjpop.df, xlab="Population",ylab="Number of Pedicels",main="Number of Pedicels")
cnjpop.num_no_fruit.year.boxplot       <- boxplot(formula=num_no_fruit~population, data=cnjpop.df, xlab="Population",ylab="Number of Pedicels w/o Fruit",main="Number of Pedicels w/o Fruit")
cnjpop.num_berries.year.boxplot       <- boxplot(formula=num_berries~population, data=cnjpop.df, xlab="Population",ylab="Number of Berries",main="Number of Berries")
cnjpop.num_aborted_flower.year.boxplot       <- boxplot(formula=num_aborted_flower~population, data=cnjpop.df, xlab="Population",ylab="Number of Aborted Flowers",main="Number of Aborted Flowers")
cnjpop.num_aborted_berries.year.boxplot       <- boxplot(formula=num_aborted_berries~population, data=cnjpop.df, xlab="Population",ylab="Number of Aborted Berries",main="Number of Aborted Berries")
cnjpop.total_berry_weight.year.boxplot       <- boxplot(formula=total_berry_weight~population, data=cnjpop.df, xlab="Population",ylab="Total Berry Weight on Upright (g)",main="Total Berry Weight on Upright (g)")
dev.off()

plot_phenotype_hists <- function(dataset, outliers.idx, phenotype, graph.cols) {
    cnjpop.nout.df <- cnjpop.df[-outliers.idx]
    #histbreaks=seq(from=floor(min(cnjpop.nout.df$berry_weight, na.rm=T)), to=ceiling(max(cnjpop.nout.df$berry_weight,na.rm=T)), length.out=21)
    cnjpop.nout.p1.tf <- grepl("CNJ04.*", cnjpop.nout.df$accession_name)
    cnjpop.nout.p1.phenotype <- cnjpop.nout.df[which(cnjpop.nout.p1.tf), phenotype]
    cnjpop.nout.p2.tf <- grepl("CNJ02.*", cnjpop.nout.df$accession_name)
    cnjpop.nout.p2.phenotype <- cnjpop.nout.df[which(cnjpop.nout.p2.tf), phenotype]
    hist(cnjpop.nout.p1.phenotype, col=graph.cols[1], ylim=c(0.0, 1.0), main=paste0("Histogram of ", phenotype, ": Population 1"), xlab="bins", freq=F)
    hist(cnjpop.nout.p2.phenotype, col=graph.cols[2], ylim=c(0.0, 1.0), main=paste0("Histogram of ", phenotype, ": Population 2"), xlab="bins", freq=F)
}

plot_phenotype <- function(lm, phenotype, graph.cols) {
    png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"p3_",phenotype,"_simple_regression.png"), bg="transparent", width=1024, height=1024, units='px')
    oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))
    plot(lm)
    dev.off()
    lm.outliers <- outlierTest(lm)
    lm.outliers.idx <- as.numeric(names(lm.outliers$rstudent))
    png(filename=paste0(DATA_PLOTS_FOLDER_PREFIX,"p3_",phenotype,"_hists.png"), bg="transparent", width=1024, height=1024, units='px')
    par(oma=c(0,0,3,0), mfrow=c(2,1))
    plot_phenotype_hists(cnjpop.df, lm.outliers.idx, phenotype, graph.cols)
    dev.off()
    par(oldpar)
    return(lm.outliers.idx)
}

#Linear regression-based outlier modeling
#Initialize the colors
graph.cols <- colorRampPalette(c(rgb(1/4,1,0,1/4),rgb(0,0,1,1/4)),alpha=T)(2)
#Consider doing regressions to find outliers
berry_weight.lm <- lm(formula=berry_weight~population+year+accession,data=cnjpop.df)
berry_weight.anova <- anova(berry_weight.lm)
berry_weight.lm.outliers.idx <- plot_phenotype(berry_weight.lm, "berry_weight", graph.cols)
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[berry_weight.lm.outliers.idx,]

berry_length.lm <- lm(formula=berry_length~population+year+accession,data=cnjpop.df)
berry_length.anova <- anova(berry_length.lm)
berry_length.lm.outliers.idx <- plot_phenotype(berry_length.lm, "berry_length", graph.cols)
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[berry_length.lm.outliers.idx,]

berry_width.lm <- lm(formula=berry_width~population+year+accession,data=cnjpop.df)
berry_width.anova <- anova(berry_width.lm)
berry_width.lm.outliers.idx <- plot_phenotype(berry_width.lm, "berry_width", graph.cols)
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[berry_width.lm.outliers.idx,]

num_seeds.lm <- lm(formula=num_seeds~population+year+accession,data=cnjpop.df)
num_seeds.anova <- anova(num_seeds.lm)
num_seeds.lm.outliers.idx <- plot_phenotype(num_seeds.lm, "num_seeds", graph.cols)
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[num_seeds.lm.outliers.idx,]

cnjpop.df$upright_length <- as.numeric(cnjpop.df$upright_length)
upright_length.lm <- lm(formula=upright_length~population+year+accession,data=cnjpop.df)
upright_length.anova <- anova(upright_length.lm)
upright_length.lm.outliers.idx <- plot_phenotype(upright_length.lm, "upright_length", graph.cols)
cnjpop.df[upright_length.lm.outliers.idx,]

cnjpop.df$secondary_growth <- as.numeric(cnjpop.df$secondary_growth)
secondary_growth.lm <- lm(formula=secondary_growth~population+year+accession,data=cnjpop.df)
secondary_growth.anova <- anova(secondary_growth.lm)
secondary_growth.lm.outliers.idx <- plot_phenotype(secondary_growth.lm, "secondary_growth", graph.cols)
cnjpop.df[secondary_growth.lm.outliers.idx,]

cnjpop.df$dry_wt_leaves <- as.numeric(cnjpop.df$dry_wt_leaves)
dry_wt_leaves.lm <- lm(formula=dry_wt_leaves~population+year+accession,data=cnjpop.df)
dry_wt_leaves.anova <- anova(dry_wt_leaves.lm)
dry_wt_leaves.lm.outliers.idx <- plot_phenotype(dry_wt_leaves.lm, "dry_wt_leaves", graph.cols)
cnjpop.df[dry_wt_leaves.lm.outliers.idx,]
View(cnjpop.df)

#Now consider throwing out weight outliers (they seem to be wildly off) and num_seeds outliers.
cnjpop.df[c(berry_weight.lm.outliers.idx, num_seeds.lm.outliers.idx),c("berry_weight", "num_seeds")] <- NA

#Now remove outliers based on regressing the variables on each other, given that they seem to be correlated
#<<<<<berry_width v. berry_length>>>>>
berry_width_v_length.lm <- lm(formula=berry_width~berry_length+population+year,data=cnjpop.df)
quartzw(file="p3_5_berry_width_v_length_simple_regression.pdf", "Simple regression plot for berry width versus berry length.")
oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))
par(oldpar)
dev.off()
berry_width_v_length.lm.outliers <- outlierTest(berry_width_v_length.lm)
berry_width_v_length.lm.outliers.idx <- as.numeric(names(berry_width_v_length.lm.outliers$rstudent))
#Throw out outliers
cnjpop.df[berry_width_v_length.lm.outliers.idx,]
cnjpop.df[berry_width_v_length.lm.outliers.idx,c("berry_width","berry_length")] <- NA

#<<<<<berry_weight v. berry_length>>>>>
berry_weight_v_length.lm <- lm(formula=berry_weight~berry_length+population+year,data=cnjpop.df)
quartzw(file="p3_6_berry_weight_v_length_simple_regression.pdf", "Simple regression plot for berry weight versus berry length.")
oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(berry_weight_v_length.lm)
par(oldpar)
dev.off()
berry_weight_v_length.lm.outliers <- outlierTest(berry_weight_v_length.lm)
berry_weight_v_length.lm.outliers.idx <- as.numeric(names(berry_weight_v_length.lm.outliers$rstudent))
#Throw out outliers
cnjpop.df[berry_weight_v_length.lm.outliers.idx,]
cnjpop.df[berry_weight_v_length.lm.outliers.idx,c("berry_weight","berry_length")] <- NA

#<<<<<num_seeds v. berry_weight>>>>>
num_seeds_v_weight.lm <- lm(formula=num_seeds~berry_weight+population+year,data=cnjpop.df)
quartzw(file="p3_7_num_seeds_v_weight_simple_regression.pdf", "Simple regression plot for num seeds versus berry weight.")
oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(num_seeds_v_weight.lm)
par(oldpar)
dev.off()
num_seeds_v_weight.lm.outliers <- outlierTest(num_seeds_v_weight.lm)
num_seeds_v_weight.lm.outliers.idx <- as.numeric(names(num_seeds_v_weight.lm.outliers$rstudent))
#Throw out outliers
cnjpop.df[num_seeds_v_weight.lm.outliers.idx,]
cnjpop.df[num_seeds_v_weight.lm.outliers.idx,c("num_seeds","berry_weight")] <- NA

###############
#Other datasets
###############
#num_peds.lm <- lm(formula=num_peds~population+year+accession+num_berries+num_no_fruit+num_aborted_flower+num_aborted_berries,data=cnjpop.df)
num_peds.lm <- lm(formula=num_peds~population+year+accession,data=cnjpop.df)
num_peds.anova <- anova(num_peds.lm)
quartzw(file="p4_1_num_peds_simple_regression.pdf", "Simple regression plot for number of pedicels.")
num_peds.lm.outliers.idx <- plot_phenotype(num_peds.lm, "num_peds", graph.cols)
dev.off()
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[num_peds.lm.outliers.idx,]
#Don't remove these outliers, as they look to be coded legitimately

num_berries.lm <- lm(formula=num_berries~population+year+accession,data=cnjpop.df)
num_berries.anova <- anova(num_berries.lm)
quartzw(file="p4_2_num_berries_simple_regression.pdf", "Simple regression plot for number of seeds.")
num_berries.lm.outliers.idx <- plot_phenotype(num_berries.lm, "num_berries", graph.cols)
dev.off()
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[num_berries.lm.outliers.idx,]
#Don't remove these outliers, as they look to be coded legitimately

num_peds_all.lm <- lm(formula=num_peds~population+year+accession+num_berries+num_no_fruit+num_aborted_flower+num_aborted_berries,data=cnjpop.df)
num_peds_all.anova <- anova(num_peds_all.lm)
quartzw(file="p4_3_num_peds_all_simple_regression.pdf", "Simple regression plot for number of pedicels modeled over all other related variables.")
num_peds_all.lm.outliers.idx <- plot_phenotype(num_peds_all.lm, "num_peds", graph.cols)
dev.off()
cnjpop.df[num_peds_all.lm.outliers.idx,]
#Don't remove these outliers, as they look to be coded legitimately

total_berry_weight.lm <- lm(formula=total_berry_weight~population+year+accession+num_berries,data=cnjpop.df)
total_berry_weight.anova <- anova(total_berry_weight.lm)
quartzw(file="p5_1_total_berry_weight_simple_regression.pdf", "Simple regression plot for total berry weight.")
total_berry_weight.lm.outliers.idx <- plot_phenotype(total_berry_weight.lm, "total_berry_weight", graph.cols)
dev.off()
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[total_berry_weight.lm.outliers.idx,c(cnjpop.cols.include.prefix,"total_berry_weight")]

cnjpop.df[num_peds.lm.outliers.idx,]
total_berry_weight_v_num_berries.lm <- lm(formula=total_berry_weight~num_berries,data=cnjpop.df)
total_berry_weight_v_num_berries.anova <- anova(total_berry_weight_v_num_berries.lm)
quartzw(file="p5_2_total_berry_weight_v_num_berries_simple_regression.pdf", "Simple regression plot for total berry weight vs. num berries.")
total_berry_weight_v_num_berries.lm.outliers.idx <- plot_phenotype(total_berry_weight_v_num_berries.lm, "total_berry_weight", graph.cols)
dev.off()
#Show which entries are outliers for cross-referencing in excel file.
cnjpop.df[total_berry_weight_v_num_berries.lm.outliers.idx,c(cnjpop.cols.include.prefix,"total_berry_weight")]

#Save the dataset with outliers removed for future access in QTL analysis.
saveRDSw(cnjpop.df, "cnjpop.noout.df.rds", compress=T)
#Read this in here if want to start from here.
cnjpop.df <- readRDSw("cnjpop.noout.df.rds")

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartzw(file="p6_phenotype_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~year+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartzw(file="p7_1_phenotype_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~num_peds+num_no_fruit+num_berries+num_aborted_flower+num_aborted_berries+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartzw(file="p7_2_phenotype_num_peds_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~num_peds+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()

quartzw(file="p7_3_phenotype_num_no_fruit_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~num_no_fruit+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()

quartzw(file="p7_4_phenotype_num_berries_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~num_berries+total_berry_weight+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()

quartzw(file="p7_5_mvariables_scatterplots_NO_outliers.pdf", title="Phenotype Scatterplots with_out_ Outliers")
pairs(formula=~num_peds+num_no_fruit+num_berries+num_aborted_flower+num_aborted_berries, data=cnjpop.df, na.action=na.omit, main="Phenotype Scatterplots with_out_ Outliers")
dev.off()


#Now that outliers have been removed, rerun linear model on all phenotypes separating years/within years, and graph boxplots of the years separated/combined
#to see if residual variation is different separated/combined
years <- unique(cnjpop.df$year)
years.num <- length(years)
berry_weight.lm.l   <- vector("list",years.num+1)
berry_length.lm.l   <- vector("list",years.num+1)
berry_width.lm.l    <- vector("list",years.num+1)
num_seeds.lm.l      <- vector("list",years.num+1)
for( i in 1:years.num ) {
    berry_weight.lm.l[i]   <- lm(formula=berry_weight~population+accession,data=cnjpop.df[which(cnjpop.df$year == years[i]),])
    berry_length.lm.l[i]   <- lm(formula=berry_length~population+accession,data=cnjpop.df[which(cnjpop.df$year == years[i]),])
    berry_width.lm.l[i]    <- lm(formula=berry_width~population+accession,data=cnjpop.df[which(cnjpop.df$year == years[i]),])
    num_seeds.lm.l[i]      <- lm(formula=num_seeds~population+accession,data=cnjpop.df[which(cnjpop.df$year == years[i]),])
}
i <- i+1
#combine all years
berry_weight.lm.l[i]   <- lm(formula=berry_weight~population+year+accession,data=cnjpop.df)
berry_length.lm.l[i]   <- lm(formula=berry_length~population+year+accession,data=cnjpop.df)
berry_width.lm.l[i]    <- lm(formula=berry_width~population+year+accession,data=cnjpop.df)
num_seeds.lm.l[i]      <- lm(formula=num_seeds~population+year+accession,data=cnjpop.df)
#Now loop through the above list and do a boxplot of residuals
for( i in 1:years.num+1 ) {
    berry_weight.lm.l[ <- 
}

#Just quickly assess remaining simple plots
quartzw(file="p11_simple_plot_of_phenotyeps.pdf", title="Phenotype Plot Matrix")
par(mfrow=c(4,2))
plot(cnjpop.df[which(cnjpop.df$population==1),"berry_weight"],ylab="Population 1 Berry Weight")
plot(cnjpop.df[which(cnjpop.df$population==2),"berry_weight"],ylab="Population 2 Berry Weight")
plot(cnjpop.df[which(cnjpop.df$population==1),"berry_width"],ylab="Population 1 Berry Width")
plot(cnjpop.df[which(cnjpop.df$population==2),"berry_width"],ylab="Population 2 Berry Width")
plot(cnjpop.df[which(cnjpop.df$population==1),"berry_length"],ylab="Population 1 Berry Length")
plot(cnjpop.df[which(cnjpop.df$population==2),"berry_length"],ylab="Population 2 Berry Length")
plot(cnjpop.df[which(cnjpop.df$population==1),"num_seeds"],ylab="Population 1 Num Seeds")
plot(cnjpop.df[which(cnjpop.df$population==2),"num_seeds"],ylab="Population 2 Num Seeds")
dev.off()

#Now write to workbook, and add style to elements that are outliers
wb <- createWorkbook()
addWorksheet(wb, dataset.sheet)
writeData(wb, 1, cnjpop.df.raw)
outlier_style <- createStyle(bgFill="red",fgFill="red")
na_style <- createStyle(bgFill="green",fgFill="green")
addStyle(wb, 1, outlier_style, rows=which(berry_weight.outliers), cols=which(colnames(cnjpop.df.raw) == "berry_weight"))
addStyle(wb, 1, na_style, rows=which(berry_weight_na.tf), cols=which(colnames(cnjpop.df.raw) == "berry_weight"))
addStyle(wb, 1, outlier_style, rows=which(berry_length.outliers), cols=which(colnames(cnjpop.df.raw) == "berry_length"))
addStyle(wb, 1, na_style, rows=which(berry_length_na.tf), cols=which(colnames(cnjpop.df.raw) == "berry_length"))
addStyle(wb, 1, outlier_style, rows=which(berry_width.outliers), cols=which(colnames(cnjpop.df.raw) == "berry_width"))
addStyle(wb, 1, na_style, rows=which(berry_width_na.tf), cols=which(colnames(cnjpop.df.raw) == "berry_width"))
addStyle(wb, 1, outlier_style, rows=which(berry_num_seed.outliers), cols=which(colnames(cnjpop.df.raw) == "num_seeds"))
addStyle(wb, 1, na_style, rows=which(num_seeds_na.tf), cols=which(colnames(cnjpop.df.raw) == "num_seeds"))
saveWorkbook(wb, file = paste(dataset.filename,"_outliers.xlsx", sep=""), overwrite = TRUE)


#Observe outliers and write to excel file

#Find the column index that corresponds to the weight
#V_Y <- vc[which(vc$grp == "year"),4]
#V_G <- vc[which(vc$grp == "accession_name"), 4]
#V_P <- vc[which(vc$grp == "population"), 4]
#Ve  <- vc[which(vc$grp == "Residual"),4]
#vc   <- as.data.frame(VarCorr(berry_weight_model_with_outliers))
##blup
#blup <- ranef(berry_weight_model_with_outliers,condVar=T,whichel="accession_name")
#
##Generate a scatterplot of some of the data subsets we care about to visually look for any relationships
#cnjpop.df$berry_length = vapply(cnjpop.df$berry_length, FUN=as.numeric, FUN.VALUE=as.numeric(NA))
#cnjpop.df$berry_width = vapply(cnjpop.df$berry_width, FUN=as.numeric, FUN.VALUE=as.numeric(NA))
#cnjpop.df$num_seeds = vapply(cnjpop.df$num_seeds, FUN=as.numeric, FUN.VALUE=as.numeric(NA))
#spm(~berry_weight+berry_length+berry_width+num_seeds, reg.line=FALSE, smooth=FALSE, spread=FALSE, span=0.5, ellipse=FALSE, levels=c(.5, .9), id.n=0, diagonal = 'density', data=cnjpop.df) 
#pairs()

#Function returns outliers from each grouping as vector of TRUE, FALSE
# 
# Arguments:
#  data - The original dataset frame to operate on.
#  variable - The dataset variable to gather statistics on.
#  subsets - Datafram of T/F vectors that specify what subsets of data to calculate statistics on.  Note: Column names correspond to the subset names that will be used in the plots/graphs.
#           NOTE: Each subset object in subsets should have an attribute 'name' that specifies a descriptor to be used in graph legends.
#  subset_mask - A logical vector of length equal to the number of obs. in data that specifies which entries to include in the data
#  threshold - The standard deviation threshold to consider the cutoff for inliers vs. outliers
#
# Returns:
#  A vector of T/F vectors specifying the outliers in the data$variable input set
assess_population_stats <- function(data, variable, subsets, subset_mask=NULL, threshold=2.57) {
data.len   <- length(data[,variable])
subset.len <- length(subsets)
outlier.tf <-  rep(FALSE, data.len)
inlier.tf  <- rep(FALSE, data.len)
subset.names <- rep("", subset.len)
subset.legends <- rep("", subset.len)
hists   <- list()
#Initialize the colors
graph.cols <- colorRampPalette(c(rgb(1/4,1,0,1/4),rgb(0,0,1,1/4)),alpha=T)(subset.len)

#Create a new graphics window
#Loop through subsets, finding outliers and generating plots
for (ix in 1:subset.len) {
subset <- subsets[[ix]]
if( !is.null(subset_mask) ) {
    subset <- subset & subset_mask
    attr(subset, "name") <- attr(subsets[[ix]],"name")
}
subset.ix <- which(subset)
subset.quantiles <- quantile(data[subset.ix,variable], probs=seq(0,1,1/10))
subset.inliers.tf <- (data[,variable] >= subset.quantiles[1] & data[,variable] <= subset.quantiles[10]) & subset
subset.inliers.ix <- which(subset.inliers.tf)
#Calculate sd and means of each, and then normalize to find outliers
subset.inliers.mean   <- mean(data[subset.inliers.ix,variable])
subset.inliers.sd     <- sd(data[subset.inliers.ix,variable])
subset.sd.normalized <- rep(subset.inliers.mean, data.len)
subset.sd.normalized[subset.ix]  <- abs(data[subset.ix,variable]-subset.inliers.mean)/subset.inliers.sd
subset.inliers.normalized.tf <- (subset.sd.normalized < threshold) & subset
outlier.tf   <- outlier.tf | ((!subset.inliers.normalized.tf) & subset)
inlier.tf    <- inlier.tf | subset.inliers.normalized.tf
subsets[[ix]]  <- subset & subset.inliers.normalized.tf
attr(subsets[[ix]], "name") <- attr(subset, "name")
attr(subsets[[ix]], "mean") <- subset.inliers.mean
attr(subsets[[ix]], "sd") <- subset.inliers.sd
subset.names[ix] <- attr(subset, "name")
mean.char <- as.character(round(subset.inliers.mean,2))
sd.char <- as.character(round(subset.inliers.sd,2))
subset.legends[ix] <- paste0(subset.names[ix],"/Mean: ",mean.char,"/SD: ",sd.char)
}
histbreaks=seq(from=floor(min(data[inlier.tf,variable])), to=ceiling(max(data[inlier.tf,variable])), length.out=21)
for (ix in 1:subset.len) {
hists[[ix]] <- hist(data[subsets[[ix]],variable], breaks=histbreaks, freq=F)
dev.off()
}
plot_filename <- paste(c(subset.names), sep="_vs_")
browser()
plot_filename <- paste0(variable, "_", plot_filename, "_histogram.pdf")
quartzw(file=plot_filename)
for (ix in 1:subset.len) {
if( ix == 1 ) {
    plot(hists[[ix]], col=graph.cols[ix], ylim=c(0.0, 1.0), main=paste0(variable, " Histogram Over All Years"), xlab=paste0(variable, " bins"), freq=F)
} else {
    plot(hists[[ix]], col=graph.cols[ix], ylim=c(0.0, 1.0), xlab=paste0(variable, " bins"), add=T, freq=F)
}
}
legend("topright",legend=subset.legends, bty='n', fill=graph.cols, col=graph.cols)
dev.off()
return(outlier.tf)
}

#Look at general population distribution, across years, etc.  Parents would also be good to look at here.
p1.tf <- grepl("CNJ04.*", cnjpop.df$accession_name)
attr(p1.tf, 'name') <-  "CNJ04"
p2.tf <- grepl("CNJ02.*", cnjpop.df$accession_name)
attr(p2.tf, 'name') <-  "CNJ02"
p1_p1.tf <- grepl("Stevens", cnjpop.df$accession_name) & (cnjpop.df$population == "1")
attr(p1_p1.tf, 'name') <-  "CNJ04-Parents-Stevens"
p1_p2.tf <- grepl("Mullica Queen", cnjpop.df$accession_name) & (cnjpop.df$population == "1")
attr(p1_p2.tf, 'name') <-  "CNJ04-Parent-Mullica Queen"
p2_p1.tf <- grepl("Mullica Queen", cnjpop.df$accession_name) & (cnjpop.df$population == "2")
attr(p2_p1.tf, 'name') <-  "CNJ02-Parent-Mullica Queen"
p2_p2.tf <- grepl("Crimson Queen", cnjpop.df$accession_name) & (cnjpop.df$population == "2")
attr(p2_p2.tf, 'name') <-  "CNJ02-Parent-Crimson Queen"
#First remove obvious outliers by eliminating the extreme quantiles, finding the average and standard deviation of the inluded quantiles, and then removing anything that exceeds 4 sd units.

outlier.threshold <- 1.96
berry_weight.outliers <- assess_population_stats(cnjpop.df,"berry_weight", list(p1.tf,p2.tf), threshold=outlier.threshold, subset_mask=(!berry_weight_na.tf))
berry_weight.outliers <- assess_population_stats(cnjpop.df,"berry_weight", list(p1.tf), threshold=outlier.threshold, subset_mask=(!berry_weight_na.tf))
berry_length.outliers <- assess_population_stats(cnjpop.df,"berry_length", list(p1.tf,p2.tf), threshold=outlier.threshold, subset_mask=(!berry_length_na.tf))
berry_width.outliers <- assess_population_stats(cnjpop.df,"berry_width", list(p1.tf,p2.tf), threshold=outlier.threshold, subset_mask=(!berry_width_na.tf))
berry_num_seed.outliers <- assess_population_stats(cnjpop.df,"num_seeds", list(p1.tf,p2.tf), threshold=outlier.threshold, subset_mask=(!num_seeds_na.tf))

#Make a copy of the cnjpop.df to cleanup outliers (set to NA)
cnjpop.df.cleaned <- cnjpop.df
cnjpop.df.cleaned$berry_weight[berry_weight.outliers] <- NA
cnjpop.df.cleaned$berry_length[berry_length.outliers] <- NA
cnjpop.df.cleaned$berry_width[berry_width.outliers] <- NA
cnjpop.df.cleaned$num_seeds[berry_num_seed.outliers] <- NA

#View the dataset as a scatterplot before doing a 'formal' outlier removal
quartz("Phenotype Scatterplots with Outliers Removed")
#Use regression on highly correlated variables to remove outliers with high residuals
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(r, digits = digits)
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
par(mfrow=c(3,2), las=1, cex=0.8)
pairs(formula=~year+berry_length+berry_width+berry_weight+num_seeds, data=cnjpop.df.cleaned, na.action=na.omit, upper.panel=panel.cor)

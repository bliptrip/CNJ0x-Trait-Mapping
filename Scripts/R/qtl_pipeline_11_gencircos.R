#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
#install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

library(dplyr)
library(qtl)
library(jsonlite)

workflow <- get0("workflow", ifnotfound="../../Workflows/1")
#Which circos traits to render.  This file is in similar format to model-traits.cfg.csv.  All it needs is the following columns: trait, label, mask.
#Any mask==TRUE fields means these fields aren't rendered in the plot.
#
#Default
circostraits.cfg="model-traits.cfg.csv"

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

#Specify to the input the traits we care about rendering in the circos plot.
circostraits.cfg.df <- read.csv(file=paste0(workflow,'/configs/',circostraits.cfg),header=T) %>% filter(is.na(mask) | (mask != "TRUE"))

source('./usefulFunctions.R')
source(paste0(workflow,"/configs/model.cfg"))

#2D LOD Profile Plots

#Need the supermap file to get marker info.
supermap.bin.df <- readRDS(file=geno_rpath2fpath(paste0(geno_consensus_file,".rds")))

extract_bin <- function(mname) {
        bin <- sub('@', ".", mname, fixed=T)
        return(bin)
}
extract_effects <- function(binname, markername, cross, trait) {
binname <- extract_bin(binname)
m.effects <- effectplot(cross, pheno.col=trait, mname1=binname, draw=FALSE)
names(m.effects$Means) <- gsub(paste0(binname,"."),"",names(m.effects$Means))
names(m.effects$SEs) <- gsub(paste0(binname,"."),"",names(m.effects$SEs))
    return(m.effects)
}

#TODO: Go through the following and figure out how I want to render effect plots for the data
#effects
#json.l[['effects']] <- vector("list",length(unlist(markers.l)))
#names(json.l[['effects']]) <- unlist(markers.l)
#k <- 1
#for( chr in chrs ) {
#    markers <- as.character(markers.l[[chr]])
#    bins    <- as.character(bins.l[[chr]])
#    for( l in 1:length(markers) ) {
#        json.l[['effects']][[markers[l]]] <- extract_effects(bins[l], markers[l], cross, trait)
#    }
#}

#trait groups, subgroups, and perform makeqtl() and fitqtl().
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
#Remove all masked entries
traits.df <- traits.df[which(is.na(traits.df$mask) | (traits.df$mask != "TRUE") ),]

exportLODProfile <- function(scan.obj, qtl_type, model, trait) {
    lodprofs         <- attr(scan.obj, "lodprofile")
    num_lod_entries  <- length(unlist(lodprofs))/ncol(lodprofs[[1]])
    model.v          <- vector("character", num_lod_entries)
    trait.v          <- rep(trait, num_lod_entries)
    chr.v            <- vector("numeric", num_lod_entries)
    lod.v            <- vector("numeric", num_lod_entries)
    position.v       <- vector("numeric", num_lod_entries)
    #nearest_marker.v <- vector("character", num_lod_entries)
    current_index    <- 1
    for( lodprof in lodprofs ) {
        lodprof.len                         <- nrow(lodprof)
        end_index                           <- current_index+lodprof.len-1
        model.v[current_index:end_index]    <- rep(model, nrow(lodprof))
        chr.v[current_index:end_index]      <- as.numeric(as.character(lodprof$chr))
        lod.v[current_index:end_index]      <- lodprof$lod
        position.v[current_index:end_index] <- lodprof$pos
        current_index                       <- end_index + 1
    }
    lodprofs <- list(method=rep(qtl_type,length(model.v)), model=model.v, trait=trait.v, chr=chr.v, lod=lod.v, position=position.v)
    return(lodprofs)
}

exportLODProfiles <- function(models, trait) {
    lodprofs.l <- vector("list",length(models))
    for( i in 1:length(models) ) {
        model <- models[i]
        trait_subfolder                 <- paste0(model,"--",trait);
        trait_subfolder_fpath           <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        trait_subsubfolder_fpath        <- file.path(paste0(trait_subfolder_fpath, "/", trait))
        scan.sw  <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.rds'))
        scan.sw.lodprofs <- exportLODProfile(scan.sw, "stepwiseqtl", model, trait)
        trait.sw.lods.df <- data.frame(scan.sw.lodprofs)
        trait.sw.merged.lods.df <- group_by(trait.sw.lods.df, method, model, trait, chr, position) %>% summarize(lod=max(lod))
        scan.one.qtl <- readRDS(paste0(trait_subsubfolder_fpath,'/scanone.qtl.rds'))
        scan.one.lodprofs <- exportLODProfile(scan.one.qtl, "scanone", model, trait)
        trait.one.lods.df <- data.frame(scan.one.lodprofs)
        trait.one.merged.lods.df <- group_by(trait.one.lods.df, method, model, trait, chr, position) %>% summarize(lod=max(lod))
        #Due to overlapping positions b/w different detected QTLs with different LOD scores (how?), I've decided to take the max lod score of either where they overlap to make the graphs look continuous.
        lodprofs.l[[i]] <- rbind(trait.sw.merged.lods.df, trait.one.merged.lods.df)
    }
    names(lodprofs.l) <- models
    lod_file <- paste0(trait,"--lodprofiles.json")
    write_json(lodprofs.l, paste0(workflow,"/configs/circos/",lod_file), auto_unbox=T, pretty=T)
    return(circosfile2path(lod_file))
}

#Consider looping through all QTL's and only keep those that are consistent across all years?
# finding QTL in at least two years
qtl.collated.df <-  read.csv(file=paste0(workflow,"/traits/qtl_collated.consensus.csv"),header=T)


# Circos.js

# Generate Karyotype File for Cranberry using consensus map
karyotype<-numeric()
k.chrs <- 1:12
k.ids <- paste0("vm",k.chrs)
k.labs <-  paste0("LG",k.chrs)
k.cols <- rep("rgb(150,150,150)", length(k.chrs))
k.lens <- (arrange(supermap.bin.df, LG) %>% group_by(LG) %>% summarize(max_consensus=max(consensus)*1000))$max_consensus
karyotype.df <- data.frame(id=k.ids, label=k.labs, color=k.cols, len=k.lens)
karyotype.json <- toJSON(karyotype.df, pretty=T)
write(karyotype.json, file = file.path(paste0(workflow,"/traits/plots/circos/karyotype.json")))

#Generate your colors based on the number of traits * number of unique models
num_traits <- nrow(select(qtl.collated.df, trait) %>% group_by(trait) %>% summarize(count=n()))
num_models <- length(unique(qtl.collated.df$model))
mycols<-matrix(paste0(hueGen(num_traits*num_models,0,340)),nrow=num_models,byrow=F)
model.cols <- colorRampPalette(c("gray20","gray39"))(num_models)

#This was for 5 traits, 3 years, 2 months -- he had a 3 diferent levels for three years, 2 sublevels for each 
#In my case, I will sort by trait, then model

build_label <- function(trait) {
    trait=camel(as.character(trait))
    label <- paste0('vm12 25000 75000 ',trait,' color=black,svgvisibility=hidden,svgclass=qtl_label,svgid=',trait)
    return(label)
}

extract_bin <- function(chr, position) {
    bin.pos <- vector("character", length(chr))
    for( i in 1:length(chr) ) {
        bin.pos[i]   <- paste0("bin_",chr,"@",position,"cM")
    }
    return(bin.pos)
}

extract_marker <- function(chr, position, nearest.marker) {
    bin.pos          <- extract_bin(chr, position)
    marker           <- vector("character", length(chr))
    for( i in 1:length(chr) ) {
        marker[i]   <- gsub(bin.pos[i], "", nearest.marker[i])
    }
    return(marker)
}

gen_color <- function(model_idx, trait) {
    #m.idx <- as.numeric(model)
    m.idx <- model_idx
    t.idx <- as.numeric(trait)
    colors <- vector("character", length(model_idx))
    for( i in 1:length(model_idx) ) {
        colors[i] <- mycols[m.idx[i],t.idx[i]]
    }
    return(colors)
}

qtl.collated.augmented.df <- qtl.collated.df %>%
                        filter(is.na(chr2) & is.na(position2)) %>% #Filter out only additive components, leaving out pairwise QTL interactions
                        arrange(trait,model) %>% 
                        filter(trait %in% circostraits.cfg.df$trait) %>%
                        mutate(trait=as.factor(trait)) %>%
                        mutate(model_idx=as.numeric(as.factor(model))) %>%
                        mutate(class=paste0("marker_", extract_marker(chr, position, nearest.marker), " ", extract_bin(chr, position), "trait_", trait)) %>%
                        mutate(color=gen_color(model_idx, trait)) %>%
                        mutate(stroke_color=gen_color(model_idx, trait)) %>%
                        mutate(model.cols=model.cols[as.numeric(model)])

#Periods in variable names aren't allowed in javascript, as they have special meaning.  Replace all periods with underscores in the column names.
colnames(qtl.collated.augmented.df) <- gsub(".", "_", colnames(qtl.collated.augmented.df), fixed=T)

qtl_files.p <- newPointer(list())
lod_files.p <- newPointer(list())
                    
append_fake_data <- function(group, key) {
    mkey		  <- key[[1]] #key is a tibble, but reduce down to scalar
    models.v      <- rep(unique(as.character(qtl.collated.augmented.df$model)), length(k.chrs))
    models_idx.v  <- rep(unique(as.numeric(qtl.collated.augmented.df$model_idx)), length(k.chrs))
    chr.v         <- unlist(lapply(k.chrs, function(x) { rep(x,num_models) } ))
    df.nrows      <- length(models.v)
    fakedata.df   <- data.frame(method=rep("fakeqtl",df.nrows), model=models.v, chr=chr.v, position=rep(-1,df.nrows),
                                    chr2=rep(NA,df.nrows), position2=rep(NA,df.nrows), nearest_marker=rep(paste0("marker--fake--trait--",mkey),df.nrows), 
									qtl_lod=rep(0,df.nrows), qtl_pvalue=rep(0,df.nrows),
                                    marker_variance=rep(0,df.nrows), model_variance=rep(0,df.nrows), interval=rep(0,df.nrows), 
									GLRpvalue=rep(0,df.nrows),GxYLRpvalue=rep(0,df.nrows),GZRpvalue=rep(0,df.nrows),GxYZRpvalue=rep(0,df.nrows),
									position_consensus=rep(-1,df.nrows),
                                    class=rep("", df.nrows), model_idx=models_idx.v, color=rep("rgb(0,0,0)",df.nrows), stroke_color=rep("rgb(0,0,0)",df.nrows), 
                                    model_cols=rep("rgb(0,0,0)",df.nrows), trait=rep(mkey, df.nrows))
    qtl_file         <- paste0(mkey, "__circos_qtl_file.csv")
    append.pointer(qtl_files.p, circosfile2path(qtl_file))
	group.df		 <- data.frame(group)
	group.df$trait	 <- mkey
	combined.df		 <- rbind(group.df,fakedata.df)
    write.csv(combined.df, file=paste0(workflow,"/configs/circos/",qtl_file), row.names=F)
    models           <- as.vector(unique(group$model))
	lodprofile_file	 <- exportLODProfiles(models, mkey)
    append.pointer(lod_files.p, lodprofile_file)
	return(combined.df)
}

qtl.collated.grouped.df <- qtl.collated.augmented.df %>% 
                                group_by(trait) %>%
                                group_map(~ append_fake_data(.x, .y)) 
ntraits <- length(unique(qtl.collated.augmented.df$trait))

qtl_files <- unlist(qtl_files.p$value)
lod_files <- unlist(lod_files.p$value)

#Write all the trait datafiles in one file for the javascript file to parse these out.
write_json(qtl_files, file.path(paste0(workflow,"/configs/circos/all_traits.json")), auto_unbox=T, pretty=T)
write_json(lod_files, file.path(paste0(workflow,"/configs/circos/all_lods.json")), auto_unbox=T, pretty=T)

rmin<-seq(.3,1.0,length.out=ntraits+1)
rmax<-rmin-(.01-diff(rmin[1:2]))

rdiff<-diff(rmin[1:2])
rmin1<-rmin
rmax1<-rmin+.7*rdiff

rmin2 <- vector("list", length(rmin1))
rmax2 <- vector("list", length(rmin1))
for(i in 1:length(rmin1)) {
    rm      <- seq(rmax1[i]+.01, rmax[i], length.out=num_models+1)
    rmin2[[i]] <- vector("numeric", num_models)
    rmax2[[i]] <- vector("numeric", num_models)
    for( j in 1:num_models ) {
        rmin2[[i]][[j]] <- rm[j]
        rmax2[[i]][[j]] <- rm[j+1]
    }
}

#TODO: Consider adding a marker density plot on outer ring

#Edit the layout here per the script preferences
layout.json <- read_json(file.path(workflow,'/configs/circos/layout.default.json'))
#layout.json$innerRadius <- 1500 - 150
#layout.json$outerRadius <- 1500 - 130
#Write back
write_json(layout.json, file.path(workflow,'/configs/circos/layout.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=T)

#Loop through the traits and generate scatter configs
scatter.configs.json <- vector("list", ntraits)
stack.configs.json <- vector("list", ntraits)
line.configs.json <- vector("list", ntraits)
for( i in 1:ntraits ) {
    #Edit the scatter layout
    scatter.config.json <- read_json(file.path(workflow,'/configs/circos/scatter.default.json'))
    scatter.config.json$innerRadius  <- rmin1[i]
    scatter.config.json$outerRadius  <- rmax1[i]
    scatter.config.json$max          <- num_models + 1
    #Generate the axes
    axis.template       <- scatter.config.json$axes[[1]]
    scatter.config.json$axes <- vector("list", num_models)
    for(j in 1:num_models) {
        scatter.config.json$axes[[j]]          <- axis.template
        scatter.config.json$axes[[j]]$position <- j
        scatter.config.json$axes[[j]]$color    <- model.cols[j]
    }
    scatter.config.json$backgrounds[[1]]$start <- 0;
    scatter.config.json$backgrounds[[1]]$end   <- num_models+1;
    scatter.config.json$backgrounds[[1]]$color <- paste0("rgb(",paste0(as.numeric(col2rgb("gray45")),collapse=","),")");
    scatter.config.json$backgrounds[[1]]$opacity <- 1;
    scatter.configs.json[[i]] <- scatter.config.json
    #Edit the stack layout
    stack.configs.inner.json <- vector("list", num_models)
    for( j in 1:num_models ) {
        stack.config.json <- read_json(file.path(workflow,'/configs/circos/stack.default.json'))
        stack.config.json$innerRadius  <- rmin2[[i]][[j]]
        stack.config.json$outerRadius  <- rmax2[[i]][[j]]
        #stack.config.json$max          <- (num_models*5) + 1
        stack.config.json$backgrounds[[1]]$start <- 0;
        #stack.config.json$backgrounds[[1]]$end   <- (num_models*5)+1;
        stack.config.json$backgrounds[[1]]$color <- paste0("rgb(",paste0(as.numeric(col2rgb("gray22")),collapse=","),")");
        stack.config.json$backgrounds[[1]]$opacity <- 1;
        stack.configs.inner.json[[j]] <- stack.config.json
    }
    stack.configs.json[[i]] = stack.configs.inner.json

    #LOD Line Configurations
    line.config.json <- read_json(file.path(workflow,'/configs/circos/line.default.json'))
    line.config.json$innerRadius  <- rmin1[i]
    line.config.json$outerRadius  <- rmax1[i]
    line.config.json$backgrounds[[1]]$color <- paste0("rgb(",paste0(as.numeric(col2rgb("gray45")),collapse=","),")");
    line.configs.json[[i]] <- line.config.json
}
#Write the axes
write_json(scatter.configs.json, file.path(workflow,'/configs/circos/scatter.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)
write_json(stack.configs.json, file.path(workflow,'/configs/circos/stack.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)
write_json(line.configs.json, file.path(workflow,'/configs/circos/line.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)

#Save image for reloading later if desired
save.image(".RData.11_gencircos")

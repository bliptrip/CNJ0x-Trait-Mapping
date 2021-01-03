#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#
#NOTE: This particular script runs on my local computer instead of the other version which is meant to be run on the UW HTCondor system.
#install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

library(dplyr)
library(qtl)
library(jsonlite)

workflow="../../Workflows/1"
#Which circos traits to render.  This file is in similar format to model-traits.cfg.csv.  All it needs is the following columns: mtraits, trait, mask.
#Any mask==TRUE fields means these fields aren't rendered in the plot.
#
#Default
circostraits.cfg="circos-traits.default.cfg.csv"

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

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

#For each qtl in the collated file, use it's position and consensus position to calculate the effects.  Store this information in the collated file?
qtl.collated.df <- read.csv(file=paste0(workflow,'/traits/qtl_collated.consensus.csv'), head=TRUE)
qtl.collated.df <- group_by(model, mtraits) %>%
qtl.collated.df <- qtl.collated.df
                for( chr in chrs ) {
                    markers <- as.character(markers.l[[chr]])
                    bins    <- as.character(bins.l[[chr]])
                    for( l in 1:length(markers) ) {
                        json.l[['effects']][[markers[l]]] <- extract_effects(bins[l], markers[l], cross, trait)
                    }
                }
#Export the qtl info to JSON files that can be parsed and displayed by D3 javascript
exportQTL2D3 <- function(traits.df) {
    #These should all be in the same mtrait category.  Select the first entry and get a breakdown of the subtraits, and for each
    #subtrait, plot it's LOD profile.
    traits                    <- unlist(strsplit(traits.df[1,"mtraits"],","))
    trait.names               <- paste0(traits,collapse="__")
    models                    <- unique(traits.df[,"model"])
    for( trait in traits ) {
        postscript(file=file.path(paste0(workflow,"/traits/plots/all-models--", trait.names, "--", trait,".eps")), horizontal=F, paper="letter", onefile=T)
        for( j in 1:length(models) ) {
            model <- models[j]
            traits.subset.df <- traits.df[which(traits.df$model == model),]
            for( i in 1:length(traits.subset.df[,1]) ) {
                json.names                <- c("phenotype","chr","lod","markerindex","markers","effects","phevals","sex","geno","individuals")
                json.l                    <- vector(mode="list", length=length(json.names))
                names(json.l)             <- json.names
                trait.cfg                 <- traits.subset.df[i,]
                trait_subfolder           <- paste0(c(trait.cfg$model,trait.names),collapse="--")
                trait_subfolder_fpath     <- file.path(paste0(workflow,"/traits"), trait_subfolder)
                trait_subsubfolder_fpath  <- file.path(paste0(trait_subfolder_fpath, "/", trait))
                cross                     <- readRDS(file=paste0(trait_subfolder_fpath,"/cross.rds"))
                json.l[["phenotype"]]     <- trait
                scan.sw  <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.RDS'))
                lodprofs <- attr(scan.sw, "lodprofile")
                #Grab list of chromosomes
                chrs <- NULL
                for( qtl in lodprofs ) {
                    chrs <- c(chrs, unique(levels(qtl$chr)))
                }
                chrs     <- unique(chrs)
                json.l[["chr"]] <- chrs
                #Grab pos + lodscores.  This requires synthesizing the separate QTLs that belong to one chromosome together 
                lods <- vector(mode="list", length=length(chrs))
                names(lods) <- chrs
                for( lodprof in lodprofs ) {
                    chr <- as.character(lodprof$chr[1])
                    if( is.null(lods[[chr]]) ) {
                        lods[[chr]] <- select(lodprof, pos, lod)
                    } else {
                        lod.temp <- merge(lods[[chr]],select(lodprof, pos, lod), by="row.names", all=TRUE)
                        lods[[chr]] <- select(lod.temp, pos.x, pos.y, lod.x, lod.y) %>% mutate(pos = apply(cbind(pos.x,pos.y),1,max,na.rm=TRUE)) %>% mutate(lod = apply(cbind(lod.x,lod.y),1,max,na.rm=TRUE)) %>% select(pos, lod)
                        rownames(lods[[chr]]) <- lod.temp$Row.names
                    }
                }
                json.l[["lod"]] <- lods
                #toJSON(lods, dataframe="columns", auto_unbox=TRUE, pretty=T)

                #Markers and marker indices
                #Karl Broman uses marker indices to index markers of interest into their location in the lod profile.  I will only include markers that are in the lod profile range.
                markers.idxs        <- vector(mode="list", length=length(chrs))
                names(markers.idxs) <- chrs
                markers.l           <- vector(mode="list", length=length(chrs))
                names(markers.l)    <- chrs
                bins.l              <- vector(mode="list", length=length(chrs))
                names(bins.l)       <- chrs
                #Get nearest marker info for each chromosome
                for( chr in chrs ) { 
                    supermap.bin.chr.df <- select(supermap.bin.df,marker,LG,consensus,binID) %>% filter(LG==chr)
                    rownames(lods[[chr]]) <- sub('.', "@", rownames(lods[[chr]]),fixed=T)
                    markers.df <- merge(lods[[chr]],supermap.bin.chr.df, by.x="row.names", by.y="binID", all.x=T)
                    markers.idx <- which(!is.na(markers.df$marker))
                    markers     <- markers.df$marker[markers.idx]
                    bins        <- markers.df$Row.names[markers.idx]
                    #markers.df  <- data.frame(matrix(markers.idx, nrow=1, ncol=length(markers.idx)))
                    #colnames(markers.df) <- markers
                    #markers.idxs[[chr]] <- markers.df
                    markers.idxs[[chr]] <- as.list(markers.idx)
                    names(markers.idxs[[chr]]) <- markers
                    markers.l[[chr]] <- markers
                    bins.l[[chr]]    <- bins
                }
                json.l[["markerindex"]] <- markers.idxs
                json.l[["markers"]] <- markers.l
                #toJSON(markers.idxs, dataframe="columns", auto_unbox=TRUE, pretty=T)

                #effects
                json.l[['effects']] <- vector("list",length(unlist(markers.l)))
                names(json.l[['effects']]) <- unlist(markers.l)
                k <- 1
                for( chr in chrs ) {
                    markers <- as.character(markers.l[[chr]])
                    bins    <- as.character(bins.l[[chr]])
                    for( l in 1:length(markers) ) {
                        json.l[['effects']][[markers[l]]] <- extract_effects(bins[l], markers[l], cross, trait)
                    }
                }

                #for each marker in marker.l and the current trait, calculate the effect sizes

                #phenvals
                json.l[['phenvals']] <- cross$pheno[trait]

                #sex

                #geno
                #preprocess
                num_markers <- 0
                for( chr in chrs ) {
                    num_markers <- num_markers + length(markers.l[[chr]])
                }
                marker.names <- vector(mode="character", length=num_markers)
                geno.l <- vector(mode="list", length=num_markers)
                idx <- 1
                for( chr in chrs ) {
                    geno      <- cross$geno[[chr]]$data
                    markers   <- markers.l[[chr]]
                    bins      <- bins.l[[chr]]
                    genos.tf  <- sub('.', '@', colnames(geno), fixed=TRUE) %in% bins
                    geno      <- geno[,genos.tf]
                    colnames(geno) <- markers
                    for( i in 1:ncol(geno) ) {
                        geno.l[[idx]] <- geno[,i]
                        marker.names[idx] <- as.character(markers[i])
                        idx <- idx + 1
                    }
                }
                names(geno.l) <- marker.names
                json.l[["geno"]] <- geno.l

                #individuals
                #Need to grab this from original dataset that generated cross file

                json.raw <- toJSON(json.l, dataframe="columns", auto_unbox=TRUE, pretty=T)
                write(json.raw, paste0(trait_subsubfolder_fpath,'/trait_qtl_info.json'))
                scan.sw  <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.RDS'))
            }
        }
        dev.off()
    }
}

#trait groups, subgroups, and perform makeqtl() and fitqtl().
traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
#Remove all masked entries
traits.df <- traits.df[which(is.na(traits.df$mask) | (traits.df$mask != "TRUE") ),]
#Apply an 'aggregate' function whereby we organize by mtrait/trait, and then within these subsets, we plot a set of lod profiles across all models (years)
#by(traits.df, INDICES=traits.df[, "mtraits"], myPlotLodProfiles)

#exportQTL2D3(traits.df)

exportLODProfiles <- function(models, mtraits, trait) {
    #These should all be in the same mtrait category.  Select the first entry and get a breakdown of the subtraits, and for each
    #subtrait, plot it's LOD profile.
    traits                    <- unlist(strsplit(mtraits,","))
    trait.names               <- paste0(traits,collapse="__")
    #Preprocess
    num_lod_entries <- 0
    lodprofs.l <- vector("list",length(models))
    for( i in 1:length(models) ) {
        model <- models[i]
        trait_subfolder                 <- paste0(model,"--",trait.names);
        trait_subfolder_fpath           <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        trait_subsubfolder_fpath        <- file.path(paste0(trait_subfolder_fpath, "/", trait))
        scan.sw  <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.RDS'))
        lodprofs.l[[i]] <- attr(scan.sw, "lodprofile")
        num_lod_entries <- num_lod_entries + (length(unlist(lodprofs.l[[i]]))/ncol(lodprofs.l[[i]][[1]]))
    }
    model.v          <- vector("character", num_lod_entries)
    mtraits.v        <- rep(trait.names, num_lod_entries)
    trait.v          <- rep(trait, num_lod_entries)
    chr.v            <- vector("numeric", num_lod_entries)
    lod.v            <- vector("numeric", num_lod_entries)
    position.v       <- vector("numeric", num_lod_entries)
    nearest_marker.v <- vector("character", num_lod_entries)
    lodprofs.data.l  <- vector("list",length(models))
    current_index    <- 1
    for( i in 1:length(models) ) {
        model     <- models[i]
        lodprofs  <- lodprofs.l[[i]]
        for( lodprof in lodprofs ) {
            lodprof.len                         <- nrow(lodprof)
            end_index                           <- current_index+lodprof.len-1
            model.v[current_index:end_index]    <- rep(model, nrow(lodprof))
            chr.v[current_index:end_index]      <- as.numeric(as.character(lodprof$chr))
            lod.v[current_index:end_index]      <- lodprof$lod
            position.v[current_index:end_index] <- lodprof$pos
            current_index                       <- end_index + 1
        }
    }
    trait.lods.df <- data.frame(model=model.v, mtraits=mtraits.v, trait=trait.v, chr=chr.v, lod=lod.v, position=position.v)
    #Due to overlapping positions b/w different detected QTLs with different LOD scores (how?), I've decided to take the max lod score of either where they overlap to make the graphs look continuous.
    trait.merged.lods.df <- group_by(trait.lods.df, model, mtraits, trait, chr, position) %>% summarize(lod=max(lod))
    lod_file <- paste0(workflow,"/configs/circos/",mtraits,"--",trait,"--lodprofiles.csv")
    write.csv(trait.merged.lods.df, file=lod_file, row.names=F)
    return(lod_file)
}

exportLODProfiles2 <- function(models, mtraits, trait) {
    #These should all be in the same mtrait category.  Select the first entry and get a breakdown of the subtraits, and for each
    #subtrait, plot it's LOD profile.
    traits                    <- unlist(strsplit(mtraits,","))
    trait.names               <- paste0(traits,collapse="__")
    #Preprocess
    lodprofs.l <- vector("list",length(models))
    for( i in 1:length(models) ) {
        model <- models[i]
        trait_subfolder                 <- paste0(model,"--",trait.names);
        trait_subfolder_fpath           <- file.path(paste0(workflow,"/traits"), trait_subfolder)
        trait_subsubfolder_fpath        <- file.path(paste0(trait_subfolder_fpath, "/", trait))
        scan.sw  <- readRDS(paste0(trait_subsubfolder_fpath,'/scansw.RDS'))
        lodprofs         <- attr(scan.sw, "lodprofile")
        num_lod_entries  <- length(unlist(lodprofs))/ncol(lodprofs[[1]])
        model.v          <- vector("character", num_lod_entries)
        mtraits.v        <- rep(trait.names, num_lod_entries)
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
        trait.lods.df <- data.frame(model=model.v, mtraits=mtraits.v, trait=trait.v, chr=chr.v, lod=lod.v, position=position.v)
        #Due to overlapping positions b/w different detected QTLs with different LOD scores (how?), I've decided to take the max lod score of either where they overlap to make the graphs look continuous.
        lodprofs.l[[i]] <- trait.merged.lods.df <- group_by(trait.lods.df, model, mtraits, trait, chr, position) %>% summarize(lod=max(lod))
    }
    names(lodprofs.l) <- models
    lod_file <- paste0(mtraits,"--",trait,"--lodprofiles.json")
    write_json(lodprofs.l, paste0(workflow,"/configs/circos/",lod_file), auto_unbox=T, pretty=T)
    return(circosfile2path(lod_file))
}

#Consider looping through all QTL's and only keep those that are consistent across all years?
# finding QTL in at least two years
qtl.collated.df <- read.csv(file=paste0(workflow,"/traits/qtl_collated.csv"),header=T)

# Circos.js

# Generate Karyotype File for Cranberry using consensus map
karyotype<-numeric()
k.chrs <- 1:12
k.ids <- paste0("vm",k.chrs)
k.labs <-  paste0("chr",k.chrs)
k.cols <- rep("rgb(82,82,82)", length(k.chrs))
k.lens <- (arrange(supermap.bin.df, LG) %>% group_by(LG) %>% summarize(max_consensus=max(consensus)*1000))$max_consensus
karyotype.df <- data.frame(id=k.ids, label=k.labs, color=k.cols, len=k.lens)
karyotype.json <- toJSON(karyotype.df, pretty=T)
write(karyotype.json, file = file.path(paste0(workflow,"/traits/plots/circos/karyotype.json")))

#Generate your colors based on the number of traits * number of unique models
num_traits <- nrow(select(qtl.collated.df, mtraits, trait) %>% group_by(mtraits,trait) %>% summarize(count=n()))
num_models <- length(unique(qtl.collated.df$model))
mycols<-matrix(paste0(hueGen(num_traits*num_models,0,340)),nrow=num_models,byrow=F)
model.cols <- colorRampPalette(c("gray20","gray39"))(num_models)

#This was for 5 traits, 3 years, 2 months -- he had a 3 diferent levels for three years, 2 sublevels for each 
#In my case, I will sort by mtraits, trait, then model

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
    return(marker[i])
}

gen_color <- function(model, mtraits_trait) {
    m.idx <- as.numeric(model)
    t.idx <- as.numeric(mtraits_trait)
    colors <- vector("character", length(model))
    for( i in 1:length(model) ) {
        colors[i] <- mycols[m.idx[i],t.idx[i]]
    }
    return(colors)
}

#Specify to the input the traits we care about rendering in the circos plot.
circostraits.cfg.df <- read.csv(file=paste0(workflow,'/configs/circos/',circostraits.cfg),header=T) %>% filter(is.na(mask) | (mask != "TRUE")) %>% mutate(mtraits_trait=paste0(mtraits,trait))

supermap.bin.chr.df <- arrange(qtl.collated.df,mtraits,trait,model) %>% 
                    mutate(mtraits_trait=paste0(mtraits,trait)) %>%
                    filter(mtraits_trait %in% circostraits.cfg.df$mtraits_trait) %>%
                    mutate(mtraits_trait=as.factor(mtraits_trait)) %>%
                    mutate(model_idx=as.numeric(model)) %>%
                    mutate(class=paste0("marker_", extract_marker(chr, position, nearest.marker), " ", extract_bin(chr, position), "trait_", trait)) %>%
                    mutate(color=gen_color(model, mtraits_trait)) %>%
                    mutate(stroke_color=gen_color(model, mtraits_trait)) %>%
                    mutate(model.cols=model.cols[as.numeric(model)])

#Periods in variable names aren't allowed in javascript, as they have special meaning.  Replace all periods with underscores in the column names.
colnames(supermap.bin.chr.df) <- gsub(".", "_", colnames(supermap.bin.chr.df), fixed=T)
                    
supermap.bin.grouped <- supermap.bin.chr.df %>% group_by(mtraits,trait)

indices <- attr(supermap.bin.grouped,"indices")
labels  <- attr(supermap.bin.grouped,"labels")
labels$mtraits <- as.character(labels$mtraits)
labels$trait <- as.character(labels$trait)
qtl_files <- vector("character", length(indices))
lod_files <- vector("character", length(indices))
#Generate a fake dataset to have all tracks display
models.v      <- rep(unique(as.character(supermap.bin.chr.df$model)), length(k.chrs))
models_idx.v  <- rep(unique(as.numeric(supermap.bin.chr.df$model)), length(k.chrs))
chr.v         <- unlist(lapply(k.chrs, function(x) { rep(x,num_models) } ))
df.nrows      <- length(models.v)
fakedata.df   <- data.frame(model=models.v, year=rep("",df.nrows), chr=chr.v, position=rep(-1,df.nrows),
                            marker_variance=rep(0,df.nrows), model_variance=rep(0,df.nrows), interval=rep(0,df.nrows), class=rep("", df.nrows),
                            model_idx=models_idx.v, color=rep("transparent",df.nrows), stroke_color=rep("transparent",df.nrows), model_cols=rep("transparent",df.nrows))
for( i in 1:length(indices) ) {
    #To avoid rendering issues where there are not at least multiple values for each tile/scatterplot in data file, inject 'hidden' elements to get around this.
    traits           <- unlist(strsplit(labels$mtraits[i],","))
    trait.names      <- paste0(traits,collapse="__")
    trait_prefix     <- paste0(c(trait.names, labels$trait[i]),collapse="--")
    fakedata.mut.df  <- mutate(fakedata.df, trait=rep(labels$trait[i], df.nrows), mtraits=rep(labels$mtraits[i], df.nrows),
                                                      mtraits_trait=rep(paste0(labels$mtraits[i],labels$trait[i]), df.nrows),
                                                      nearest_marker=rep(paste0("marker--fake--trait--",labels$trait[i]), df.nrows))
    qtl_file		 <- paste0(trait_prefix, "__circos_qtl_file.csv")
    qtl_files[i]     <- circosfile2path(qtl_file)
    write.csv(rbind(data.frame(supermap.bin.grouped[indices[[i]]+1,]),fakedata.mut.df), file=paste0(workflow,"/configs/circos/",qtl_file), row.names=F)
    models           <- as.vector(unique(supermap.bin.grouped[indices[[i]]+1, "model"])$model)
    lod_files[i]     <- exportLODProfiles2(models, labels$mtraits[i], labels$trait[i])
}

#Write all the trait datafiles in one file for the javascript file to parse these out.
write_json(qtl_files, file.path(paste0(workflow,"/configs/circos/all_traits.json")), auto_unbox=T, pretty=T)
write_json(lod_files, file.path(paste0(workflow,"/configs/circos/all_lods.json")), auto_unbox=T, pretty=T)

rmin<-seq(.3,1.0,length.out=length(indices)+1)
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

#Edit the layout here per the script preferences
layout.json <- read_json(file.path(workflow,'/configs/circos/layout.default.json'))
layout.json$innerRadius <- 1500 - 150
layout.json$outerRadius <- 1500 - 130
#Write back
write_json(layout.json, file.path(workflow,'/configs/circos/layout.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=T)

#Loop through the traits and generate scatter configs
scatter.configs.json <- vector("list", length(indices))
stack.configs.json <- vector("list", length(indices))
line.configs.json <- vector("list", length(indices))
for( i in 1:length(indices) ) {
    #Edit the scatter layout
    scatter.config.json <- read_json(file.path(workflow,'/configs/circos/scatter.default.json'))
    scatter.config.json$innerRadius  <- rmin1[i]
    scatter.config.json$outerRadius  <- rmax1[i]
    scatter.config.json$max          <- num_models + 1
    #Generate the axes
    axis.template       <- scatter.config.json$axes[[1]]
    scatter.config.json$axes <- vector("list", num_models)
    for(j in 1:num_models) {
        scatter.config.json$axes[[j]]         <- axis.template
        scatter.config.json$axes[[j]]$start   <- j
        scatter.config.json$axes[[j]]$color   <- model.cols[j]
    }
    scatter.config.json$backgrounds[[1]]$start <- 0;
    scatter.config.json$backgrounds[[1]]$end   <- num_models+1;
    scatter.config.json$backgrounds[[1]]$color <- col2rgb("gray15");
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
        stack.config.json$backgrounds[[1]]$color <- col2rgb("gray15");
        stack.config.json$backgrounds[[1]]$opacity <- 1;
        stack.configs.inner.json[[j]] <- stack.config.json
    }
    stack.configs.json[[i]] = stack.configs.inner.json

    #LOD Line Configurations
    line.config.json <- read_json(file.path(workflow,'/configs/circos/line.default.json'))
    line.config.json$innerRadius  <- rmin1[i]
    line.config.json$outerRadius  <- rmax1[i]
    line.configs.json[[i]] <- line.config.json
}
#Write the axes
write_json(scatter.configs.json, file.path(workflow,'/configs/circos/scatter.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)
write_json(stack.configs.json, file.path(workflow,'/configs/circos/stack.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)
write_json(line.configs.json, file.path(workflow,'/configs/circos/line.configs.json'), simplifyVector=TRUE, auto_unbox=TRUE, pretty=TRUE)

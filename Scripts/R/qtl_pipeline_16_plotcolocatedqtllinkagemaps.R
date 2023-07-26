#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(conflicted)
library(RColorBrewer)
library(tidyverse)
devtools::load_all("~/software/bliptrip/LinkageMapView")

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

cnj02_workflow="../../Workflows/9"
cnj04_workflow="../../Workflows/10"

source('./usefulFunctions.R')

all_qtl = read_csv("../../Data/publication/tables/all.raw.consensus.csv")

model_config = read_csv("../../Data/publication/tables/model-traits.cfg.csv") %>%
                mutate(trait=factor(trait,levels=unique(trait),ordered=TRUE),
                       label_short=factor(label_short,levels=unique(label_short),ordered=TRUE))
model_config_trait_levels = levels(model_config$trait)
model_config_short_levels = levels(model_config$label_short)

genetic_map <-  read_csv("../../Data/genetic_data/DerivedData/consensusMapAll2_withGenes.csv") %>%
                    rename(position=consensus) %>%
                    filter(!is.na(LG)) %>%
                    mutate(chr=paste0("lg",LG)) %>%
                    select(marker,position,chr)

genetic_map_lg_summary <- genetic_map %>%
                            group_by(chr) %>%
                            summarize(min=min(position),
                                      max=max(position))

#Normally this would be a genetic marker map -- Replace with relevant individual trait QTLs
columns_highlighted = c("trait","chr","model","position","marker_variance","study")

qtls_cnj02 = read_csv(paste0(cnj02_workflow,"/traits/qtl_collated.consensus.csv")) %>%
                filter((chr2 == '') | is.na(chr2)) %>%
                rename(qtl_lod="qtl.lod",
                       marker_variance="marker.variance",
                       nearest_marker="nearest.marker",
                       interval_left="interval.left",
                       interval_right="interval.right") %>%
                mutate(study="c2") %>%
                mutate(trait=factor(trait,levels=model_config_trait_levels,ordered=TRUE)) %>%
                mutate(trait=factor(model_config_short_levels[as.integer(trait)],levels=model_config_short_levels,ordered=TRUE)) %>%
                mutate(model=as.character(model)) %>%
                select(any_of(columns_highlighted))

h2_cnj02_raw = read_csv(paste0(cnj02_workflow,"/traits/h2.csv"),show_col_types=FALSE)
h2_cnj02 = h2_cnj02_raw %>%
            group_by(label_short) %>%
            summarize(h2=mean(h2))


h2_cnj04_raw = read_csv(paste0(cnj04_workflow,"/traits/h2.csv"),show_col_types=FALSE)
h2_cnj04 = h2_cnj04_raw %>%
            group_by(label_short) %>%
            summarize(h2=mean(h2))

qtls_cnj04 = read_csv(paste0(cnj04_workflow,"/traits/qtl_collated.consensus.csv")) %>%
                filter((chr2 == '') | is.na(chr2)) %>%
                rename(qtl_lod="qtl.lod",
                       marker_variance="marker.variance",
                       nearest_marker="nearest.marker",
                       interval_left="interval.left",
                       interval_right="interval.right") %>%
                mutate(study="c4") %>%
                mutate(trait=factor(trait,levels=model_config_trait_levels,ordered=TRUE)) %>%
                mutate(trait=factor(model_config_short_levels[as.integer(trait)],levels=model_config_short_levels,ordered=TRUE)) %>%
                mutate(model=as.character(model)) %>%
                select(any_of(columns_highlighted))

qtls_schlautman = read_csv("../../Data/publication/tables/schlautman2015.qtls.csv") %>%
                    mutate(study="s") %>%
                    mutate(model=as.character(model)) %>%
                    select(columns_highlighted)

qtls_diazGarciaImage = read_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.qtls.csv") %>%
                        mutate(study="dgi") %>%
                        mutate(model=as.character(model)) %>%
                        select(columns_highlighted)

qtls_diazGarciaMassive = read_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.qtls.csv") %>%
                            mutate(study="dgm") %>%
                            mutate(model=as.character(model)) %>%
                            select(columns_highlighted)

colocated_cnj02 = read_csv("../../Data/publication/tables/cnj02_qtl_collated.consensus.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 3 &
                          mean_marker_variance >= 10) %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj02, "../../Data/publication/tables/cnj02_qtl_collated.consensus.filtered.csv")

colocated_cnj02_nayay = read_csv("../../Data/publication/tables/cnj02_qtl_collated.consensus.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 2 &
                          mean_marker_variance >= 10 &
                          grepl("all-years",model))  %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj02_nayay, "../../Data/publication/tables/cnj02_qtl_collated.consensus.nayay.filtered.csv")

colocated_cnj02_nay = read_csv("../../Data/publication/tables/cnj02_qtl_collated.consensus.nay.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 2 &
                          mean_marker_variance >= 10) %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj02_nay, "../../Data/publication/tables/cnj02_qtl_collated.consensus.nay.filtered.csv")
    
min_mean_marker_variance <- min(colocated_cnj02$mean_marker_variance)
max_mean_marker_variance <- max(colocated_cnj02$mean_marker_variance)
med_mean_marker_variance <- median(colocated_cnj02$mean_marker_variance)
median(colocated_cnj02$mean_marker_variance)
max_model_count <- max(colocated_cnj02$model_count)

colocated_cnj04 = read_csv("../../Data/publication/tables/cnj04_qtl_collated.consensus.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 3 &
                          mean_marker_variance >= 10) %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj04, "../../Data/publication/tables/cnj04_qtl_collated.consensus.filtered.csv")

colocated_cnj04_nayay = read_csv("../../Data/publication/tables/cnj04_qtl_collated.consensus.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 2 &
                          mean_marker_variance >= 10 &
                          grepl("all-years",model))  %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj04_nayay, "../../Data/publication/tables/cnj04_qtl_collated.consensus.nayay.filtered.csv")

colocated_cnj04_nay = read_csv("../../Data/publication/tables/cnj04_qtl_collated.consensus.nay.csv", show_col_types=FALSE) %>%
                   filter(model_count >= 2 &
                          mean_marker_variance > 10) %>%
                   arrange(desc(mean_marker_variance))
write_csv(colocated_cnj04_nay, "../../Data/publication/tables/cnj04_qtl_collated.consensus.nay.filtered.csv")
                          
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_cnj04$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_cnj04$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_cnj04$model_count))
                           
colocated_cnj02_grouped = read_csv("../../Data/publication/tables/cnj02_qtl_collated.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (model_count >= 4)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_cnj02_grouped, "../../Data/publication/tables/cnj02_qtl_collated.consensus.grouped.filtered.csv")

#all_traits <- unique(unlist(strsplit(colocated_cnj02_grouped$trait,"+",fixed=TRUE)))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_cnj02_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_cnj02_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_cnj02_grouped$model_count))

colocated_cnj04_grouped = read_csv("../../Data/publication/tables/cnj04_qtl_collated.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (model_count >= 4) &
                                   (trait_count >= 1)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_cnj04_grouped, "../../Data/publication/tables/cnj04_qtl_collated.consensus.grouped.filtered.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_cnj04_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_cnj04_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_cnj04_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_cnj04_grouped$model_count))

colocated_diazGarciaImage_grouped = read_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.grouped.consensus.qtls.csv") %>%
                            filter((mean_marker_variance >= 8) &
                                   (model_count >= 2)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_diazGarciaImage_grouped, "../../Data/publication/tables/diazGarcia2018ImagePhenotyping_qtl_collated.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_diazGarciaImage_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_diazGarciaImage_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_diazGarciaImage_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_diazGarciaImage_grouped$model_count))

colocated_diazGarciaMassive_grouped = read_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.grouped.consensus.qtls.csv") %>%
                            filter(mean_marker_variance >= 10) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_diazGarciaMassive_grouped, "../../Data/publication/tables/diazGarcia2018MassivePhenotyping_qtl_collated.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_diazGarciaMassive_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_diazGarciaMassive_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_diazGarciaMassive_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_diazGarciaMassive_grouped$model_count))

colocated_cnj0x_grouped = read_csv("../../Data/publication/tables/cnj0x_qtl_collated.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (model_count >= 4) &
                                   (population_count >= 2)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_cnj0x_grouped, "../../Data/publication/tables/cnj0x_qtl_collated.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_cnj0x_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_cnj0x_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_cnj0x_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_cnj0x_grouped$model_count))

colocated_diazGarciaImage_cnj0x_grouped = read_csv("../../Data/publication/tables/diazGarciaImagePhenotyping2018_cnj0x.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (study_count >= 2) &
                                   (trait_count >= 3)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_diazGarciaImage_cnj0x_grouped, "../../Data/publication/tables/diazGarciaImagePhenotyping2018_cnj0x.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_diazGarciaImage_cnj0x_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_diazGarciaImage_cnj0x_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_diazGarciaImage_cnj0x_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_diazGarciaImage_cnj0x_grouped$model_count))

colocated_diazGarciaMassive_cnj0x_grouped = read_csv("../../Data/publication/tables/diazGarciaMassivePhenotyping2018_cnj0x.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (study_count >= 2) &
                                   (population_count >= 2)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_diazGarciaMassive_cnj0x_grouped, "../../Data/publication/tables/diazGarciaMassivePhenotyping2018_cnj0x.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_diazGarciaMassive_cnj0x_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_diazGarciaMassive_cnj0x_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_diazGarciaMassive_cnj0x_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_diazGarciaMassive_cnj0x_grouped$model_count))

colocated_schlautman_cnj0x_grouped = read_csv("../../Data/publication/tables/schlautman2015_cnj0x.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (study_count >= 2)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_schlautman_cnj0x_grouped, "../../Data/publication/tables/schlautman2015_cnj0x.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_schlautman_cnj0x_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_schlautman_cnj0x_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_schlautman_cnj0x_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_schlautman_cnj0x_grouped$model_count))

colocated_schlautman_diazGarciaImage_cnj0x_grouped = read_csv("../../Data/publication/tables/schlautman2015_diazGarciaImagePhenotyping2018_cnj0x.grouped.consensus.csv") %>%
                            filter((mean_marker_variance >= 10) &
                                   (study_count >= 2) &
                                   (trait_count >= 3)) %>%
                            arrange(desc(mean_marker_variance))
write_csv(colocated_schlautman_diazGarciaImage_cnj0x_grouped, "../../Data/publication/tables/schlautman2015_diazGarciaImagePhenotyping2018_cnj0x.filtered.consensus.csv")

#all_traits <- unique(c(all_traits,unlist(strsplit(colocated_schlautman_diazGarciaImage_cnj0x_grouped$trait,"+",fixed=TRUE))))
min_mean_marker_variance <- min(c(min_mean_marker_variance,colocated_schlautman_diazGarciaImage_cnj0x_grouped$mean_marker_variance))
max_mean_marker_variance <- max(c(max_mean_marker_variance,colocated_schlautman_diazGarciaImage_cnj0x_grouped$mean_marker_variance))
max_model_count <- max(c(max_model_count,colocated_schlautman_diazGarciaImage_cnj0x_grouped$model_count))

#Color based on single trait
#Brightness based on number of models associated with colocated QTL
#Saturation based on QTL variance explained -- Need to get min and max variance explained and scale/interpolate based on this

blend_hues <- function(traits, tlevels, cm_rgb) {
    unlist(lapply(strsplit(traits,"+",fixed=TRUE), function(x) {   a <- cm_rgb %>% filter(tlevels %in% x);
                                                            return(mean(rgb2hsv(a$r,a$g,a$b)[1,])) }))
}

generate_brightness <- function(nmodels, max_nmodels, base=0.375, ceil=0.75) {
    bdiff = ceil-base
    binc = bdiff/(max_nmodels-1)
    return(base + (nmodels-1) * binc)
}

saturation_model = lm(saturation ~ mv, tibble(mv=c(min_mean_marker_variance,max_mean_marker_variance), saturation=c(0.25, 1.00)))
generate_saturation <- function(mean_marker_variance,mm=saturation_model) {
    predict(mm, list(mv=mean_marker_variance))
}

colocated_all <- colocated_cnj02 %>%
                    bind_rows(colocated_cnj04) %>%
                    bind_rows(colocated_cnj02_grouped) %>%
                    bind_rows(colocated_cnj04_grouped) %>%
                    bind_rows(colocated_cnj0x_grouped) %>%
                    bind_rows(colocated_diazGarciaImage_grouped) %>%
                    bind_rows(colocated_diazGarciaMassive_grouped) %>%
                    bind_rows(colocated_diazGarciaImage_cnj0x_grouped) %>%
                    bind_rows(colocated_diazGarciaMassive_cnj0x_grouped) %>%
                    bind_rows(colocated_schlautman_cnj0x_grouped) %>%
                    bind_rows(colocated_schlautman_diazGarciaImage_cnj0x_grouped)

all_traits <- factor(unique(unlist(lapply(colocated_all$trait,function(x) { strsplit(x,'+',fixed=T) }))))
trait_levels = levels(all_traits)
colormap.blind <- read_tsv("../../Data/publication/tables/24.color.blindness.palette.txt", skip=10)
colormap.main  <- colormap.blind %>%
                    filter(type == 'main') %>%
                    head(n=length(trait_levels))
colormap.alt   <- colormap.blind %>%
                    filter(type == 'alt') %>%
                    head(n=nrow(colormap.main)-length(trait_levels))
colormap <- colormap.main %>%
                bind_rows(colormap.alt) %>%
                mutate(tcol=rgb(r/255.,g/255.,b/255.))
if( nrow(colormap) < length(trait_levels) ) {
    warning("Number of traits (n=",length(trait_levels),") exceeds colormap size (m=",nrow(colormap),").  Some colors will be repeated.")
    colorr <- TRUE
} else {
    colorr <- FALSE
}
set.seed(0xEEF7331) #For consistent results in sample.int()
trait_rgb = colormap[sample.int(n=nrow(colormap),size=length(trait_levels),replace=colorr),]

generate_linkagemapplots_primitive <- function(qtl_collocated, all_qtl_raw, filename, lgsums, lg.col="white", alpha=0.4, jitterloci=2, lgw=0.25, render.type='png', sectcolindex=NULL, ...) {
    qtl_collocated <- qtl_collocated %>%
                            mutate(hue=blend_hues(trait,trait_levels,trait_rgb),
                                   brightness=generate_brightness(model_count, max_model_count),
                                   saturation=generate_saturation(mean_marker_variance)) %>%
                             mutate(color=hsv(hue,saturation,brightness))
    qtl_collocated_raw_ids <- as.numeric(unlist(sapply(qtl_collocated$id, function(x) { strsplit(x,"+",fixed=T) })))
    #Filter out only relevant QTL
    qtl_raw  <-   all_qtl_raw %>%
                    filter(id %in% qtl_collocated_raw_ids) %>%
                    nest_by(trait) %>%
                    mutate(trait=factor(trait,levels=trait_levels)) %>%
                    mutate(segcol=trait_rgb$tcol[as.integer(trait)]) %>%
                    unnest(data) %>%
                    ungroup()

    qmap_add <- qtl_raw %>%
                    mutate(group=paste0("lg",chr),
                           locus=trait,
                           position=as.numeric(position))

    qmap_sectcoldf <- qtl_raw %>%
                        mutate(chr=paste0("lg",chr),
                               s=as.numeric(interval_left),
                               e=as.numeric(interval_right),
                               col=rgb(t(col2rgb(segcol))/255,alpha=alpha))
    if(!is.null(sectcolindex)) {
        indexlevels = as.vector(unique(arrange(qtl_raw, pick(sectcolindex))[,sectcolindex]))[[sectcolindex]]
        qmap_add <- qmap_add %>%
                    mutate(index=factor(qtl_raw[[sectcolindex]], levels=indexlevels, ordered=TRUE)) %>%
                    select(group,position,locus,segcol,index) %>%
                    arrange(group,position)
        qmap_sectcoldf <- qmap_sectcoldf %>%
                                mutate(index=factor(qtl_raw[[sectcolindex]], levels=indexlevels, ordered=TRUE)) %>%
                                select(chr,s,e,col,index) %>%
                                arrange(chr,s)
    } else {
        qmap_add <- qmap_add %>%
                    select(group,position,locus,segcol) %>%
                    arrange(group,position)
        qmap_sectcoldf <- qmap_sectcoldf %>%
                                select(chr,s,e,col) %>%
                                arrange(chr,s)
    }

    qmap_sectcoldf <- data.frame(qmap_sectcoldf)

    qtldf <- qtl_collocated %>%
                rename(so="interval_left",
                       si="position",
                       eo="interval_right",
                       col="color") %>%
                mutate(ei=si,
                    qtl=trait,
                    chr=paste0("lg",chr)) %>%
                select(chr,qtl,so,si,ei,eo,col)

    lgs = unique(qtldf$chr)
    #Add linkage group lower/upper bounds so that whole linkage group is displayed
    qmap_add_lgbounds <- lgsums %>% 
                            filter(chr %in% lgs) %>%
                            pivot_longer(c("min","max")) %>%
                            rename(group=chr,
                                   locus=name,
                                   position=value) %>%
                            mutate(segcol=hsv(0,0,0,0)) %>%
                            select(group,position,locus,segcol)
    if(!is.null(sectcolindex)) {
        qmap_add_lgbounds$index = indexlevels[1]
    }
    qmap_add <- data.frame(qmap_add %>% bind_rows(qmap_add_lgbounds)) #Add in first and last segments of chromosomes so that the whole LG is shown
    if(!is.null(sectcolindex)) {
        qmap_add <- qmap_add %>% mutate(index=factor(index,levels=indexlevels, ordered=TRUE))
    }

    maxpos <- ceiling(max(all_qtl_raw$position))+10
    at.axis <- seq(0, maxpos)

    ## put labels on ruler at every 5 cM
    axlab <- vector()
    for (lab in 0:maxpos) {
        if (!lab %% 10) {
            axlab <- c(axlab, lab)
        }
        else {
            axlab <- c(axlab, NA)
        }
    }
    lmv.linkage.plot(qmap_add,
                    filename,
                    mapthese=lgs[sort(sapply(lgs,function(x) { as.numeric(regmatches(x,regexpr("[0-9]+$",x))) }),index.return=T)$ix],
                    ruler=TRUE,
                    lims=c(floor(min(lgsums$min))-10,ceiling(max(lgsums$max))+10),
                    lg.col=lg.col,
                    rsegcol=FALSE,
                    segcol="segcol",
                    qtldf=qtldf,
                    pdf.pointsize=8,
                    cex.axis = 1,
                    at.axis = at.axis, 
                    labels.axis = axlab,
                    sectcoldf=qmap_sectcoldf,
                    jitterloci=jitterloci,
                    lgw=lgw,
                    render.type=render.type,
                    ...
                    )
}

generate_linkagemapplots <- function(qtl_collocated, all_qtl_raw, filename, lgsums, lg.col="white", alpha=0.3, jitterloci=2, lgw=0.25, render.type='png', ...) {
    filename = paste0('../../Data/publication/figures/',filename,'.grouped.consensus.',render.type)
    generate_linkagemapplots_primitive(qtl_collocated, all_qtl_raw, filename, lgsums, lg.col, alpha, jitterloci, lgw, render.type, ...)
}

oldpar <- generate_linkagemapplots(colocated_cnj02_grouped, all_qtl, "cnj02_qtl_collated", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.2, lg.minheight=5.3, lgperrow=3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_cnj04_grouped, all_qtl, "cnj04_qtl_collated", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.2, lg.minheight=5.3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_cnj0x_grouped, all_qtl, "cnj0x_qtl_collated", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.25, lg.minheight=5.3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_diazGarciaImage_grouped, all_qtl, "diazGarcia2018ImagePhenotyping", genetic_map_lg_summary, lg.col="aliceblue", lg.minheight=5.3, pdf.width=5.5, lgperrow=3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_diazGarciaMassive_grouped, all_qtl, "diazGarcia2018MassivePhenotyping", genetic_map_lg_summary, jitterloci=1, lg.col="aliceblue", alpha=0.4, lg.minheight=5.3, pdf.width=5, lgperrow=4, noloci=T)
oldpar <- generate_linkagemapplots(colocated_diazGarciaImage_cnj0x_grouped, all_qtl, "diazGarcia2018ImagePhenotyping_cnj0x", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.2, lg.minheight=5.3, lgperrow=3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_diazGarciaMassive_cnj0x_grouped, all_qtl, "diazGarcia2018MassivePhenotyping_cnj0x", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.2, lg.minheight=5.3, noloci=T)
oldpar <- generate_linkagemapplots(colocated_schlautman_cnj0x_grouped, all_qtl, "schlautman_cnj0x", genetic_map_lg_summary, lg.col="aliceblue", alpha=0.1, lg.minheight=5.3, lgperrow=2, noloci=T)
oldpar <- generate_linkagemapplots(colocated_schlautman_diazGarciaImage_cnj0x_grouped, all_qtl, "schlautman_diazGarciaImagePhenotyping_cnj0x", lg.col="aliceblue", alpha=0.1, genetic_map_lg_summary, lg.minheight=5.3, lgperrow=4, noloci=T)

save.image(".RData.16_plotcolocatedqtllinkagemaps.R")

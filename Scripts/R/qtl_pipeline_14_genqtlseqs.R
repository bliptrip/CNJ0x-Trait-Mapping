#!/usr/bin/env RScript

#NOTE: This particular script generates informative tables for QTLs.
#

# loading libraries
library(seqinr)
library(tidyverse)
library(vcfR)

source('./usefulFunctions.R')

#Defaults (can be overridden with command-line invocation)
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
cnj02_workflow  <- get0("workflow", ifnotfound="../../Workflows/9")
cnj04_workflow  <- get0("workflow", ifnotfound="../../Workflows/10")
consensus_map   <- get0("consensus_map", ifnotfound="../../Data/genetic_data/RawData/consensusMapAll2.csv") #Consensus map loaded in order to find alternative markers in case something is not found -- uses markers in same bin
h2_threshold    <- get0("h2_threshold", ifnotfound=0.4)
qtl_variance_threshold    <- get0("qtl_variance_threshold", ifnotfound=8.0) #Only QTL with marker variance above this threshold will be used
genome            <- get0("genome", ifnotfound="../../Data/genetic_data/RawData/GCA_000775335.2_ASM77533v2_genomic.fna")
cnj0x_vcf_file          <- get0("cnj0x_vcf_file", ifnotfound="../../Data/genetic_data/RawData/cnj0x.all.mergedSNPs.vcf.gz")
gryg_vcf_file           <- get0("gryg_vcf_file", ifnotfound="../../Data/genetic_data/RawData/gryg.all.mergedSNPs.vcf.gz")
snp_seq_halfwidth       <- get0("snp_seq_halfwidth", ifnotfound=150)

#In case we override the workflow on the command-line
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

generateConsensusMapNearestMarker <- function(LG.filter, qtl, consensus.map, max_distance=5.0) {
    consensus.map.lg <- consensus.map %>% filter(LG == LG.filter)
    qtl.m <- cbind(as.matrix(qtl, ncol=1, nrow=length(qtl)), matrix(-1, ncol=1, nrow=length(qtl)))
    markers.pos <- rbind(matrix(1, nrow=1, ncol=nrow(consensus.map.lg)), matrix(consensus.map.lg$position, nrow=1, ncol=nrow(consensus.map.lg)))
    distances.m <- abs(qtl.m %*% markers.pos)
    distances.sorted.indices.m <- apply(distances.m, FUN=function(x) {xs = base::sort(x, index.return=T); return(xs$ix) }, MARGIN=1)
    closest_markers <- apply(distances.sorted.indices.m, 
                             FUN=function(i) 
                             {
                                markers.sorted <- consensus.map.lg[i, ]
                                markers.sorted.closest <- markers.sorted[1,] %>% 
                                    select(!LG) %>% 
                                    rename(marker_position='position')
                                return(tibble(markers.sorted.closest))
                             },
                             MARGIN=2)
    return(closest_markers)
}

genome.2014.fa <- read.fasta(file=genome,
                             seqtype="DNA",
                             whole.header=TRUE)
names(genome.2014.fa) <- gsub("^(.+)(scaffold_[[:digit:]]+)( [[:alpha:]]+)?(,.+)","\\2",names(genome.2014.fa), fixed=FALSE) #Replace names with only scaffold names

VGetSequence <- function(x) { if(is.null(x)) NA else getSequence(x, as.string=T) }
VGetAnnot <- function(x) { if(is.null(x)) NA else getAnnot(x) }

generateVariants <- function(vcf) {
    vcf.df <- as.tibble(vcf@fix) %>% 
                        mutate(marker=paste0("scaffold_",CHROM,"_",POS)) %>% 
                        rename(SNP_id="ID",
                            SNP_ref="REF",
                            SNP_alt="ALT",
                            SNP_qual="QUAL",
                            SNP_filter="FILTER",
                            SNP_info="INFO",
                            SNP_vcf_offset="POS") %>%
                        select(!CHROM)
    return(vcf.df)
}

getSNPSeqBounds <- function(lower, upper, strlength) {
    bounds <- list()
    if( lower < 1 ) {
        bounds$lower <- 1
        bounds$upper <- upper - lower + 1
        if( bounds$upper > strlength ) {
            bounds$upper = strlength
        }
    } else if( upper > strlength ) {
        bounds$upper <- strlength
        bounds$lower <- lower - (upper - strlength)
        if( bounds$lower < 1 ) {
            bounds$lower <- 1
        }
    } else {
        bounds$lower <- lower
        bounds$upper <- upper
    }
    return(bounds)
}

getSNPSeqLower <- function(lower, upper, strlength) {
    bounds <- getSNPSeqBounds(lower, upper, strlength)
    return(bounds$lower)
}

vGetSNPSeqLower <- Vectorize(getSNPSeqLower)

getSNPSeqUpper <- function(lower, upper, strlength) {
    bounds <- getSNPSeqBounds(lower, upper, strlength)
    return(bounds$upper)
}

vGetSNPSeqUpper <- Vectorize(getSNPSeqUpper)

getSNPSeq <- function(ssequence, soffset, slower, supper, sref, salt) {
    ssequence.split <- str_split(ssequence, "", n=Inf)
    for( i in 1:length(ssequence.split) ) {
        ssequence.split[[i]][soffset[i]] <- paste0("[",sref[i],"/",salt[i],"]")
        ssequence.split[[i]] <- ssequence.split[[i]][slower[i]:supper[i]]
    }
    ssequence.split <- lapply(ssequence.split, str_flatten)
    return(unlist(ssequence.split))
}

cnj0x_vcf <- read.vcfR(cnj0x_vcf_file)
cnj0x_vcf.df <- generateVariants(cnj0x_vcf)
str(cnj0x_vcf.df)
gryg_vcf <- read.vcfR(gryg_vcf_file)
gryg_vcf.df <- generateVariants(gryg_vcf)

#Read in consensus map, filter out scaffolds only found in genome, and then bin to select all markers in same bin
consensus.map <- read_csv(consensus_map)
consensus.map.binned <- consensus.map %>% 
                            select(marker,LG,consensus) %>%
                            filter(grepl("scaffold",marker)) %>% #Only look at marker bins that have scaffold_* designators, indicating they have positional info. in ref genome
                            mutate(scaffold = gsub("(scaffold_[^_]+)_.+", "\\1", marker, fixed=FALSE),
                                   SNP_offset = gsub("(scaffold_[^_]+)_(.+)", "\\2", marker, fixed=FALSE)) %>%
                            filter(!is.null(genome.2014.fa[scaffold])) %>% # Remove all markers not found in the genome
                            mutate(annot = sapply(genome.2014.fa[scaffold], VGetAnnot, simplify=T)) %>%
                            mutate(type = trim(gsub("^(.+)(scaffold_[[:digit:]]+)( [[:alpha:]]+)?(,.+)","\\3", annot, fixed=FALSE))) %>% #Exctract type, in case it is labeled in mitochondrial or chloroplast DNA -- don't want these!
                            filter(type == "") %>% #If no type designator, then it is nuclear DNA, and we only care about these scaffolds/markers.
                            mutate(id = gsub("^>(JOTO[[:alnum:].]+).+","\\1", annot, fixed=FALSE),
                                   scaffold_seq = as.character(sapply(genome.2014.fa[scaffold], VGetSequence, simplify=T))) %>%
                            rename(position='consensus') %>%
                            select(!type) %>%
                            select(!annot)
cnj0x.consensus.map.binned <- consensus.map.binned %>%
                                inner_join(cnj0x_vcf.df, by="marker") %>%
                                mutate(SNP_offset = as.numeric(SNP_offset)) %>%
                                mutate(SNP_lower = SNP_offset - snp_seq_halfwidth,
                                       SNP_upper = SNP_offset + snp_seq_halfwidth,
                                       scaffold_seq_len = str_length(scaffold_seq)) %>%
                                mutate(SNP_seq_lower = vGetSNPSeqLower(SNP_lower, SNP_upper, scaffold_seq_len),
                                       SNP_seq_upper = vGetSNPSeqUpper(SNP_lower, SNP_upper, scaffold_seq_len)) %>%
                                mutate(SNP_seq = getSNPSeq(scaffold_seq, SNP_offset, SNP_seq_lower, SNP_seq_upper, SNP_ref, SNP_alt))

gryg.consensus.map.binned <- consensus.map.binned %>%
                                inner_join(gryg_vcf.df, by="marker") %>%
                                mutate(SNP_offset = as.numeric(SNP_offset)) %>%
                                mutate(SNP_lower = SNP_offset - snp_seq_halfwidth,
                                       SNP_upper = SNP_offset + snp_seq_halfwidth,
                                       scaffold_seq_len = str_length(scaffold_seq)) %>%
                                mutate(SNP_seq_lower = vGetSNPSeqLower(SNP_lower, SNP_upper, scaffold_seq_len),
                                       SNP_seq_upper = vGetSNPSeqUpper(SNP_lower, SNP_upper, scaffold_seq_len)) %>%
                                mutate(SNP_seq = getSNPSeq(scaffold_seq, SNP_offset, SNP_seq_lower, SNP_seq_upper, SNP_ref, SNP_alt))

qtl.Schlautman.csv <- read_csv("../../Data/publication/tables/schlautman2015.qtls.csv",
                            col_names=c("trait","model","LG","nearest_marker","position","LOD","LOD2_min","LOD2_max","qtl_variance","MQ_Effect","CQ_Effect","Interaction_Effect"),
                            skip=1) %>%
                    mutate(population = "cnj02",
                           method = "MapQTLv6.0",
                           interval_2.0LOD = round((LOD2_max - LOD2_min), digits=2),
                           model_variance=NA) %>%
                    select(!LOD) %>%
                    select(!LOD2_max) %>%
                    select(!LOD2_min) %>%
                    select(!MQ_Effect) %>%
                    select(!CQ_Effect) %>%
                    select(!Interaction_Effect) %>%
                    filter(qtl_variance > qtl_variance_threshold) %>%
                    left_join(consensus.map %>% select(marker, consensus, LG) %>% rename(consensus_LG='LG', nearest_marker='marker'), by="nearest_marker") %>%
                    select(!LG) %>%
                    select(!position) %>%
                    rename(LG='consensus_LG',
                           position='consensus') %>%
                    mutate(map = NA) %>%
                    select(trait,model,LG,position,qtl_variance,model_variance,interval_2.0LOD,population,method,map) %>%
                    group_by(LG) %>%
                    mutate(consensus_data=generateConsensusMapNearestMarker(LG, position, cnj0x.consensus.map.binned)) %>%
                    ungroup() %>%
                    unnest(cols=consensus_data)
write_csv(file="../../Data/genetic_data/DerivedData/schlautman2015.qtls.csv", qtl.Schlautman.csv)


qtl.DG.GiNA.csv <- read_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.qtls.csv") %>%
                    select(!sequence_position) %>%
                    mutate(population = "gryg",
                           method = "stepwiseqtl",
                           Year = as.numeric(paste0("20",Year))) %>%
                    rename(trait = "Trait",
                           model = "Year",
                           LG = "Chromosome",
                           position = "Position (cM)",
                           nearest_marker = "Nearest marker",
                           qtl_variance = "Marker variance",
                           model_variance = "Model variance",
                           interval_1.5LOD = "1.5 LOD interval") %>%
                    filter(qtl_variance > qtl_variance_threshold) %>%
                    mutate(map="Cranberry-Composite_map-F1") %>%
                    select(trait,model,LG,position,qtl_variance,model_variance,interval_1.5LOD,population,method,map) %>%
                    group_by(LG) %>%
                    mutate(consensus_data=generateConsensusMapNearestMarker(LG, position, gryg.consensus.map.binned)) %>%
                    ungroup() %>%
                    unnest(cols=consensus_data)
write_csv(file="../../Data/genetic_data/DerivedData/diazGarcia2018ImagePhenotyping.supplemental.qtls.csv",qtl.DG.GiNA.csv)
                    
                    

qtl.DG.chemistry.csv <- read_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.qtls.csv",
                                 col_names=c("trait","LG","position","nearest_marker","lod","interval_1.5LOD","qtl_variance","model_variance","population","phenotyping_method"),
                                 skip=1) %>%
                    filter(phenotyping_method == "digital") %>%
                    mutate(model= "all-years",
                           method = "scanone") %>%
                    mutate(position = as.double(round(position, digits=2))) %>%
                    mutate(map="Cranberry-Composite_map-F1") %>%
                    select(trait,model,LG,position,qtl_variance,model_variance,interval_1.5LOD,population,method,map)

qtl.DG.chemistry.cnj0x.csv <- qtl.DG.chemistry.csv %>%
                    filter(population != "gryg") %>%
                    group_by(LG) %>%
                    mutate(consensus_data=generateConsensusMapNearestMarker(LG, position, cnj0x.consensus.map.binned)) %>%
                    ungroup() %>%
                    unnest(cols=consensus_data)

qtl.DG.chemistry.gryg.csv <- qtl.DG.chemistry.csv %>%
                    filter(population == "gryg") %>%
                    group_by(LG) %>%
                    mutate(consensus_data=generateConsensusMapNearestMarker(LG, position, gryg.consensus.map.binned)) %>%
                    ungroup() %>%
                    unnest(cols=consensus_data)

qtl.DG.chemistry.merged.csv <- bind_rows(qtl.DG.chemistry.cnj0x.csv,qtl.DG.chemistry.gryg.csv)
write_csv(file="../../Data/genetic_data/DerivedData/diazGarcia2018MassivePhenotyping.supplemental.qtls.csv",qtl.DG.chemistry.merged.csv)


generateQTLScaffolds <- function(workflow, map.binned, exclude) {
    qtl.h2.csv <- read_csv(paste0(cnj02_workflow, "/traits/h2.csv")) %>%
                            filter(model == "all-years") %>%
                            filter(h2 > h2_threshold) %>%
                            filter(!(trait %in% exclude))
    qtl.cnj0x.csv <- read_csv(paste0(workflow, "/traits/qtl_collated.csv"),
                                    col_names=c("method","model","trait","LG","position","chr2","position2","nearest_marker","qtl_lod","qtl_variance","qtl_pvalue","model_variance","interval_left","interval_marker_left","interval_1.5LOD","interval_right","interval_marker_right","GLRpvalue","GxYLRpvalue","GZRpvalue","GxYZRpvalue"),
                                    skip=1) %>%
                        filter(method == "stepwiseqtl") %>%
                        filter(model == "all-years") %>%
                        right_join(qtl.h2.csv, by=c("model","trait")) %>%
                        select(!trait) %>%
                        rename(trait="label_short") %>%
                        mutate(model= "all-years",
                            population="cnj02") %>%
                        group_by(model,trait) %>%
                        arrange(desc(qtl_variance)) %>%
                        filter(qtl_variance > 8.0) %>%
                        ungroup() %>%
                        mutate(map="Cranberry-Composite_map-F1") %>%
                        select(trait,model,LG,position,qtl_variance,model_variance,interval_1.5LOD,population,method,map) %>%
                        mutate(position = as.double(round(position, digits=2))) %>%
                        group_by(LG) %>%
                        mutate(consensus_data=generateConsensusMapNearestMarker(LG, position, map.binned)) %>%
                        ungroup() %>%
                        unnest(cols=consensus_data)
        return(qtl.cnj0x.csv)
}

exclude_traits <- c("BBI_UTBM","berry_skin","calyx_diam","calyx_lobe_form","calyx_lobe_size","chimera_solidity","chimera_tortuosity","chimera_umccLogX","chimera_umccLogY")
qtl.AM.cnj02.csv <- generateQTLScaffolds(cnj02_workflow, cnj0x.consensus.map.binned, exclude_traits)
write_csv(file="../../Data/genetic_data/DerivedData/maule2022Upright.cnj02.qtls.csv",qtl.AM.cnj02.csv)
qtl.AM.cnj04.csv <- generateQTLScaffolds(cnj04_workflow, cnj0x.consensus.map.binned, exclude_traits)
write_csv(file="../../Data/genetic_data/DerivedData/maule2022Upright.cnj04.qtls.csv",qtl.AM.cnj04.csv)


#Now convert and merge datasets for submission for SNP array
convertToSNPArrayFormat <- function(data) {
    data %>%
        rename(SNPid="SNP_id",
               REF="SNP_ref",
               ALT="SNP_alt",
               Chromosome="LG",
               Position="SNP_offset",
               Sequence="SNP_seq",
               Population="population",
               Trait="trait",
               Model="model",
               GeneticMap="map",
               QTLmethod="method",
               QTLposition="position",
               QTLvariance="qtl_variance",
               ModelVariance="model_variance",
               MarkerScaffold="scaffold",
               ScaffoldSequence="scaffold_seq",
               MarkerGeneticPosition="marker_position",
               GenomeRefId="id",
               SNPqual="SNP_qual",
               SNPfilter="SNP_filter",
               SNPinfo="SNP_info") %>%
        mutate(Genome="ASM77533v2") %>%
        select(SNPid,Chromosome,Position,REF,ALT,Sequence,Population,Trait,Model,GeneticMap,QTLmethod,QTLposition,QTLvariance,ModelVariance,MarkerScaffold,ScaffoldSequence,MarkerGeneticPosition,Genome,GenomeRefId,SNPqual,SNPfilter,SNPinfo)
}

qtl.combined.tsv <- rbind(convertToSNPArrayFormat(qtl.Schlautman.csv),
                          convertToSNPArrayFormat(qtl.DG.GiNA.csv),
                          convertToSNPArrayFormat(qtl.DG.chemistry.merged.csv),
                          convertToSNPArrayFormat(qtl.AM.cnj02.csv),
                          convertToSNPArrayFormat(qtl.AM.cnj04.csv))

write_tsv(file="../../Data/genetic_data/DerivedData/zalapa.qtl.tsv",qtl.combined.tsv)

save.image(paste0(workflow,"/.RData.14.genqtlseqs"))
load(paste0(workflow,"/.RData.14.genqtlseqs"))

# generador de huecolors
# 
hueGen<-function(n,from,to,s="100%",l="50%"){
    hsl<-paste0('hsl(',round(seq(from,to,length.out=n),0),',',s,',',l,')')
    return(hsl)
}


#Specify the directory prefix for storing data files.
DATA_REL_PATH <- get0("DATA_REL_PATH", ifnotfound='../..')
DATA_FOLDER_PREFIX <- get0("DATA_FOLDER_PREFIX", ifnotfound=paste0(DATA_REL_PATH,'/Data/phenotypic data/DerivedData/cleanup_data.R.output'))
DATA_ROBJS_FOLDER_PREFIX <- get0("DATA_ROBJS_FOLDER_PREFIX", ifnotfound=paste0(DATA_FOLDER_PREFIX,"/Robjs"))
DATA_PLOTS_FOLDER_PREFIX <- get0("DATA_PLOTS_FOLDER_PREFIX", ifnotfound=paste0(DATA_FOLDER_PREFIX,"/plots"))
GDATA_FOLDER_PREFIX <- get0("GDATA_FOLDER_PREFIX", ifnotfound=paste0(DATA_REL_PATH,"/Data/genetic_data/RawData"))
GDDATA_FOLDER_PREFIX <- get0("GDDATA_FOLDER_PREFIX", ifnotfound=paste0(DATA_REL_PATH,"/Data/genetic_data/DerivedData"))
PDATA_FOLDER_PREFIX <- get0("PDATA_FOLDER_PREFIX", ifnotfound=paste0(DATA_REL_PATH,"/Data/phenotypic data/RawData"))
PDDATA_FOLDER_PREFIX <- get0("PDDATA_FOLDER_PREFIX", ifnotfound=paste0(DATA_REL_PATH,"/Data/phenotypic data/DerivedData"))

# Returns string without leading white space
trim.leading <- function (x)  sub("^\\s+", "", x)

# Returns string without trailing white space
trim.trailing <- function (x) sub("\\s+$", "", x)

# Returns string without leading or trailing white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#Wrapper functions for readRDS() and writeRDS()
saveRDSw <- function(object, file, compress=FALSE) {
    saveRDS(object, paste0(DATA_ROBJS_FOLDER_PREFIX,"/",file), compress=compress)
}

readRDSw <- function(file) {
    return(readRDS(paste0(DATA_ROBJS_FOLDER_PREFIX,"/",file)))
}

quartzw <- function(file, title) {
    quartz(title=title, file=paste0(DATA_PLOTS_FOLDER_PREFIX,"/",file), type="pdf")
}

postscriptw <- function(file, title, width=10, height=10) {
    postscript(title=title, file=paste0(DATA_PLOTS_FOLDER_PREFIX,"/",file), width=width, height=height, horizontal=F)
}

write.csvw <- function(df, file) {
    write.csv(df,paste0(DATA_FOLDER_PREFIX,"/",file))
}

read.csvw <- function(file, row.names=NULL) {
    return(read.csv(paste0(DATA_FOLDER_PREFIX,"/",file), row.names=row.names))
}

geno_rpath2fpath <- function(rpath) {
    return(paste0(GDATA_FOLDER_PREFIX, "/", rpath))
}

geno_dpath2fpath <- function(rpath) {
    return(paste0(GDDATA_FOLDER_PREFIX, "/", rpath))
}

pheno_rpath2fpath <- function(rpath) {
    return(paste0(PDATA_FOLDER_PREFIX, "/", rpath))
}

pheno_dpath2fpath <- function(rpath) {
    return(paste0(PDDATA_FOLDER_PREFIX, "/", rpath))
}

is.empty <- function(mstr) {
    return(is.na(mstr) || (trimws(mstr,which="both") == ""))
}

trait_is_unmasked <- function(trait.cfg) {
    return(is.na(trait.cfg$mask) || (trait.cfg$mask != "TRUE"))
}

mtrait_subfolder <- function(trait.cfg) {
    model                    <- as.character(trait.cfg$model)
    year                     <- as.numeric(trait.cfg$year)
    traits                   <- unlist(strsplit(trait.cfg$mtraits,","))
    trait.names              <- paste0(traits,collapse="__")
    trait_subfolder          <- paste0(c(trait.cfg$model,trait.names),collapse="--")
    return(trait_subfolder)
}


generate_cross_file <- function(trait.cfg, blups, gData, superMap.df, file.path) {
    print(paste0("Population Analysis Set: ", trait.cfg["model"]))
    #First make sure to use individuals only in common in the blups data and in the marker data
    geno.intersect<-intersect(names(blups),rownames(gData))
    y<-as.data.frame(blups[geno.intersect]) #Convert to data.frame to deal with issues in setting colnames in univariate analysis (blups aren't a data.frame in this case)
    colnames(y)<-trait.cfg$trait
    gData.sub <- gData[geno.intersect,]

    #Second, we are only interested in markers in common with the consensus map
    geno.intersect.sub<-intersect(colnames(gData.sub),superMap.df$marker)
    gData.sub<-gData.sub[,geno.intersect.sub]
    superMap.sub<-superMap.df[geno.intersect.sub,]

    gData.sub<-rbind(superMap.sub$LG,superMap.sub$consensus,gData.sub)
    colnames(gData.sub)<-superMap.sub$marker

    #Invert the gData.sub, bin the data, re-invert, and then write cross file.
    gData.inv.sub <- t(gData.sub)
    colnames(gData.inv.sub) <- c("LG","consensus",colnames(gData.inv.sub)[3:ncol(gData.inv.sub)])
    #Bin the marker data
    gData.bin.df  <- condense_map_bins(data.frame(gData.inv.sub))
    #Subset the genotype marker calls based on the binned results
    gData.inv.sub <- gData.inv.sub[rownames(gData.bin.df),]
    gData.sub     <- t(gData.inv.sub)

    gData.sub<-data.frame(rbind('','',y),gData.sub)
    rownames(gData.sub)[1:2]<-c('chr','pos')

    #Now write the cross file as a csv in the appropriate folder
    write.csv(gData.sub,file=paste0(file.path,"/cross.csv"),row.names=FALSE)
    cross <- read.cross(format = "csv", file=paste0(file.path,"/cross.csv"), genotypes = NULL)
    cross <- jittermap(cross)
    cross <- calc.genoprob(cross,step=0,map.function="kosambi")
    saveRDS(cross, file=paste0(file.path,"/cross.rds"), compress=TRUE)
}

trait.abbrev.map.df <- read.csv(pheno_dpath2fpath("trait_to_abbreviation_map.csv"), header=F, row.names=1, stringsAsFactors=F)
rename_traits <- function(traits, abbrevmap) {
    #If you don't convert the trait to a character (if it's a factor), this function will fail to correctly map!
    traits  <- as.character(traits)
    abbrevs <- vector("character", length(traits))
    for( i in 1:length(traits) ) {
        trait        <- traits[i]
        if( trait %in% rownames(abbrevmap) ) {
            abbrevs[i]   <- abbrevmap[trait, 1]
        } else {
            abbrevs[i]   <- trait
        }
    }
    return(abbrevs)
}

circosfile2path <- function(filename) {
		#For now, this is the relative path to the location of the interactive circos qtl html file/javascript set
		return(paste0("configs/circos/",filename))
}


model_to_name <- function(trait.cfg.tb,model,trait) {
    model_names <- vector(mode="character", length=length(model))
    for( i in 1:length(model) ) {
        fstring <- paste0("trait.cfg.tb %>% filter((model == \"",model[i],"\") & (trait == \"",trait[i],"\"))")
        mt <- eval(parse(text=fstring))
        model_names[i] <- mt$model_label[1]
    }
    return(model_names)
}

trait_to_name <- function(trait.cfg.tb,model,trait) {
    trait_names <- vector(mode="character", length=length(model))
    for( i in 1:length(model) ) {
        fstring <- paste0("trait.cfg.tb %>% filter((model == \"",model[i],"\") & (trait == \"",trait[i],"\"))")
        mt <- eval(parse(text=fstring))
        trait_names[i] <- mt$label[1]
    }
    return(trait_names)
}

siginfo <- function(x){
    y="NS"
    if(!is.na(x)) {
        if(x >= 0 & x < 0.001){y="***"}
        else if(x >= 0.001 & x < 0.01){y="**"}
        else if(x >= 0.01 & x < 0.05){y="*"}
        else if(x >= 0.05 & x < 0.1){y="."}
    }
    return(y)
}

round.digits <- function(values,digits) {
	rvalue <- gsub("NA", "", format(round(values,digits=digits),nsmall=digits))
}

signif.digits <- function(values,digits) {
	rvalue <- as.numeric(gsub("NA", "", format(round(values,digits=digits),nsmall=digits)))
}


									
#Function for looping through all unmasked traits in the configs/model-traits.cfg.csv file.


loopThruTraits <- function(workflow, loopFunCallback, loopArgs=NULL) {
        #Loop over all mmers, trait groups, subgroups, and perform makeqtl() and fitqtl().
        traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
        for( i in 1:length(traits.df[,1]) ) {
            trait.cfg       <- traits.df[i,]
            if ( trait_is_unmasked(trait.cfg) ) {
                    model <- trait.cfg$model
                    trait <- trait.cfg$trait
                    trait_subfolder          <- paste0(c(model,trait),collapse="--")
                    trait_subfolder_fpath    <- file.path(paste0(workflow,"/traits"), trait_subfolder)
                    loopFunCallback(trait.cfg, trait_subfolder_fpath, loopArgs)
            }
        }
}

#Define a poiner object that uses the environment to pass variables by reference.  It uses generic S3 methods to
#generate the constructor with the envrionment and more.
newPointer<-function(inputValue){ 
    object=new.env(parent=globalenv()) 
    object$value=inputValue 
    class(object)='pointer'

    return(object) 
} 

copy<-function (object, ...) { # create S3 generic 
    UseMethod("copy")  
}  
copy.pointer<-function(object1,object2=NULL,...){ 
    if (is.null(object2)) {  
        object2 <- new.env(parent = globalenv())  
            class(object2) <- class(object1) 
            nullFlag <- TRUE  
    }  
    elements <- names(object1)  
        for (index in 1:length(elements)) {  
            assign(elements[index], get(elements[index], env = object1,  
                        inherits = FALSE), env = object2)  
        }  
    if (nullFlag)  
    { return(object2)  
    } else {  
        return(NULL)  
    }  
} 

updatePointerValue<-function (object, ...) { # create S3 generic 
    UseMethod("updatePointerValue")  
}  
updatePointerValue.pointer<-function(object,newValue){ # create S3 method 
    if (!is(object, "pointer")) { stop(" 'object' argument must be of class 'pointer' .") }  
    object$value<-newValue 
        return(NULL) 
} 

append.pointer<-function(object, value) {
    if (!is(object, "pointer")) { stop(" 'object' argument must be of class 'pointer' .") }
    object.class <- class(object$value)
    if( object.class == "data.frame" ) {
        pos <- nrow(object$value) + 1
        object$value[pos,] <- value
    } else if( object.class == "list" ) {
        pos <- length(object$value) + 1
        object$value[[pos]] <- value
    }
    return(object)
}

assign.pointer<-function(object, value, pos) {
    if (!is(object, "pointer")) { stop(" 'object' argument must be of class 'pointer' .") }  
    object.class <- class(object$value)
    if( object.class == "data.frame" ) {
        object$value[pos,] <- value
    } else if( object.class == "list" ) {
        object$value[[pos]] <- value
    }
    return(object)
}

#The following functions can only be used within knitr renderings
is_word_output = function() {
    return(grepl("docx",knitr::opts_knit$get("rmarkdown.pandoc.to")))
}

is_pdf_output = function() {
    return(grepl("latex",knitr::opts_knit$get("rmarkdown.pandoc.to")))
}

is_html_output = function() {
    return(grepl("html",knitr::opts_knit$get("rmarkdown.pandoc.to")))
}

#Shamelessly taken from http://www.sthda.com/english/wiki/correlation-matrix-formatting-and-visualization
flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
    )
}

# # markers
#                   '<plot>',
#                   'type=text',
#                   'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/markers.txt',
#                   'color=black',
#                   'r1=0.99r',
#                   'r0=.80r',
#                   'show_links=yes',
#                   'link_dims=8p,8p,8p,4p,8p',
#                   'link_thickness=2p',
#                   'link_color=black',
#                   'label_size=13p',
#                   'label_font=default',
#                   'padding=2p',
#                   'rpadding=2p',
#                   'max_snuggle_distance=3r',
#                   'label_snuggle=yes',
#                   '</plot>',




# generador de huecolors
# 
hueGen<-function(n,from,to){
    hsl<-paste0('hsl(',round(seq(from,to,length.out=n),0),',100%,50%',')')
    return(hsl)
}


#Specify the directory prefix for storing data files.
DATA_FOLDER_PREFIX <- "../../Data/phenotypic data/DerivedData/cleanup_data.R.output"
DATA_ROBJS_FOLDER_PREFIX <- paste0(DATA_FOLDER_PREFIX,"/Robjs")
DATA_PLOTS_FOLDER_PREFIX <- paste0(DATA_FOLDER_PREFIX,"/plots")
GDATA_FOLDER_PREFIX <- "../../Data/genetic_data/RawData"
GDDATA_FOLDER_PREFIX <- "../../Data/genetic_data/DerivedData"
PDATA_FOLDER_PREFIX <- "../../Data/phenotypic data/RawData"
PDDATA_FOLDER_PREFIX <- "../../Data/phenotypic data/DerivedData"

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

circosfile2path <- function(filename) {
		#For now, this is the relative path to the location of the interactive circos qtl html file/javascript set
		return(paste0("../../../configs/circos/",filename))
}

#Function for looping through all unmasked traits in the configs/model-traits.cfg.csv file.
loopThruTraits <- function(workflow, loopFunCallback, loopArgs=NUL) {
        #Loop over all mmers, trait groups, subgroups, and perform makeqtl() and fitqtl().
        traits.df <- read.csv(file=paste0(workflow,"/configs/model-traits.cfg.csv"),header=T,stringsAsFactors=F)
        for( i in 1:length(traits.df[,1]) ) {
            trait.cfg       <- traits.df[i,]
            if ( trait_is_unmasked(trait.cfg) ) {
                    model                    <- as.character(trait.cfg$model)
                    year                     <- as.numeric(trait.cfg$year)
                    traits                   <- unlist(strsplit(trait.cfg$mtraits,","))
                    trait.names              <- paste0(traits,collapse="__")
                    trait_subfolder          <- paste0(c(trait.cfg$model,trait.names),collapse="--")
                    trait_subfolder_fpath    <- file.path(paste0(workflow,"/traits"), trait_subfolder)
                    #Read in the model result file
                    traits  <- unlist(strsplit(trait.cfg$mtraits,","))
                    loopFunCallback(trait.cfg, trait.names, traits, trait_subfolder_fpath, loopArgs)
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
    pos <- nrow(object$value) + 1
    object$value[pos,] <- value
    return(object)
}

assign.pointer<-function(object, value, pos) {
    if (!is(object, "pointer")) { stop(" 'object' argument must be of class 'pointer' .") }  
    object$value[pos,] <- value
    return(object)
}

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

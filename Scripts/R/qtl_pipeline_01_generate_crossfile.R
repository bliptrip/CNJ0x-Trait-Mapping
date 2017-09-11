#This is a QTL pipeline script given to me by Luis, but I've adapted to work with the Vorsa upright datasets from 2011-2014
#

# loading libraries
source('./usefulFunctions.R')


library(lme4)
library(qtl)
library(sommer)
library(IRanges)
library(GenomicRanges) 
library(intervals)
library(snow)
library(doSNOW)
library(lattice)
library(RColorBrewer)

require(openxlsx)

#Read the mixed model analysis results from the R-data file for later analysis -- Start from here if we want to save time
cnjpop.mmer.p2 <- readRDSw('cnjpop.mmer2.p2.rds')


#Convert consensus map to bins
superMap.df<-read.table(geno_rpath2fpath('consensusMapAll2.csv'),header=T,sep=',')
superMap.df<-superMap.df[,c('marker','LG','consensus')]
superMap.df<-superMap.df[order(superMap.df[,2],superMap.df[,3]),]
superMap.df$binID<-NA
superMap.bin.df<-numeric()
for (LG in unique(superMap.df$LG)){
  f<-which(superMap.df$LG==LG)
  mybins<-unique(superMap.df$consensus[f])
  mybinsID<-paste0('bin_',LG,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap.df$consensus[f]==mybins[j])
    superMap.df$binID[f[f2]]<-mybinsID[j]      
    superMap.bin.df<-rbind(superMap.bin.df,superMap.df[f[f2[1]],])
  }
}
rownames(superMap.bin.df)<-superMap.bin.df$marker


## genetic analysis
geno<-read.table(geno_rpath2fpath('CNJ02_AllASMapData.csv'),header=T,sep=',')  


matrixK<-matrix(NA,nrow=nrow(geno),ncol=ncol(geno)-6)
g<-which(geno[,7:ncol(geno)]=='ac')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='ad')
matrixK[g]<-3  
g<-which(geno[,7:ncol(geno)]=='bc')
matrixK[g]<-2
g<-which(geno[,7:ncol(geno)]=='bd')
matrixK[g]<-4
g<-which(geno[,7:ncol(geno)]=='ll')
matrixK[g]<-5
g<-which(geno[,7:ncol(geno)]=='lm')
matrixK[g]<-6
g<-which(geno[,7:ncol(geno)]=='nn')
matrixK[g]<-7
g<-which(geno[,7:ncol(geno)]=='np')
matrixK[g]<-8
g<-which(geno[,7:ncol(geno)]=='hh')
matrixK[g]<-1
g<-which(geno[,7:ncol(geno)]=='hk')
matrixK[g]<-10
g<-which(geno[,7:ncol(geno)]=='kk')
matrixK[g]<-4

matrixK[1:10,1:10]
gData<-t(matrixK)
rownames(gData)<-colnames(geno)[7:ncol(geno)]
colnames(gData)<-geno$X
gData[1:10,1:10]


#Generate the QTL cross files for reading into R/qtl
for (mmer in cnjpop.mmer.p2) {
    print(paste0("Population Analysis Set: ",mmer$description))
    for(analysis in mmer$analyses) {
        traits <- analysis$traits
        geno.intersect<-intersect(rownames(analysis$blups),rownames(gData))
        y<-as.data.frame(analysis$blups[geno.intersect,traits]) #Convert to data.frame to deal with issues in setting colnames in univariate analysis (blups aren't a data.frame in this case)
        colnames(y)<-traits

        gData.sub <- gData[geno.intersect,]
        geno.intersect.sub<-intersect(colnames(gData.sub),superMap.bin.df$marker)
        gData.sub<-gData.sub[,geno.intersect.sub]
        superMap.sub<-superMap.bin.df[geno.intersect.sub,]

        gData.sub<-rbind(superMap.sub$LG,superMap.sub$consensus,gData.sub)

        colnames(gData.sub)<-superMap.sub$binID

        gData.sub<-data.frame(rbind('','',y),gData.sub)
        rownames(gData.sub)[1:2]<-c('chr','pos')

        qtlfile=paste0('qtl/',mmer$description,"__",paste(traits,collapse="__"),"__QTL.csv")
        write.csv(file=geno_dpath2fpath(qtlfile),gData.sub,row.names = FALSE)
    }
}

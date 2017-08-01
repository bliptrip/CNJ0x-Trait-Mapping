#Only install the following packages when running for the first time
install.packages(c("lme4"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("qtl"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("sommer"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#IRanges is a part of bioconductor package: 
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")
biocLite("GenomicRanges")

install.packages(c("intervals"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
#The following may need to be run as a superuser for the first time on a macosx.  Also, one may manually need to install openmpi.
install.packages(c("Rmpi"), repos = "http://mirror.las.iastate.edu/CRAN/", configure.args="--with-Rmpi-include=/opt/openmpi/include --with-Rmpi-libpath=/opt/openmpi/lib --with-Rmpi-type=OPENMPI", dependencies=TRUE, verbose=TRUE)
install.packages(c("snow"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)
install.packages(c("doSNOW"), repos = "http://mirror.las.iastate.edu/CRAN/", dependencies = TRUE)

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

# creating a list to save phenotypes, r2 and other important data

#AM - Unecessary.  Quoting out.
#phenoData<-list()


# reading CNJ02 and getting BLUPS

jagCN02<-read.csv('chemicalData/CNJ02_data.csv',header=T,sep=',')
#jagCN02<-read.csv("/crdata4/luis4/QTLcolorServer_copy/QTLcolorPaper/chemicalData/CNJ02_data.csv",header=T,sep=',')


# removing outliers
#plot(jagCN02$mfw)
f<-which(jagCN02$mfw<.5)
jagCN02<-jagCN02[-f,]
#plot(jagCN02$brix)
f<-which(jagCN02$brix>20)
jagCN02<-jagCN02[-f,]
#plot(jagCN02$ta)
f<-which(jagCN02$ta>5)
jagCN02<-jagCN02[-f,]
#plot(jagCN02$pac)
#plot(jagCN02$tacy)




# getting everything by separate
## separating data by month

jag<-list()
mymonths<-unique(jagCN02$Month)
myyears<-unique(jagCN02$Year)
mytraits<-c('mfw','tacy','brix','ta','pac')
#6 different sampling times: 3 years, 2 months/year = 6 times
completeGenotypes<-names(which(table(jagCN02$Genotype)==6))

cont<-1
for (i in 1:2){
  for (j in 1:3){
    f<-which(jagCN02$Year==myyears[[j]] & jagCN02$Month==mymonths[[i]] & (jagCN02$Genotype %in% completeGenotypes))
    jag[[cont]]<-jagCN02[f,c(1,4,5,6,7,8)]
    cont<-cont+1
  }
}





phenoData<-list()



geno<-read.table('genetic_data/RawData/CNJ02_AllASMapData.csv',header=T,sep=',')  
rownames(geno)<-geno$X
geno<-geno[,-c(1:6)]
geno2<-atcg1234(t(geno))



# there are problems in 
# AM - Note: Trait 3 is brix, trait 5 is PAC
# trait 3, sep, 2013
# trait 5, sep, 2011
# trait 3, oct, 2011 and 2013
# trait 5, oct, 2011



for (i in 1:5){
    #September
    phenoA<-cbind(jag[[1]][,i+1],jag[[2]][,i+1],jag[[3]][,i+1])
    #October
    phenoB<-cbind(jag[[4]][,i+1],jag[[5]][,i+1],jag[[6]][,i+1])

    #Brix - Remove 2013 data for September, 2011 and 2013 for October
    if (i == 3){
        phenoA<-phenoA[,-3]
        phenoB<-phenoB[,-c(1,3)]
    }

    #AM - PAC - Remove 2011 data for both months
    if (i == 5){
        phenoA<-phenoA[,-1]
        phenoB<-phenoB[,-1]
    }


    if (i==3){
        rownames(phenoA)<-jag[[1]]$Genotype
        #AM - Why not rownames?
        names(phenoB)<-jag[[1]]$Genotype
    } else{
        rownames(phenoA)<-jag[[1]]$Genotype
        rownames(phenoB)<-jag[[1]]$Genotype
    }

    f<-intersect(rownames(phenoA),rownames(geno2))
    phenoA2<-phenoA[f,]
    geno3<-geno2[f,]
    A<-A.mat(geno3)


    Za <- diag(nrow(A)) 
    ETA.A <- list(add=list(Z=Za,K=A)) 
    repo1<-mmer(Y=phenoA2,Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=TRUE)
    blups1<-repo1$u.hat$add
    var1<-repo1$var.comp

    if (i==3) {
        f<-intersect(names(phenoB),rownames(geno2))
        phenoA2<-phenoB[f]
    } else {
        f<-intersect(rownames(phenoB),rownames(geno2))
        phenoA2<-phenoB[f,]
    }
    geno3<-geno2[f,]
    A<-A.mat(geno3)


    Za <- diag(nrow(A)) 
    ETA.A <- list(add=list(Z=Za,K=A)) 
    repo1<-mmer(Y=phenoA2,Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=TRUE)
    blups2<-repo1$u.hat$add
    var2<-repo1$var.comp

    phenoData[[i]]<-list(sep=list(blups=as.matrix(blups1),vcov=var1,trait=mytraits[i]),oct=list(blups=as.matrix(blups2),vcov=var2,trait=mytraits[i]))

}

# calculating h2 and correlations
for (i in 1:5){
    for (j in 1:2){
        #For brix and october, there is no correlation matrix calculated.
        if (i==3 & j==2) {
            #AM edit: Had to add 'component' reference as was getting error in sum function.
            phenoData[[i]][[j]]$h2=phenoData[[i]][[j]]$vcov[1,'component']/sum(phenoData[[i]][[j]]$vcov[,'component'])
            names(phenoData[[i]][[j]]$h2)<-myyears[i]
            phenoData[[i]][[j]]$gcor=NA
            names(phenoData[[i]][[j]]$gcor)<-myyears[i]
            names(phenoData[[i]][[j]]$blups)<-rownames(phenoA2)
        } else {
            gv<-diag(phenoData[[i]][[j]]$vcov$add)
            rv<-diag(phenoData[[i]][[j]]$vcov$Residual)
            phenoData[[i]][[j]]$h2<-gv/(gv+rv)

            x<-cor(phenoData[[i]][[j]]$blups,use='complete.obs')
            goodY<-as.numeric(unlist(regmatches(colnames(x),gregexpr('[0-9]+',colnames(x)))))
            phenoData[[i]][[j]]$gcor<-x
            colnames(phenoData[[i]][[j]]$gcor)<-myyears[goodY]
            rownames(phenoData[[i]][[j]]$gcor)<-myyears[goodY]
            names(phenoData[[i]][[j]]$h2)<-myyears[goodY]
            rownames(phenoData[[i]][[j]]$blups)<-rownames(phenoA2)
        }
    }
}

superMap<-read.table('genetic_data/RawData/consensusMapAll2.csv',header=T,sep=',')
superMap<-superMap[,c('marker','LG','consensus')]
superMap<-superMap[order(superMap[,2],superMap[,3]),]
superMap$binID<-NA
superMap2<-numeric()
for (i in 1:12){
  f<-which(superMap$LG==i)
  mybins<-unique(superMap$consensus[f])
  mybinsID<-paste0('bin_',i,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap$consensus[f]==mybins[j])
    superMap$binID[f[f2]]<-mybinsID[j]      
    superMap2<-rbind(superMap2,superMap[f[f2[1]],])
  }
}

rownames(superMap2)<-superMap2$marker

## genetic analysis
geno<-read.table('genetic_data/RawData/CNJ02_AllASMapData.csv',header=T,sep=',')  


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



for (ii in 1:5){
    for (i in 1:2){
        if (ii==3 & i==2) {
            f<-intersect(names(phenoData[[ii]][[i]]$blups),rownames(gData))
            y2<-phenoData[[ii]][[i]]$blups[f]
            names(y2)<-myyears[ii]
        } else {
            goodY<-as.numeric(unlist(regmatches(colnames(phenoData[[ii]][[i]]$blups),gregexpr('[0-9]+',colnames(phenoData[[ii]][[i]]$blups)))))
            f<-intersect(rownames(phenoData[[ii]][[i]]$blups),rownames(gData))
            y2<-phenoData[[ii]][[i]]$blups[f,]
            colnames(y2)<-myyears[goodY]
        }
        gData2<-gData[f,]
        gData2[1:10,1:10]

        f<-intersect(colnames(gData2),superMap2$marker)

        gData3<-gData2[,f]
        superMap3<-superMap2[f,]

        gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
        gData4[1:10,1:10]

        colnames(gData4)<-superMap3$binID

        if (ii==3 & i==2){
            gData5<-data.frame(c('','',y2),gData4)
        }else{
            gData5<-data.frame(rbind('','',y2),gData4)
        }
        rownames(gData5)[1:2]<-c('chr','pos')
        gData5[1:10,1:20]

        write.csv(file='temp/testQTL.csv',gData5,row.names = FALSE)

        superF <- (read.cross(format = "csv", file='temp/testQTL.csv',genotypes = NULL))

        phenoData[[ii]][[i]]$cross<-numeric()
        phenoData[[ii]][[i]]$cross <- calc.genoprob(superF,step=0,map.function="kosambi") 

    }
}


sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties,additive.only=FALSE)
  return(x1)
}

#For vaccininum server, set cores to 32
n.cores <- 32

cl <- snow::makeCluster(n.cores, spec="MPI")
registerDoSNOW(cl)
for (i in 1:5){
    for (j in 1:2){
        superOut <- foreach(ph = 1:ncol(phenoData[[i]][[j]]$cross$pheno),.packages = 'qtl') %dopar% perroQTL(phenoData[[i]][[j]]$cross,ph,max.qtl=10,penalties=sixo)
        phenoData[[i]][[j]]$scantwo<-superOut
    }
    print(i)
}

snow::stopCluster(cl)

for (tri in 1:5){
    for (ii in 1:2){
        modelRes<-vector("list", 4) 
        for (tri2 in 1:ncol(phenoData[[tri]][[ii]]$cross$pheno)){
            rob1<-phenoData[[tri]][[ii]]$scantwo[[tri2]]
            if (length(rob1)==0){

            }else{
                superQTL<-makeqtl(phenoData[[tri]][[ii]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
                md1<-fitqtl(phenoData[[tri]][[ii]]$cross,pheno.col = tri2,superQTL,formula = formula(rob1),get.ests=F)
                x<-summary(md1)
                modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                        drop=as.data.frame(x$result.drop),
                        model=rob1,
                        inter=numeric())
                for (j in 1:length(rob1$chr)){
                    mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
                    modelRes[[tri2]]$inter[j]<-mylod
                }
            }
        } 
        phenoData[[tri]][[ii]]$qtldata<-modelRes
    }
}

# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

jaguar1<-data.frame(trait=numeric(),year=numeric(),date=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:5){
  for (i2 in 1:2){
  t1<-mytraits[i]
  t2<-as.character(mymonths[i2])
  for (j in 1:ncol(phenoData[[i]][[i2]]$cross$pheno)){
    t3<-as.character(myyears[j])
    nQTL<-length(phenoData[[i]][[i2]]$qtldata[[j]]$model$chr)
    if (nQTL==0){

      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]][[i2]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]][[i2]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]][[i2]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]][[i2]]$qtldata[[j]]$drop))
      if (nQTL==1){
        xK<-phenoData[[i]][[i2]]$qtldata[[j]]$variance
      }else{
      x2<-phenoData[[i]][[i2]]$qtldata[[j]]$drop[f2,]
      xK<-x2$'%var'
    }

      jaguar1[cont,]<-c(t1,t3,t2,mychr,mypos,mymarker,xK,phenoData[[i]][[i2]]$qtldata[[j]]$variance,phenoData[[i]][[i2]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)

jaguar1$date2<-paste0(jaguar1$date,jaguar1$year)
mydates<-as.character(unique(jaguar1$date2))


# circos


# ideogram
ideo1<-numeric()
for (i in 1:12){
  ideo1<-c(ideo1,paste0('chr - c',i,' chr',i,' 0 ',max(superMap3$consensus[which(superMap3$LG==i)])*1000," vvdgrey"))
}
write(ideo1, file = '/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogram.txt')




mycols<-matrix(paste0(hueGen(30,0,340),'_a2'),nrow=6,byrow=F)
mycols2<-matrix(hueGen(30,0,340),nrow=6,byrow=F)


for (ii in 1:5){

mark1<-matrix(NA,nrow=length(which(jaguar1$trait==mytraits[ii])),ncol=1)
cont<-1
for (i in 1:nrow(jaguar1)){
  if (jaguar1$trait[i]==mytraits[ii]){
   mark1[cont,]<-paste0('c',jaguar1$chr[i]," ",round(jaguar1$position[i]*1000,0),' ',round(jaguar1$position[i]*1000,0)+1,' ',which(jaguar1$date2[i]==mydates),' color=',mycols[which(jaguar1$date2[i]==mydates),ii],',stroke_color=',mycols2[which(jaguar1$date2[i]==mydates),ii],',glyph_size=',round(jaguar1$marker.variance[i]*3,0))
  cont<-cont+1
}
}

filename<-paste0('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat',ii,'.txt')
write(mark1, file = filename)

mark2<-matrix(NA,nrow=length(which(jaguar1$trait==mytraits[ii])),ncol=1)

cont<-1
for (i in 1:nrow(jaguar1)){
  if (jaguar1$trait[i]==mytraits[ii]){
   minInt<-round(jaguar1$position[i]+(1000*(jaguar1$position[i]-.5*jaguar1$interval[i])),0)
   maxInt<-round(jaguar1$position[i]+(1000*(jaguar1$position[i]+.5*jaguar1$interval[i])),0)
   mark2[cont,]<-paste0('c',jaguar1$chr[i]," ",minInt,' ',maxInt,' ',which(jaguar1$date2[i]==mydates),' color=',mycols2[which(jaguar1$date2[i]==mydates),ii])
   cont<-cont+1
}
}


filename<-paste0('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line',ii,'.txt')
write(mark2, file = filename)


}


back1<-c('<backgrounds>',
          '<background>',
          'color=vdgrey',
          'y0=0',
          '</background>',
          '<background>',
          'y0=3.5',
          'y1=7',
          'color=dgrey',
          '</background>',
          '</backgrounds>')

back2<-c('<backgrounds>',
          '<background>',
          'y0=0',
          'y1=7',
          'color=vdgrey',
          '</background>',
          '</backgrounds>')



rmin<-seq(.3,.85,length.out=5)
rmax<-rmin-(.01-diff(rmin[1:2]))

rmin1<-rmin
rmax1<-rmin+.7*diff(rmin[1:2])

rmin2<-rmax1+0.01
rmax2<-rmax


rmin1<-round(rmin1,3)
rmin2<-round(rmin2,3)
rmax1<-round(rmax1,3)
rmax1<-round(rmax1,3)

circosStarter<-c('<<include etc/colors_fonts_patterns.conf>>',
                  '<<include /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ticks.conf>>',
                  'karyotype = /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogram.txt',
                  'chromosomes_units=1000',
                  'chromosomes_display_default=yes',
                  '<<include /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/ideogramLalo.conf>>',
                  '<image>',
                  'radius* = 1500p',
                  '<<include etc/image.conf>>',
                  '</image>',
                  'chromosomes_radius=c1:1r,c2:1r,c3:1r,c4:1r,c5:1r,c6:1r,c7:1r,c8:1r,c9:1r,c10:1r,c11:1r,c12:1r',
                  'chromosomes_order=c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12',
                  '<<include etc/housekeeping.conf>>',
                  '<plots>',
                  

                  '<plot>',
                  paste0('r0=',rmin1[1],'r'),
                  paste0('r1=',rmax1[1],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat1.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line1.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[1],'r'),
                  paste0('r1=',rmax2[1],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                  '<plot>',
                  paste0('r0=',rmin1[2],'r'),
                  paste0('r1=',rmax1[2],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat2.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line2.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[2],'r'),
                  paste0('r1=',rmax2[2],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                                    '<plot>',
                  paste0('r0=',rmin1[3],'r'),
                  paste0('r1=',rmax1[3],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat3.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line3.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[3],'r'),
                  paste0('r1=',rmax2[3],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',



                  '<plot>',
                  paste0('r0=',rmin1[4],'r'),
                  paste0('r1=',rmax1[4],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat4.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line4.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[4],'r'),
                  paste0('r1=',rmax2[4],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',


                  '<plot>',
                  paste0('r0=',rmin1[5],'r'),
                  paste0('r1=',rmax1[5],'r'),
                  'type=scatter',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/scat5.txt',
                  'glyph=circle',
                 # 'glyph_size=40',
                  'stroke_color=black',
                  'stroke_thickness=1',
                  'orientation=out',
                  'min=0',
                  'max=7',
                  back1,
                  '</plot>',


# line plot                                    
                  '<plot>',
                  'type=tile',
                  'file=/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/line5.txt',
                  'thickness=6',
                  paste0('r0=',rmin2[5],'r'),
                  paste0('r1=',rmax2[5],'r'),
                  'orientation=center',
                  'layers=7',
                  'margin=0u',
                  back2,
                  '</plot>',




                  '</plots>')


write(circosStarter, file = '/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/circos.conf')


cmd <- paste("perl /Users/luisdiaz/Documents/chamba/BIOINFO_TOOLS/circos-0.67-7/bin/circos -conf ~/Documents/chamba/PHENOTYPING/QTLcolorPaper/circosFiles/circos.conf")
system(cmd)




saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataCNJ02.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1CNJ02.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3CNJ02.rds')


## a better exploration of color

tmp1<-scantwo(phenoData[[2]]$sep$cross,pheno.col=1,chr=3)
plot(clean(tmp1))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[2]]$oct$cross,mname1=find.marker(phenoData[[2]]$oct$cross,3,i),draw=F)
print(i)
}

mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)







####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################'
####################################################################################################################################################################################'
############################################################################     P I C T U R E S      ##############################################################################'
####################################################################################################################################################################################'
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES





digitalCN<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markCNJ_WISC.csv',header=TRUE,sep=',')



mypops<-gregexpr('CNJ02',digitalCN$id3)
tt<-sapply(mypops, "[[", 1)

f02<-which(tt==1)
f04<-which(tt==-1)


mlx<-lm(digitalCN$g.mean[f02]~digitalCN$id3[f02])
yy02<-unique(predict(mlx))
yyRes02<-tapply(digitalCN$g.mean[f02],digitalCN$id3[f02],sd)
yyRes02<-yyRes02[1:166]

names(yy02)<-names(yyRes02)
cor(yy02,yyRes02,use='complete.obs')

digData02<-data.frame(color=yy02,colovar=yyRes02)
f<-which(digData02[,1]>0.51)
digData02<-digData02[-f,]
gt1<-gsub(pattern = "-",replacement = "_",rownames(digData02))
gt2<-gsub(pattern = ".JPG",replacement = "",gt1)
rownames(digData02)<-gt2


mlx<-lm(digitalCN$g.mean[f04]~digitalCN$id3[f04])
yy04<-unique(predict(mlx))
yyRes04<-tapply(digitalCN$g.mean[f04],digitalCN$id3[f04],sd)
yyRes04<-yyRes04[167:235]

names(yy04)<-names(yyRes04)
cor(yy04,yyRes04,use='complete.obs')
digData04<-data.frame(color=yy04,colovar=yyRes04)
gt1<-gsub(pattern = "-",replacement = "_",rownames(digData04))
gt2<-gsub(pattern = ".JPG",replacement = "",gt1)
rownames(digData04)<-gt2



superMap<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/consensusMapAll2.csv',header=T,sep=',')
superMap<-superMap[,c('marker','LG','consensus')]
superMap<-superMap[order(superMap[,2],superMap[,3]),]
superMap$binID<-NA
superMap2<-numeric()
for (i in 1:12){
  f<-which(superMap$LG==i)
  mybins<-unique(superMap$consensus[f])
  mybinsID<-paste0('bin_',i,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap$consensus[f]==mybins[j])
    superMap$binID[f[f2]]<-mybinsID[j]      
    superMap2<-rbind(superMap2,superMap[f[f2[1]],])
  }
}

rownames(superMap2)<-superMap2$marker



# cnj02
geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/CNJ02_AllASMapData.csv',header=T,sep=',')  



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



f<-intersect(rownames(digData02),rownames(gData))
y2<-digData02[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]


write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))



phenoData<-list()
phenoData[[1]]<-list(cross=calc.genoprob(superF,step=0,map.function="kosambi"),blups=digData02)

phenoData[[1]]$blups<-phenoData[[1]]$blups*-1
#### cnj04

geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/CNJ04_AllASMapData.csv',header=T,sep=',')  



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



f<-intersect(rownames(digData04),rownames(gData))
y2<-digData04[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]


write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))



phenoData[[2]]<-list(cross=calc.genoprob(superF,step=0,map.function="kosambi"),blups=digData04)


phenoData[[2]]$blups<-phenoData[[2]]$blups*-1





sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties)
  return(x1)
}



n.cores <- 6

cl <- snow::makeCluster(n.cores)
registerDoSNOW(cl)
for (i in 1:2){
superOut <- foreach(ph = 1:2,.packages = 'qtl') %dopar% perroQTL(phenoData[[i]]$cross,ph,max.qtl=5,penalties=sixo)
phenoData[[i]]$scantwo<-superOut
print(i)
}

snow::stopCluster(cl)







for (tri in 1:2){
  modelRes<-vector("list", 4) 
  for (tri2 in 1:2){
    rob1<-phenoData[[tri]]$scantwo[[tri2]]
    if (length(rob1)==0){
      
    }else{
      superQTL<-makeqtl(phenoData[[tri]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
      md1<-fitqtl(phenoData[[tri]]$cross,pheno.col = tri,superQTL,formula = formula(rob1),get.ests=F)
      x<-summary(md1)
      modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                      drop=as.data.frame(x$result.drop),
                      model=rob1,
                      inter=numeric())
      for (j in 1:length(rob1$chr)){
        mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
        modelRes[[tri2]]$inter[j]<-mylod
      }
    }
  }
  phenoData[[tri]]$qtldata<-modelRes
}


phenoData[[1]]$trait<-'cnj02'
phenoData[[2]]$trait<-'cnj04'

# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

myt1<-colnames(phenoData[[1]]$blups)
jaguar1<-data.frame(population=numeric(),trait=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:2){
  t1<-phenoData[[i]]$trait
  for (j in 1:2){
    nQTL<-length(phenoData[[i]]$qtldata[[j]]$model$chr)
    if (nQTL==1){
     jaguar1[cont,]<-c(t1,myt1[j],as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr),phenoData[[i]]$qtldata[[j]]$model$pos,paste0(x$binID[idx],'_',x$marker[idx]),rep(phenoData[[i]]$qtldata[[j]]$variance,2),phenoData[[i]]$qtldata[[j]]$inter)
      cont<-cont+1
      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]]$qtldata[[j]]$drop))
      x2<-phenoData[[i]]$qtldata[[j]]$drop[f2,]
      jaguar1[cont,]<-c(t1,myt1[j],mychr,mypos,mymarker,x2$'%var',phenoData[[i]]$qtldata[[j]]$variance,phenoData[[i]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)


saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_bothPop.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1Pictures_bothPop.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3_bothPop.rds')




## a better exploration of color

tmp1<-scantwo(phenoData[[2]]$cross,pheno.col=1,chr=3)
plot(clean(tmp1))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[2]]$cross,mname1=find.marker(phenoData[[2]]$cross,3,i),draw=F)
print(i)
}

mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(0.35,.46),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)





##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
##### GRYG #######################################################################################################################################
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES
# N O T E: I moved the files from QTL_project to /Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/rawdata_PICTURES


GRYG1<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2014.csv',header=TRUE,sep=',')
GRYG1$year<-'a'
GRYG2<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2015.csv',header=TRUE,sep=',')
GRYG2$year<-'b'
colnames(GRYG2)<-colnames(GRYG1)
GRYG3<-read.csv('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTL_project/markGRYG_2016.csv',header=TRUE,sep=',')
GRYG3$year<-'c'
colnames(GRYG3)<-colnames(GRYG1)



GRYG1$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG1$id3)))
GRYG2$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG2$id3)))
GRYG3$id3<-as.factor(as.numeric(gsub('.JPG','',GRYG3$id3)))

phenoData<-list()

phenoData[[1]]<-list(rawData=GRYG1)
phenoData[[2]]<-list(rawData=GRYG2)
phenoData[[3]]<-list(rawData=GRYG3)





for (i in 1:3){
  mdl1<-lmer(b.mean.1~(1|id3),data=phenoData[[i]]$rawData)
  blup1<-ranef(mdl1)$id3
  myColVar<-tapply(phenoData[[i]]$rawData$b.mean.1,list(myID=phenoData[[i]]$rawData$id3),sd)
  f<-intersect(names(myColVar),rownames(myColVar))
  blupsX<-as.data.frame(cbind(blup1[f,],myColVar[f]))
  colnames(blupsX)<-c('color','colorvar')
  phenoData[[i]]$blups<-blupsX
}

pairs(phenoData[[3]]$blups)

f<-intersect(rownames(phenoData[[1]]$blups),rownames(phenoData[[2]]$blups))
f<-intersect(f,rownames(phenoData[[3]]$blups))



for (i in 1:3){
  phenoData[[i]]$blups2<-phenoData[[i]]$blups[f,]
  rownames(phenoData[[i]]$blups2)<-paste0('P',rownames(phenoData[[i]]$blups2))
}


# in the following code I will combine years into a single list element. 

phenoDataX<-phenoData


# doing multivariate




geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/Gryg_AllASMapData.csv',header=T,sep=',')  
rownames(geno)<-geno$X
geno<-geno[,-c(1:6)]
geno2<-atcg1234(t(geno))


phenoData<-list()
mytraits<-c('color','colorvar')
myyears<-c('2014','2015','2016')
for (i in 1:2){

phenoA<-cbind(phenoDataX[[1]]$blups2[,i],phenoDataX[[2]]$blups2[,i],phenoDataX[[3]]$blups2[,i])

rownames(phenoA)<-rownames(phenoDataX[[1]]$blups2)


f<-intersect(rownames(phenoA),rownames(geno2))
phenoA2<-phenoA[f,]
geno3<-geno2[f,]
A<-A.mat(geno3)


Za <- diag(nrow(A)) 
ETA.A <- list(add=list(Z=Za,K=A)) 
repo1<-mmer(Y=phenoA2,Z=ETA.A,silent=FALSE,MVM=TRUE,EIGEND = FALSE,draw=TRUE)
blups1<-repo1$u.hat$add
var1<-repo1$var.comp


phenoData[[i]]<-list(blups=as.matrix(blups1),vcov=var1,trait=mytraits[i])

}


for (i in 1:2){
  phenoData[[i]]$blups<-phenoData[[i]]$blups*-1
}


# calculating h2 and correlations

for (i in 1:2){
  
    gv<-diag(phenoData[[i]]$vcov$add)
    rv<-diag(phenoData[[i]]$vcov$Residual)
    phenoData[[i]]$h2<-gv/(gv+rv)
    x<-cor(phenoData[[i]]$blups,use='complete.obs')
    phenoData[[i]]$gcor<-x
    colnames(phenoData[[i]]$gcor)<-myyears
    rownames(phenoData[[i]]$gcor)<-myyears
    names(phenoData[[i]]$h2)<-myyears
    rownames(phenoData[[i]]$blups)<-rownames(phenoA2)
    colnames(phenoData[[i]]$blups)<-myyears
}




superMap<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/consensusMapAll2.csv',header=T,sep=',')
superMap<-superMap[,c('marker','LG','consensus')]
superMap<-superMap[order(superMap[,2],superMap[,3]),]
superMap$binID<-NA
superMap2<-numeric()
for (i in 1:12){
  f<-which(superMap$LG==i)
  mybins<-unique(superMap$consensus[f])
  mybinsID<-paste0('bin_',i,'@',mybins,'cM')
  
  for (j in 1:length(mybins)){
    f2<-which(superMap$consensus[f]==mybins[j])
    superMap$binID[f[f2]]<-mybinsID[j]      
    superMap2<-rbind(superMap2,superMap[f[f2[1]],])
  }
}

rownames(superMap2)<-superMap2$marker

## genetic analysis
geno<-read.table('/Users/luisdiaz/Documents/chamba/PHENOTYPING/QTLcolorPaper/genetic_data/Gryg_AllASMapData.csv',header=T,sep=',')  

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

for (ii in 1:2){


f<-intersect(rownames(phenoData[[ii]]$blups),rownames(gData))
y2<-phenoData[[ii]]$blups[f,]

gData2<-gData[f,]
gData2[1:10,1:10]


f<-intersect(colnames(gData2),superMap2$marker)

gData3<-gData2[,f]
superMap3<-superMap2[f,]

gData4<-rbind(superMap3$LG,superMap3$consensus,gData3)
gData4[1:10,1:10]

colnames(gData4)<-superMap3$binID


gData5<-data.frame(rbind('','',y2),gData4)
rownames(gData5)[1:2]<-c('chr','pos')
gData5[1:10,1:20]

#write.csv(file='/crdata4/luis4/QTLcolorServer_copy/testQTL',gData5,row.names = FALSE)
write.csv(file='temp/testQTL',gData5,row.names = FALSE)


superF <- (read.cross(format = "csv", file='temp/testQTL',genotypes = NULL))
#superF <- (read.cross(format = "csv", file='/crdata4/luis4/QTLcolorServer_copy/testQTL',genotypes = NULL))


phenoData[[ii]]$cross<-numeric()
phenoData[[ii]]$cross <- calc.genoprob(superF,step=0,map.function="kosambi") 

}




sixo<-c(4.081565, 9.150947, 6.473634)

perroQTL<-function(cross,phenotype,max.qtl,penalties){
  x1<-stepwiseqtl(cross,pheno.col=phenotype,max.qtl=max.qtl,method='hk',penalties=penalties)
  return(x1)
}



n.cores <- 6

cl <- snow::makeCluster(n.cores)
registerDoSNOW(cl)
for (i in 1:2){
superOut <- foreach(ph = 1:3,.packages = 'qtl') %dopar% perroQTL(phenoData[[i]]$cross,ph,max.qtl=10,penalties=sixo)
phenoData[[i]]$scantwo<-superOut
print(i)
}

snow::stopCluster(cl)





for (tri in 1:2){
  modelRes<-vector("list", 4) 
  for (tri2 in 1:3){
    rob1<-phenoData[[tri]]$scantwo[[tri2]]
    if (length(rob1)==0){
      
    }else{
      superQTL<-makeqtl(phenoData[[tri]]$cross,chr=rob1$chr,pos=rob1$pos,what = 'prob')
      md1<-fitqtl(phenoData[[tri]]$cross,pheno.col = tri2,superQTL,formula = formula(rob1),get.ests=F)
      x<-summary(md1)
      modelRes[[tri2]]<-list(variance=x$result.full[1,5],
                      drop=as.data.frame(x$result.drop),
                      model=rob1,
                      inter=numeric())
      for (j in 1:length(rob1$chr)){
        mylod<-diff(lodint(rob1,qtl.index = j)[c(1,3),2])
        modelRes[[tri2]]$inter[j]<-mylod
      }
    }
  }
  phenoData[[tri]]$qtldata<-modelRes
}



# making a huge list of QTLs and fushion those that are the same based on a 5cM interval.

jaguar1<-data.frame(trait=numeric(),date=numeric(),chr=numeric(),position=numeric(),nearest.marker=numeric(),marker.variance=numeric(),model.variance=numeric(),interval=numeric())
cont<-1
for (i in 1:2){
  t1<-phenoData[[i]]$trait
  for (j in 1:3){
    nQTL<-length(phenoData[[i]]$qtldata[[j]]$model$chr)
    if (nQTL==1){
     jaguar1[cont,]<-c(t1,myt1[j],as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr),phenoData[[i]]$qtldata[[j]]$model$pos,paste0(x$binID[idx],'_',x$marker[idx]),rep(phenoData[[i]]$qtldata[[j]]$variance,2),phenoData[[i]]$qtldata[[j]]$inter)
      cont<-cont+1
      }else{
    for (k in 1:nQTL){
      mypos<-phenoData[[i]]$qtldata[[j]]$model$pos[k]
      mychr<-as.numeric(phenoData[[i]]$qtldata[[j]]$model$chr[k])

      f<-which(superMap3$LG==mychr)
      x<-superMap3[f,]
      idx<-which.min(abs(x$consensus-mypos))
      mymarker<-paste0(x$binID[idx],'_',x$marker[idx])
      f2<-intersect(phenoData[[i]]$qtldata[[j]]$model$name[k],rownames(phenoData[[i]]$qtldata[[j]]$drop))
      x2<-phenoData[[i]]$qtldata[[j]]$drop[f2,]
      jaguar1[cont,]<-c(t1,myyears[j],mychr,mypos,mymarker,x2$'%var',phenoData[[i]]$qtldata[[j]]$variance,phenoData[[i]]$qtldata[[j]]$inter[k])
      cont<-cont+1
    }
  }
  }
}


jaguar1$chr<-as.numeric(jaguar1$chr)
jaguar1$position<-as.numeric(jaguar1$position)
jaguar1$marker.variance<-as.numeric(jaguar1$marker.variance)
jaguar1$model.variance<-as.numeric(jaguar1$model.variance)
jaguar1$interval<-as.numeric(jaguar1$interval)


saveRDS(phenoData,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_GRYG_3years.rds')
saveRDS(jaguar1,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/jaguar1Pictures_GRYG_3years.rds')
saveRDS(superMap3,'~/Documents/chamba/PHENOTYPING/QTLcolorPaper/superMap3_GRYG_3years.rds')



## a better exploration of color

phenoData<-readRDS('~/Documents/chamba/PHENOTYPING/QTLcolorPaper/phenoDataPictures_GRYG_3years.rds')

tmp1<-scantwo(phenoData[[1]]$cross,pheno.col=1,chr=3)

postscript(file="~/Documents/chamba/PHENOTYPING/QTLcolorPaper/scantwogryg.eps", 
            width=10, 
            height=10, 
            horizontal=TRUE) 
plot(clean(tmp1))
dev.off()


plot(scanone(phenoData[[1]]$cross,pheno.col=1,chr=3))


tmp2<-list()
for (i in 1:90){
tmp2[[i]]<-effectplot(phenoData[[1]]$cross,mname1=find.marker(phenoData[[1]]$cross,3,i),draw=F)
print(i)
}
postscript(file="~/Documents/chamba/PHENOTYPING/QTLcolorPaper/effectplotGRYG.eps", 
            width=20, 
            height=10, 
            horizontal=TRUE) 
mycolsX<-sample(brocolors('crayons'),4)
plot(rep(1,4),tmp2[[1]]$Means,xlim=c(0,100),ylim=c(-.03,.03),col=mycolsX,pch=20,cex=2,ylab='allele effect',xlab='marker position (cM), chr 3')
for (i in 1:4){
lines(c(1,1),c(tmp2[[1]]$Means[i]-tmp2[[1]]$SEs[i],tmp2[[1]]$Means[i]+tmp2[[1]]$SEs[i]),col=mycolsX[i])
}
for (i in 2:90){
points(rep(i,4),tmp2[[i]]$Means,xlim=c(0,100),ylim=c(-5,5),col=mycolsX,pch=20,cex=2)
for (j in 1:4){
lines(c(i,i),c(tmp2[[i]]$Means[j]-tmp2[[i]]$SEs[j],tmp2[[i]]$Means[j]+tmp2[[i]]$SEs[j]),col=mycolsX[j])
}

}
abline(h=0,col='black')
legend('top',legend=c('AC','BC','AD','BD'),col=mycolsX,pch=20,horiz=T)

dev.off()




rm(list=ls())
dir()
library(lmerTest)
library(lattice)
library(calibrate)

# loading the data
#datar <- read.csv("CQxMQ_PhenotypicData_EduardoUsed.csv",header=T)
#datar <- read.csv("CQxMQ_PhenotypicData_EduardoUsed.csv",header=T)
datar <- read.csv("CNJ02_PlotData.csv",header=T)

head(datar)
tail(datar)

################################################################################
# use NNA function to get the nearest neighbour effect--------------------------
################################################################################

Tacy<-datar[,c(1:5,10)]
Tacy[383,6]<-NA
boxplot(Tacy$Tacy~Tacy$Year+Tacy$Month)
plot(Tacy$Tacy[Tacy$Month=="September"],Tacy$Tacy[Tacy$Month=="October"])
Tacy <- NNA(datar=Tacy, plantid="Plant", repetition="Month", covar="Year", response="Tacy")


Brix<-datar[,c(1:5,11)]
Brix[563,6]<-NA
Brix[259,6]<-8.8
boxplot(Brix$Brix~Brix$Covariate_Year+Brix$Month)
plot(Brix$Brix[(Brix$Year=="y2011")],Brix$Brix[Brix$Year=="y2012"])
Brix<- NNA(datar=Brix, plantid="Plant", repetition="Month", covar="Year", response="Brix")

Tacid<-datar[,c(1:5,12)]
Tacid[704,6]<-NA
Tacid[751,6]<-NA
boxplot(Tacid$TA~Tacid$Year+Tacid$Month)
Tacid <- NNA(datar=Tacid, plantid="Plant", repetition="Month", covar="Year", response="TA")

PACy<-datar[,c(1:5,13)]
PACy[209,6]<-1.2
boxplot(PACy$PA~PACy$Year+PACy$Month)
PACy <- NNA(datar=PACy, plantid="Plant", repetition="Month", covar="Year", response="PAC")
boxplot(PACy$PA~PACy$Year+PACy$Month)

################################################################################
# Run Models for TACY, BRIX, TA, PAC
################################################################################

#Create three matrices to store the variance components and significance levels
#for the random effects models.

#response~(1|Geno)+(1|Year)+(1|Geno:Year)+(1|NS)+(1|EW)+(1|NS:EW)
#uses month as a rep
yearMatrix<-matrix(" ",6,4)
rownames(yearMatrix)<-c("Genotype","Year","Genotype x Year","North/South", "East/West", "Residual")
colnames(yearMatrix)<-c("Tacy","Brix","TA","PAC")

#response~(1|Geno)+(1|Month)+(1|Geno:Month)+(1|NS)+(1|EW)+(1|NS:EW)
#uses the year as a rep
monthMatrix<-matrix(" ",6,4)
rownames(monthMatrix)<-c("Genotype","Month","Genotype x Month","North/South", "East/West", "Residual")
colnames(monthMatrix)<-c("Tacy","Brix","TA","PAC")

#response(Oct or Sept only)~(1|Geno)+(1|NS)+(1|EW)+(1|NS:EW)
#uses the year as a rep
septOctMatrix<-matrix(" ",4,8)
rownames(septOctMatrix)<-c("Genotype","North/South", "East/West", "Residual")
colnames(septOctMatrix)<-c("Tacy","Brix","TA","PAC","Tacy","Brix","TA","PAC")



# TACY -------------------------------------------------------------------------

#Get Tacy Heatmaps of model with only Genotype Effect
TacyHeat<-Heatplots1(x = Tacy)
print(c(TacyHeat[[7]],TacyHeat[[8]],TacyHeat[[9]],layout=c(1,3)))
#Get Tacy Heatmaps of model with only Genotype Effect
TacyHeat<-Heatplots2(x = Tacy)
print(c(TacyHeat[[1]],TacyHeat[[2]],TacyHeat[[3]],layout=c(1,3)))

#Get Matrix Spearman Correlations between Tacy months and years
TacyMatSP<-CorMatrixMonthYear(x = Tacy,y = 6,z = 4,l = 5,n = 3)
TacyMatbyYear<-CorMatrixMonthYear(x = Tacy,y = 6,z = 4,l = "nope",n = 3)
TacyMatbyMonth<-CorMatrixMonthYear(x = Tacy,y = 6,z = 5,l = "nope",n = 1)
TacyMatSP$Spearman
TacyMatbyYear$Spearman
head(TacyMatSP$theMatrix)
image(TacyMatSP$theMatrix[order(TacyMatSP$theMatrix[,1]),])
write.csv(TacyMatSP$Spearman,file="TacySpearman.csv")
write.csv(TacyMatbyYear$Spearman,file="TacySpearmanByYear.csv")
write.csv(TacyMatbyMonth$Spearman,file="TacySpearmanByMonth.csv")
#Tacy Models--------------------------------------------------------------------
#Month as rep
TacyMod2<-lmer(Tacy~(1|Genotype)+(1|Covariate_Year)+(1|Genotype:Covariate_Year)+(1|NS)+(1|EW),data=Tacy[complete.cases(Tacy),])
T1<-summary(TacyMod2)
T2<-rand(TacyMod2)
plot(fitted(TacyMod2),resid(TacyMod2),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(TacyMod2), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(TacyMod2))
qqline(resid(TacyMod2))
yearMatrix[,1]<-byYear(T1,T2)

#Year as rep
TacyMod3<-lmer(Tacy~(1|Genotype)+(1|Month)+(1|Genotype:Month)+(1|NS)+(1|EW),data=Tacy[complete.cases(Tacy),])
T3<-summary(TacyMod3)
T4<-rand(TacyMod3)
plot(fitted(TacyMod3),resid(TacyMod3),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(TacyMod3), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(TacyMod3))
qqline(resid(TacyMod3))
monthMatrix[,1]<-byBothMonths(T3,T4)

#September model
TacyModSept<-lmer(Tacy~(1|Genotype)+(1|NS)+(1|EW),data=Tacy[complete.cases(Tacy),][Tacy$Month=="September",])
T5<-summary(TacyModSept)
T6<-rand(TacyModSept)
plot(fitted(TacyModSept),resid(TacyModSept),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(TacyModSept), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(TacyModSept))
qqline(resid(TacyModSept))
septOctMatrix[,1]<-oneMonth(T5,T6)

#October model
TacyModOct<-lmer(Tacy~(1|Genotype)+(1|NS)+(1|EW),data=Tacy[complete.cases(Tacy),][Tacy$Month=="October",])
T7<-summary(TacyModOct)
T8<-rand(TacyModOct)
septOctMatrix[,5]<-oneMonth(T7,T8)


TacyBLUPs<-data.frame(coef(TacyMod2)$Genotype,coef(TacyMod3)$Genotype,coef(TacyModSept)$Genotype,coef(TacyModOct)$Genotype)
colnames(TacyBLUPs)<-c("TacyYear","TacyMonth","TacySept","TacyOct")
#write.csv(TacyBLUPs,file="TacyBLUPs.csv")


################################################################################
#Brix -------------------------------------------------------------------------
################################################################################

#Get Brix HeatMaps
BrixHeat<-Heatplots1(x = Brix)
print(c(BrixHeat[[7]],BrixHeat[[7]],BrixHeat[[9]],layout=c(1,3)))

#Get Matrix Spearman Correlations between Brix months and years

BrixMatSP<-CorMatrixMonthYear(x = Brix,y = 6,z = 4,l = 5,n = 3)
BrixMatbyYear<-CorMatrixMonthYear(x = Brix,y = 6,z = 4,l = "nope",n = 3)
BrixMatbyMonth<-CorMatrixMonthYear(x = Brix,y = 6,z = 5,l = "nope",n = 1)
BrixMatSP$Spearman
BrixMatbyYear$Spearman
head(BrixMatSP$theMatrix)
write.csv(BrixMatSP$Spearman,file="BrixSpearman.csv")
write.csv(BrixMatbyYear$Spearman,file="BrixSpearmanByYear.csv")
write.csv(BrixMatbyMonth$Spearman,file="BrixSpearmanByMonth.csv")
#Brix Models

#Brix using Month as rep: No GenotypexYear interaction! Use this model for BLUPs!
BrixMod4<-lmer(Brix~(1|Genotype)+(1|Covariate_Year)+(1|Genotype:Covariate_Year)+(1|NS)+(1|EW),data=Brix[complete.cases(Brix),])
B1<-summary(BrixMod4)
B2<-rand(BrixMod4)
plot(fitted(BrixMod4),resid(BrixMod4),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(BrixMod4), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(BrixMod4))
qqline(resid(BrixMod4))
yearMatrix[,2]<-byYear(B1,B2)

#Brix using year as rep: No GenotypeXMonth interaction so can use month as rep!
BrixMod3<-lmer(Brix~(1|Genotype)+(1|Month)+(1|Genotype:Month)+(1|NS)+(1|EW),data=Brix[complete.cases(Brix),])
B3<-summary(BrixMod3)
B4<-rand(BrixMod3)
plot(BrixMod3)
qqnorm(resid(BrixMod3))
qqline(resid(BrixMod3))
monthMatrix[,2]<-byBothMonths(B3,B4)


#September Brix, no significant Genotype variance
BrixModSept<-lmer(Brix~(1|Genotype)+(1|NS)+(1|EW),data=Brix[complete.cases(Brix),][Brix$Month=="September",])
B5<-summary(BrixModSept)
B6<-rand(BrixModSept)
plot(BrixModSept)
qqnorm(resid(BrixModSept))
qqline(resid(BrixModSept))
septOctMatrix[,2]<-oneMonth(B5,B6)

#October Brix, significant Genotype variance! use for QTL mapping also.
BrixModOct<-lmer(Brix~(1|Genotype)+(1|NS)+(1|EW),data=Brix[complete.cases(Brix),][Brix$Month=="October",])
B7<-summary(BrixModOct)
B8<-rand(BrixModOct)
plot(BrixModOct)
qqnorm(resid(BrixModOct))
qqline(resid(BrixModOct))
septOctMatrix[,6]<-oneMonth(B7,B8)

BrixBLUPs<-data.frame(coef(BrixMod4)$Genotype,coef(BrixMod3)$Genotype,coef(BrixModSept)$Genotype,coef(BrixModOct)$Genotype)
colnames(BrixBLUPs)<-c("BrixYear","BrixMonth","BrixSept","BrixOct")
write.csv(BrixBLUPs,file="BrixBLUPs.csv")

################################################################################
# TA----------------------------------------------------------------------------
################################################################################

#Get Brix HeatMaps
Tacider<-Heatplots1(x = Tacid)
print(c(Tacider[[7]],Tacider[[7]],Tacider[[9]],layout=c(1,3)))

#Get Matrix Spearman Correlations between Brix months and years
TacidMatSP<-CorMatrixMonthYear(x = Tacid,y = 6,z = 4,l = 5,n = 3)
TacidMatbyYear<-CorMatrixMonthYear(x = Tacid,y = 6,z = 4,l = "nope",n = 3)
TacidMatbyMonth<-CorMatrixMonthYear(x = Tacid,y = 6,z = 5,l = "nope",n = 1)
TacidMatSP$Spearman
TacidMatbyYear$Spearman
head(TacidMatSP$theMatrix)
write.csv(TacidMatSP$Spearman,file="TacidSpearman.csv")
write.csv(TacidMatbyYear$Spearman,file="TacidSpearmanByYear.csv")
write.csv(TacidMatbyMonth$Spearman,file="TacidSpearmanByMonth.csv")

#TA Models------------------------------------------------------------------------

#TA using Month as rep: No GenotypexYear interaction! Use this model for BLUPs and QTLs!
TAMod4<-lmer(TA~(1|Genotype)+(1|Covariate_Year)+(1|Genotype:Covariate_Year)+(1|NS)+(1|EW),data=Tacid[complete.cases(Tacid),])
TA1<-summary(TAMod4)
TA2<-rand(TAMod4)
plot(fitted(TAMod4),resid(TAMod4),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(TAMod4), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(TAMod4))
qqline(resid(TAMod4))
yearMatrix[,3]<-byYear(TA1,TA2)

#TA using year as rep: No GenotypeXMonth interaction and barely significant Month so can use month as rep!
TAMod3<-lmer(TA~(1|Genotype)+(1|Month)+(1|Genotype:Month)+(1|NS)+(1|EW),data=Tacid[complete.cases(Tacid),])
TA3<-summary(TAMod3)
TA4<-rand(TAMod3)
plot(TAMod3)
qqnorm(resid(TAMod3))
qqline(resid(TAMod3))
monthMatrix[,3]<-byBothMonths(TA3,TA4)

#September TA, significant Genotype variance use for QTL mapping
TAModSept<-lmer(TA~(1|Genotype)+(1|NS)+(1|EW),data=Tacid[complete.cases(Tacid),][Tacid$Month=="September",])
TA5<-summary(TAModSept)
TA6<-rand(TAModSept)
plot(TAModSept)
qqnorm(resid(TAModSept))
qqline(resid(TAModSept))
septOctMatrix[,3]<-oneMonth(TA5,TA6)

#October TA, significant Genotype variance! use for QTL mapping also.
TAModOct<-lmer(TA~(1|Genotype)+(1|NS)+(1|EW),data=Tacid[complete.cases(Tacid),][Tacid$Month=="October",])
TA7<-summary(TAModOct)
TA8<-rand(TAModOct)
plot(TAModOct)
qqnorm(resid(TAModOct))
qqline(resid(TAModOct))
septOctMatrix[,7]<-oneMonth(TA7,TA8)

TABLUPs<-data.frame(coef(TAMod4)$Genotype,coef(TAMod3)$Genotype,coef(TAModSept)$Genotype,coef(TAModOct)$Genotype)
colnames(TABLUPs)<-c("TAYear","TAMonth","TASept","TAOct")
write.csv(TABLUPs,file="TABLUPs.csv")

################################################################################
# PAC --------------------------------------------------------------------------
################################################################################
#Get PAC HeatMaps
PACyer<-Heatplots1(x = PACy)
print(c(PACyer[[7]],PACyer[[7]],PACyer[[9]],layout=c(1,3)))

#Get Matrix Spearman Correlations between PAC months and years
PACMatSP<-CorMatrixMonthYear(x = PACy,y = 6,z = 4,l = 5,n = 3)
PACMatbyYear<-CorMatrixMonthYear(x = PACy,y = 6,z = 4,l = "nope",n = 3)
PACMatbyMonth<-CorMatrixMonthYear(x = PACy,y = 6,z = 5,l = "nope",n = 1)
PACMatSP$Spearman
PACMatbyYear$Spearman
head(PACMatSP$theMatrix)
write.csv(PACMatSP$Spearman,file="PACSpearman.csv")
write.csv(PACMatbyYear$Spearman,file="PACSpearmanByYear.csv")
write.csv(PACMatbyMonth$Spearman,file="PACSpearmanByMonth.csv")
#PAC Models--------------------------------------------------------------------

#PAC using Month as rep: No GenotypexYear interaction! Use this model for BLUPs and QTLs!
PACyMod4<-lmer(PAC~(1|Genotype)+(1|Covariate_Year)+(1|Genotype:Covariate_Year)+(1|NS)+(1|EW),data=PACy[complete.cases(PACy),])
P1<-summary(PACyMod4)
P2<-rand(PACyMod4)
plot(fitted(PACyMod4),resid(PACyMod4),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(PACyMod4), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(PACyMod4))
qqline(resid(PACyMod4))
yearMatrix[,4]<-byYear(P1,P2)

#PAC using year as rep: No GenotypeXMonth interaction so can use month as rep!
PACyMod3<-lmer(PAC~(1|Genotype)+(1|Month)+(1|Genotype:Month)+(1|NS)+(1|EW),data=PACy[complete.cases(PACy),])
P3<-summary(PACyMod3)
P4<-rand(PACyMod3)
plot(fitted(PACyMod3),resid(PACyMod3),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(PACyMod3), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(PACyMod3))
qqline(resid(PACyMod3))
monthMatrix[,4]<-byBothMonths(P3,P4)

#September PAC, significant Genotype variance use for QTL mapping
PACyModSept<-lmer(PAC~(1|Genotype)+(1|NS)+(1|EW),data=PACy[complete.cases(PACy),][PACy$Month=="September",])
P5<-summary(PACyModSept)
P6<-rand(PACyModSept)
plot(fitted(PACyModSept),resid(PACyModSept),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(PACyModSept), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(PACyModSept))
qqline(resid(PACyModSept))
septOctMatrix[,4]<-oneMonth(P5,P6)

#October PAC, significant Genotype variance! use for QTL mapping also.
PACyModOct<-lmer(PAC~(1|Genotype)+(1|NS)+(1|EW),data=PACy[complete.cases(PACy),][PACy$Month=="October",])
P7<-summary(PACyModOct)
P8<-rand(PACyModOct)
plot(fitted(PACyModOct),resid(PACyModOct),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
hist(resid(PACyModOct), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")
qqnorm(resid(PACyModOct))
qqline(resid(PACyModOct))
septOctMatrix[,8]<-oneMonth(P7,P8)



PACyBLUPs<-data.frame(coef(PACyMod4)$Genotype,coef(PACyMod3)$Genotype,coef(PACyModSept)$Genotype,coef(PACyModOct)$Genotype)
colnames(PACyBLUPs)<-c("PACyYear","PACyMonth","PACySept","PACyOct")
write.csv(PACyBLUPs,file="PACyBLUPs.csv")



write.csv(yearMatrix,file="VarianceComponentsYearModel.csv")
write.csv(monthMatrix,file="VarianceComponentsMonthModel.csv")
write.csv(septOctMatrix,file="VarianceComponentsSeptOctModel.csv")







# models to compare
mo0 <- lmer(TY~(1|Genotype), data=ready)
mo2 <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year), data=ready)
mo3 <- lmer(TY~(1|Genotype) + (1|Covariate_Year) + (1|Genotype:Covariate_Year) + (1|NS) + (1|Covariate_Year:NS), data=ready)
#mo3 <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year) + (1|EW)+ (1|NS), data=ready)
mo4 <- lmer(TY~(1|Plant) + (1|EW), data=ready)
rand(mo3)
rand(mo4)
summary(mo3)
# add residuals after the best model fitting
ready$resmo0  <- NA
ready$resmo3  <- NA
f <- as.numeric(names(resid(mo3)))
g <- as.numeric(names(resid(mo0)))
ready[f,"resmo3"] <- resid(mo3)
ready[g,"resmo0"] <- resid(mo0)
my_palette <- colorRampPalette(c("red", "white","red"))(n = 599)
col_breaks = c(seq(-500,-100,length=200),seq(-99,200,length=200),seq(201,600,length=200))
levelplot(resmo3 ~Col*Row,col.regions = my_palette, at=col_breaks, main="residuals distribution for adjusted model", data=ready[(ready[,4]=="y2011")&(ready[,5]=="September"),])
levelplot(resmo0 ~Col*Row,col.regions = my_palette, at=col_breaks, main="residuals distribution for unadjusted model", data=ready[(ready[,4]=="y2011")&(ready[,5]=="September"),])

attach(ready)
## Compare BLUP to line averages on a scatterplot
lmean = tapply(TY, Plant, na.rm=T, mean) # applying the function mean to the variable TY(yield), grouping by LINE
plot(x=ranef(mo3)$Plant[,1], y=lmean, col="blue", main="Linear relationship between BLUPs and regular means",xlab="BLUP",ylab="Averages") # in x we have the BLUP's, in y the regular averages of the Plants
#dotplot(ranef(mo3)[[3]])
# compare model 1 that only uses the Plant effect against 
plot(x=ranef(mo3)$Plant[,1], y=ranef(mo2)$Plant[,1], col="blue", main="Linear relationship between BLUPs and regular means",xlab="BLUP",ylab="Averages") # in x we have the BLUP's, in y the regular averages of the Plants
# compare regular model against the NNA adjustment
plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for different mixed models")
points(x=fitted(mo2),y=resid(mo2), col="blue",pch=4)
points(x=fitted(mo3),y=resid(mo3), col="red",pch=4)
legend("topleft", legend=c("G","G:Y","G:Y + NNA"), lty=c(1,1,1), col=c("black","blue","red"))
# in case is not easy to see
layout(matrix(1:3,1,3))
#plot(1,1)
plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model",ylim=c(-300,700), main="Residual vs Fitted values for dummy model")
plot(x=fitted(mo2),y=resid(mo2),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", ylim=c(-300,700), main="Residual vs Fitted values for G*Y model")
plot(x=fitted(mo3),y=resid(mo3),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", ylim=c(-300,700), main="Residual vs Fitted values for G*Y + NNA")

plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Year model using months as rep")
density(resid(mo0), xlab="residuals", ylab="frequency", main="Residual densitys for Year model using months as rep")

plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for Month model using years as rep")
density(resid(mo0), xlab="residuals", ylab="frequency", main="Residual density Month model using years as rep")

plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for September Model")
density(resid(mo0), xlab="residuals", ylab="frequency", main="Residual density for September")


plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", main="Residual vs Fitted values for October Model")
density(resid(mo0), xlab="residuals", ylab="frequency", main="Residual density for October")
































#### comparing old summaries and new summaries
summary(mold <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year), data=ready))
summary(mox <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year) + (1|EW)  + 
                      (1|Year:NS), data=ready))
rand(mox)
(ranifs <- coef(mox)[[3]])
write.csv(ranifs, "nna3.csv")
################################
################################
################################
data.tre <- ready
#################
data.tre$r <- data.tre$Row
data.tre$r2 <- data.tre$Row^2
data.tre$r3 <- data.tre$Row^3
data.tre$c <- data.tre$Col
data.tre$c2 <- data.tre$Col^2
data.tre$c3 <- data.tre$Col^3
# trend analysis
data.tre$r <- lm(data.tre$r ~ 1)$res
data.tre$r2 <- lm(data.tre$r2 ~ 1 + data.tre$r)$res
data.tre$r3 <- lm(data.tre$r3 ~ 1 + data.tre$r + data.tre$r2)$res
data.tre$c <- lm(data.tre$c ~ 1)$res
data.tre$c2 <- lm(data.tre$c2 ~ 1 + data.tre$c)$res
data.tre$c3 <- lm(data.tre$c3 ~ 1 + data.tre$c + data.tre$c2)$res

layout(matrix(1,1,1))
with(data.tre, interaction.plot(Col,Row,res, col=2:5))
plot(res~Row, data=data.tre)
plot(res~Col, data=data.tre)
anova(lm(res ~ r*r2*r3*c*c2*c3, data=data.tre))
# linear effects for row and column exist

trend <- lmer(TY ~ (1|Plant) + (1|Year) + (1|Plant:Year) + (1|Row) + (1|Col),control=lmerControl(optimizer="bobyqa"), data=data.tre)
summary(trend)
rand(trend)

################################
################################
################################
################################
# checking correlations

vars <- names(ready)[6:9]


mo0 <- lmer(TY ~ (1|Plant), data=ready[which(ready$Year == "y2011"),])
y11 <- coef(mo0)$Plant
mo0 <- lmer(TY~(1|Plant), data=ready[which(ready$Year == "y2012"),])
y12 <- coef(mo0)$Plant
mo0 <- lmer(TY~(1|Plant), data=ready[which(ready$Year == "y2013"),])
y13 <- coef(mo0)$Plant

layout(matrix(1:4,2,2))
v <- which(rownames(y11) %in% rownames(y12))
gg <- cor(y11[v,1],y12[,1])
plot(y11[v,1],y12[,1], xlab="2011",ylab="2012", col="blue")
legend("topleft", legend=paste("cor=",gg))

v <- which(rownames(y11) %in% rownames(y13))
gg <- cor(y11[v,1],y13[,1])
plot(y11[v,1],y13[,1], xlab="2011",ylab="2013",col="red")
legend("topleft", legend=paste("cor=",gg))

v <- which(rownames(y12) %in% rownames(y13))
gg <- cor(y12[v,1],y13[,1])
plot(y12[v,1],y13[,1],xlab="2012",ylab="2013",col="green")
legend("topleft", legend=paste("cor=",gg))

###############################
##############################
# comparison with old models
data <- data.tre

for(i in 2:8){
  nb <- i
  data$iblock <- 0
  head(data)
  #################################### VERTICAL BLOCKS, i.e. more blocks based on 54 cols
  div <-round(max(data$Col)/nb)
  s <- seq(1,max(data$Col),by=div) # where each block start
  for(i in 1:length(s)){ # loop to add the new variable block
    r <- seq(s[i],s[i]+div-1)
    v <- which(data$Col %in% r) # adding a block designation based location in the field, all plants in row1-3 assigned to block 1 i.e.
    data$iblock[v]=i
  }
  
  mor <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year) + 
                (1|iblock)  + 
                (1|iblock:Year), data=data, REML=F)
  print(paste(i,"blocks"))
  print(paste("error variance=",attr(summary(mor)$varcor, 'sc')^2))
  print(summary(mor)$AICtab)
}
mox2 <- update(mox, REML=F)
print(paste("error variance=",attr(summary(mox2)$varcor, 'sc')^2))
print(summary(mox2)$AICtab)

layout(matrix(1:4,2,2))
plot(1,1)
plot(fitted(mo0),resid(mo0),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model",ylim=c(-300,700), main="Residual vs Fitted values for dummy model")
plot(x=fitted(mor),y=resid(mor),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", ylim=c(-300,700), main="Residual vs Fitted values for G*Y + block")
plot(x=fitted(mo2),y=resid(mo2),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", ylim=c(-300,700), main="Residual vs Fitted values for G*Y model")
plot(x=fitted(mo3),y=resid(mo3),col="black", pch=4, xlab="Fitted values", ylab="Residual in the model", ylim=c(-300,700), main="Residual vs Fitted values for G*Y + NNA")

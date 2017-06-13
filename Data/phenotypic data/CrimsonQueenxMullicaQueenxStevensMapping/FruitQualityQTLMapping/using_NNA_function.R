rm(list=ls())
#setwd("~/PAPERS TO DO/Brandon analysis")
setwd("C:/Users/zalapalab/Dropbox/Brandon")
dir()
library(lmerTest)
library(lattice)
# loading the data
#datar <- read.csv("CQxMQ_PhenotypicData_EduardoUsed.csv",header=T)
#datar <- read.csv("CQxMQ_PhenotypicData_EduardoUsed.csv",header=T)
datar <- read.csv("CNJ02_PlotData.csv",header=T)

head(datar)
tail(datar)

######################################################
# use NNA function to get the nearest neighbour effect
vars <- names(datar)[6:13]
for(i in 1:length(vars)){
  ho <- vars[i]
  ready <- NNA(datar=datar, plantid="Plant", repetition="Month", covar="Year", response=paste(ho))
  print(head(ready))
}

# heatplots to asses spatial variation
layout(matrix(1,1,1))
library(lattice)
library(lmerTest)
levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main="2011 residuals distribution", data=ready[which(ready[,4] == "y2011"),])
levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main="2012 residuals distribution", data=ready[which(ready[,4] == "y2012"),])
levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main="2013 residuals distribution", data=ready[which(ready[,4] == "y2013"),])
# models to compare
mo0 <- lmer(TY~(1|Plant), data=ready)
mo2 <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year), data=ready)
mo3 <- lmer(TY~(1|Plant) + (1|Year) + (1|Plant:Year) + (1|NS) + (1|Year:NS), data=ready)
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
levelplot(resmo3 ~Col*Row,col.regions = my_palette, at=col_breaks, main="residuals distribution for adjusted model", data=ready)
levelplot(resmo0 ~Col*Row,col.regions = my_palette, at=col_breaks, main="residuals distribution for unadjusted model", data=ready)

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

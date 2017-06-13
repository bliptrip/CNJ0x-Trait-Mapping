dir()
CQxMQ<-read.csv("CNJ02_PlotData.csv")
str(CQxMQ)
head(CQxMQ)

Tacy<-lm(CQxMQ$Tacy~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(Tacy)
summary(Tacy)

Brix<-lm(CQxMQ$Brix~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(Brix)
summary(Brix)

TA<-lm(CQxMQ$TA~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(TA)
summary(TA)

PAC<-lm(CQxMQ$PAC~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(PAC)
summary(PAC)

TY<-lm(CQxMQ$TY~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(TY)
summary(TY)

SFY<-lm(CQxMQ$SFY~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(SFY)
summary(SFY)

PFR<-lm(CQxMQ$PFR~CQxMQ$Month*CQxMQ$Year, data=CQxMQ)
anova(PFR)
summary(PFR)

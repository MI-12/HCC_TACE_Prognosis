
'This is the source R code of the submission, which including the:
    1. Univariate multinomial logistic regression (MLR) for radiomics feature selection
    2. Multivariate MLR to create the HCC radiomics signature
    3. Patient characteristics evaluation and significant characteristics identification
    4. Preoperative and postoperative models construction
    5. Calibration plot, decision curve analysis, and clinical impact curve in this study'

"Notice: When it comes to the path that users need to read their own data, as well as the variables defined,
         they are all represented in this file by'...'. Please modify these to your own."

rm(list = ls())
gc()
library("rJava")
library("xlsxjars")
library("xlsx")
library("survival")
library("glmnet")
library("SDMTools")
library("pROC")
library("readxl")
library("nnet")
library("mlogit")
library("survminer")
library("survivalROC")
library("rms")
library("survival")
library("survIDINRI")
library("rmda")

workbook<-"Load your data.xlsx"

Sta = read_excel(workbook,Number)


#Your Training and Validation datasets
Train<-data.frame(Sta['Validation data',])
Verify<-data.frame(Sta['Training data',])

#Calculate the ROC of LD level prediction for each feature
roc.feature<-vector(length=ncol('Number of Your Training Features'))
for (i in 1:ncol('Number of Your Training Features'))
  
  mult.cere<-multinom(Train$LD.score~TrainingFeature[,i],data=Train)
  N<-fitted(mult.cere)
  S<-apply(N, 1, function(x){max(x)})
  RocofTrain<-roc(TrainVariable1$ExpModRes,S)
  rocis<-RocofTrain$auc
  roci1<-sapply(rocis,as.numeric)
  roc.feature[i]<-roci1
  
mult.cere<-multinom(TrainVariable1$ExpModRes~'The Significant Feature'+'The Significant Feature'+'The Significant Feature',data=Train)
N<-fitted(mult.cere)
S<-apply(N, 1, function(x){max(x)})
RocofTrain<-roc(TrainVariable1$ExpModRes,S)
RocofTrain$auc
  
plot(RocofTrain, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
    grid.col=c("green", "red"), max.auc.polygon=TRUE,
    auc.polygon.col="skyblue", print.thres=TRUE)
  

#Prognosis models construction
TrainVariable<-Train[('Your Training Variables')]

#size
Size<-TrainVariable$Size

#metastasis
Meta<-TrainVariable$Meta

#capsule
Envolope<-TrainVariable$capsule

#radiomicssignature
Signature<-TrainVariable$Signature

#LD Score
LD<-TrainVariable1$LD.Score

#age
Age<-TrainVariable1$Age

#Gender
Gender<-TrainVariable1$Gender

#AFP
AFP<-as.numeric(TrainVariable1$AFP)

#Karnofsky Score
Karno<-TrainVariable1$Karnofsky

#Child-Puge
Child<-TrainVariable1$Child.Puge

dataTrain<-data.frame(Size,Meta,Envolope,Signa,LD,Age,Gender,AFP,Karno,Child)
ddlist<-datadist(dataTrain)
options(datadist='ddlist')

TTP<-TrainVariable1$TTP
Status<-TrainVariable1$status
S<-Surv(TTP,Status)


#univariate Cox regression analysis to obtain the significant patient characteristics
Uni <- cph(S ~ 'Each of the characteristics' x=T, y=T,surv=TRUE,data = dataTrain,time.inc='The cut-off time')
surv<-Survival(f)
surv<-function(x) surv(365,lp=x) #'The cut-off time = 365 days'
rcorrcens(Surv(TTP,Status)~predict(Uni,data = TrainVariable))


#For preoperative prognosis model
Pref <- cph(S ~  Size+Meta+Envolope, x=T, y=T,surv=TRUE,data = dataTrain,time.inc='The cut-off time')
surv<-Survival(Pref)
surv1<-function(x) surv(180,lp=x) ##'The cut-off time = half a year'
surv2<-function(x)surv(365,lp=x) ##'The cut-off time = one year'
nom<-nomogram(Pref,fun=list(surv1,surv2),lp=F,funlabel=c('6-months TTP','1-year TTP'),
              fun.at = c(.05,seq(.1,.5,by = .05),.9))
plot(nom,xfrac = .25,label.every = 2,lmgp = .2,cex.axis=1.0)
validate(Postf,method="boot",B=1000,dxy=T)
s<-rcorrcens(Surv(TTP,Status)~predict(Pref),data = TrainVariable)

#For postoperative prognosis model
Postf <- cph(S ~  Size+Meta+Envolope+LD, x=T, y=T,surv=TRUE,data = dataTrain,time.inc='The cut-off time')
surv<-Survival(Postf)
surv1<-function(x) surv(180,lp=x) ##'The cut-off time = half a year'
surv2<-function(x)surv(365,lp=x) ##'The cut-off time = one year'
nom<-nomogram(Postf,fun=list(surv1,surv2),lp=F,funlabel=c('6-months TTP','1-year TTP'),
              fun.at = c(.05,seq(.1,.5,by = .05),.9))
plot(nom,xfrac = .25,label.every = 2,lmgp = .2,cex.axis=1.0)
validate(Postf,method="boot",B=1000,dxy=T)
s<-rcorrcens(Surv(TTP,Status)~predict(Postf),data = TrainVariable)



#Validation of the models
VerifyVariable<-Verify[('Your Verify Variables')]

#size
Size<-VerifyVariable$Size

#metastasis
Meta<-VerifyVariable$Metastasis

#envolope
Envolope<-VerifyVariable$capsule

#radiomics signature
Signa<-VerifyVariable$Signature

#LD score
LD<-VerifyVariable$LD.score

#age
Age<-VerifyVariable$Age

#Gender
Gender<-VerifyVariable$Gender

#AFP
AFP<-VerifyVariable$AFP

#Karnofsky score
Karno<-VerifyVariable$Karnofsky

#Child-Puge
Child<-VerifyVariable$Child.pugh
   
dverify <-data.frame(BCLC,Size,Meta,Envolope,Signa,LD,Age,Gender,AFP,Karno,Child)
ddlistverify<-datadist(dverify)
options(datadist='ddlistverify')
TTP<-VerifyVariable$TTP
Status<-VerifyVariable$status
fev4 <- cph(Surv(TTP, Status) ~ predict(Pref/Postf, newdata=dverify), x=T, y=T, surv=T, data=d4, time.inc='The cut-off time')
#validate(fev3, method="boot", B=1000, dxy=T)
rcorrcens(Surv(TTP, Status) ~ predict(Pref/Postf, newdata=dverify), data = dverify)

     
##Calibration Plot  
calPref<-calibrate(Pref,cmethod ='KM',method ="boot",u = 'Cut-off time',m = 50,B = 300)
calPostf<-calibrate(Postf,cmethod ='KM',method ="boot",u = 'Cut-off time',m = 50,B = 300)
par(mar = c(8,5,3,2),cex = 1.0)
plot(calPref,lwd=2,lty = 1,errbar.col = c(rgb(0,118,192,maxColorValue = 255)),xlim = c(0,1),ylim = c(0,1),
     xlab = "Nomogram Predicted Survival",ylab = "Actual Survival",cex.axis=1.9,cex.lab=1.9,
     col = c(rgb(192,98,83,maxColorValue = 255)))

lines(calPref[,c('mean.predicted','KM')],type = 'b',lwd = 2,col = c(rgb(0,0,255,maxColorValue = 255)),pch = 16)
mtext(" ")
box(lwd = 1)
abline(0,1,lty = 3,lwd = 2,col = c(rgb(0,118,192,maxColorValue = 255)))
lines(calPostf[,c('mean.predicted','KM')],type = 'b',lwd = 2,col = c(rgb(255,0,0,maxColorValue = 255)),pch=16)
title("Calibration plots of the models",cex.main = 2, font.main= 1, col.main= "black")


##Decision Curve Analysis
data<-dataTrain
PDLog<-(Train$TTP < 'Cut-off time')
PDLog<-as.numeric(PDLog)
glmData<-PDLog
glmData<-data.frame(glmData)
Preglmmodel = glm(PDLog~Size+Meta+Envolope,family=binomial(link="logit"),data=data)
glmData$pred<-predict(Preglmmodel,type="response")
glmData<-data.frame(glmData)
Pre=dca(data=glmData,outcome="glmData",predictors="pred",smooth="TRUE", xstop=0.90)

Postglmmodel = glm(PDLog~Size+Meta+Envolope+LD,family=binomial(link="logit"),data=data)
glmData$pred2 = predict(Postglmmodel,type="response")
Post=dca(data=glmData, outcome="glmData", predictors="pred2", smooth="TRUE", xstop=0.90)
dca(data=glmData, outcome="glmData", predictors=c("pred","pred2"),xstop=0.9) 
col=c("black","red","green","blue","cyan","orange")
plot(Pre$net.benefit.threshold, Pre$net.benefit.none, type = "l", lwd=2,
     xlim=c(0,.90), ylim=c(-.05, .90), xlab = "Threshold Probability",
     ylab = "Net Benefit", font.lab = 1,cex.lab =1.9,col.lab="black",cex.axis=1.9)
lines(Pre$net.benefit$threshold, Pre$net.benefit$none, type="l", col=col[1], lwd=2)
lines(Pre$net.benefit$threshold, Pre$net.benefit$all, type="l", col=col[6], lwd=2)
lines(Pre$net.benefit$threshold, Pre$net.benefit$pred, type="l", col=col[2], lwd=2)
lines(Post$net.benefit$threshold, Post$net.benefit$pred2, type="l", col=col[4], lwd=2)
title("Decision Curve Analysis",cex.main = 2, font.main= 1, col.main= "black")


#Clinical Impact Curve
PDLog<-(sample.matrixtrain$TTP < 'Cut-off time')
PDLog<-as.numeric(PDLog)
data<-data.frame(data,PDLog)
baseline.model <- decision_curve(PDLog~Size+Meta+Envolope,data = data,family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),confidence.intervals= 0.95,study.design = 'cohort',population.prevalence= 0.3)
plot_clinical_impact(baseline.model, population.size = 1000,cost.benefit.axis = F, n.cost.benefits= 8,col = c('red','blue'), confidence.intervals= T,legend.position = "none",font.lab = 1.9,cex.lab =1.9,col.lab="black",cex.axis=1.9,cex.main = 1.9)
title("Clinical Impact Curve",cex.main = 2, font.main= 1, col.main= "black")










#This uses the tutorial for running a boosted conditional logistic regression model
#It uses methodology described in Lessman et al. 2012 
#(https://www.sciencedirect.com/science/article/pii/S0377221711009714)
#which uses stacking to combine a boosted regression with a clogit to improve predictions 
#compared to a conventional conditional logistic regression

#load libraries
library(amt)
library(lubridate)
library(dismo)
library(survival)
library(raster)
library(sf)
library(tidyverse)
library(udpipe)
library(sp)


########################################
###Format data##########################
########################################
#Need get data in correct format, and create the false steps (see other code)
Data.Path <- "X"
setwd(Data.Path)
ISSFStepsGb<-subset(ISSFStepsG,ISSFStepsG$group!=2)
groups <- ISSFStepsGb[, c("group")]
groups2<-unique(groups)

table(groups2)
table(ISSFStepsG$group)

for(i in 1:length(groups2)) { 
  Dat_T<-subset(ISSFStepsGb,group==groups2[[i]][1])
  IDT<-table(Dat_T$id)
  Dat_T<-droplevels(Dat_T)
  NLCDT<-as.data.frame(table(Dat_T$nlcd2,Dat_T$case_))
  NLCDT2<-subset(NLCDT,NLCDT$Var2=="TRUE"&NLCDT$Freq>0)
  NLCDT2$nlcd2<-NLCDT2$Var1
  Count<-count(NLCDT2)
  if(Count$n>1) {
    
    Dat_T2<-merge(Dat_T,NLCDT2, by="nlcd2")  
    table(Dat_T2$Freq)
    Dat_T2<-droplevels(Dat_T2)
    table(Dat_T2$nlcd2,Dat_T2$case_)
    contrasts(Dat_T2$nlcd2) =contr.sum(as.numeric(Count$n))
tr4<-Dat_T2

tr4$case_[which(tr4$case_==FALSE)]=0
tr4$case_[which(tr4$case_==TRUE)]=1

#make a dataframe
tr5=as.data.frame(tr4)

#set gbm.x and y
gbmx=c(1,8)
gbmy=16

#Boosted Regression

#binomial family
gbm.fisher=gbm.step(tr5, gbmx, gbmy, tree.complexity = 1,
										learning.rate = 0.001, bag.fraction = 0.75, n.folds = 10,
										family = "bernoulli", n.trees = 50, max.trees = 5000)

#make the predictions
fisherpreds=predict(gbm.fisher, tr5,type="response")

#bind preds to original df
tr5<-cbind(tr5,fisherpreds)

#Stack boosted predictions with clogit model
#use predicted probabilities from the boosted model to predict outcome in a clogit model
#this is slight change to the amt wrapper for clogit, uses the same underlying function
#when using the amt wrapper, predict doesn't work-- so need to get the one directly from survival
m.bst = survival::clogit(case_ ~ fisherpreds + strata(set), tr5)
m.clog = survival::clogit(case_ ~ nlcd2*sl_ + strata(set), tr5)

#use predict function
preds=predict(m.bst,tr5,type="expected")
preds.clog=predict(m.clog,tr5,type="expected")

#bind results to dataframe
tr5=cbind(tr5,preds) 
tr5=cbind(tr5,preds.clog)

#This rounds to 1 or 0 to see if predicts correctly 
tr5$preds2<-round(10*tr5$preds)
tr5$preds.clog2<-round(10*tr5$preds.clog)
TPred<-table(tr5$case_,tr5$preds2)
TPred.clog<-table(tr5$case_,tr5$preds.clog2)

#Model comparison

#run a glm between preds and expected 
#highest estimate will indicate best model
glm1=glm(tr5$fisherpreds~tr5$case_) #boosting only
glm2=glm(tr5$preds~tr5$case_) #boosting plus clogit
glm3=glm(tr5$preds.clog~tr5$case_) #clogit only

#Save outputs
sum1T<-summary(glm1)$coefficients
sum2T<-summary(glm2)$coefficients
sum3T<-summary(glm3)$coefficients
sum3T<-summary(m.bst)$coefficients
SumName1 <- paste0("SumA",i,".csv")
SumName2 <- paste0("SumB",i,".csv")
SumName3 <- paste0("SumC",i,".csv")
LandName<- paste0("LT",i,".csv")
IDName<- paste0("IDT",i,".csv")
PredName<- paste0("Pred_",i,".csv")
PredCName<- paste0("PredC_",i,".csv")
write.csv(sum1T, file=SumName1)
write.csv(sum2T, file=SumName2)
write.csv(sum3T, file=SumName3)
write.csv(NLCDT, file=LandName)
write.csv(IDT, file=IDName)
write.csv(TPred, file=PredName)
  write.csv(TPred.clog, file=PredCName)


sumclogT<-summary(m.clog)$coefficients
SumNameC <- paste0("SumClog",i,".csv")
write.csv(sumclogT, file=SumNameC)

#compare MSE and R2 between each of the three models
MSE.bst=sum((tr5$case_-tr5$fisherpreds)^2)/nrow(tr5)
MSE.bst
MSE.bst.clogit=sum((tr5$case_-tr5$preds)^2)/nrow(tr5)
MSE.bst.clogit
MSE.clog=sum(na.omit((tr5$case_-tr5$preds.clog)^2))/nrow(tr5)
MSE.clog

End1<-rbind(MSE.bst,MSE.bst.clogit)
End1<-rbind(End1,MSE.clog)

End1<-as.data.frame(End1)
EndName<- paste0("D",i,".csv")
write.csv(End1, file=EndName)}

  else {}
}

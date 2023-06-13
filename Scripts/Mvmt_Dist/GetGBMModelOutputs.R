##Identifying optimal parameters for gradient boosted models

#created using tutorials:
#https://bradleyboehmke.github.io/HOML/gbm.html
#https://rspatial.org/raster/sdm/9_sdm_brt.html

###################
##Process Outline##
###################

#input: pigsums.csv

################
##Script Setup##
################

#Load libraries
library(tidyverse)
library(dplyr)
library(dismo)
library(gbm)
library(gpboost)
library(job)
library(SHAPforxgboost)
library(caret)
library(ggplot2)
library(glmnet)

#set working directories
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt"
#remote:
#home<-"C:\\Users\\kayleigh.chalkowski\\OneDrive - USDA\\Projects\\StatPigMvmt\\"
setwd(home)
#gpbloaded=gpb.load(filename = paste0(home,"/Outputs/ModelComparison/gpbmod.json"))
  
#read data
pigsums<-read.csv(paste0(home,"/Data/pigsums.csv"))
pigswsite<-read.csv(paste0(home,"/Data/PigsDisplaceJan.csv"))

pigswsite[pigswsite$id=="33348"&pigswsite$region=="Mitchell",]$id<-"33348.m"
pigswsite[pigswsite$id=="33348"&pigswsite$region=="Nate Snow",]$id<-"33348.ns"
pigswsite[pigswsite$id=="33356"&pigswsite$region=="Mitchell",]$id<-"33356.m"
pigswsite[pigswsite$id=="33356"&pigswsite$region=="Nate Snow",]$id<-"33356.ns"
pigswsite[pigswsite$id=="33358"&pigswsite$region=="Mitchell",]$id<-"33358.m"
pigswsite[pigswsite$id=="33358"&pigswsite$region=="Nate Snow",]$id<-"33358.ns"

pigsums$id<-as.character(pigsums$id)
pigsites=pigswsite[,c(2,3)]
pigsums2=left_join(pigsums,pigsites,by="id")
pigsums2=pigsums2[!duplicated(pigsums2),]

#region got duplicated
pigsums2<-pigsums2[,-1]
pigsums2<-pigsums2[,-c(2)]
colnames(pigsums2)[ncol(pigsums2)]<-"region"

#make sure all correct classes
pigsums2[,c(2,3,6:48)] <- lapply(pigsums2[,c(2,3,6:48)],as.numeric)
pigsums2[,c(1,4,5,49)] <- lapply(pigsums2[,c(1,4,5,49)],as.factor)

#33348
#33356
#33358
#pigsums2[pigsums2$id=="33348",]
#pigswsite[pigswsite$id=="33348",]$Random
#pigswsite[pigswsite$id=="33356",]$Random
#tail(pigswsite[pigswsite$id=="33358",]$region)
#nrow(pigsums)
#length(unique(pigsums$id))
#length(unique(pigswsite$id))
#pigsums[pigsums$id=="33358",]


#############
##Tidy Data##
#############

#source functions
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/DropUnimportantVars.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/Optimize_Alpha_Lasso.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/run_Lasso_function.R"))

#repname="17FEB23_Runs"
###############################
##Hyperparameter optimization##
###############################
library(job)
pigsums=pigsums2

job::job({
#for testing
ko_t=10 #inner k-fold cross validations
ki_t=3 #inner k-fold cross validations
#response="sl_"
#pigsums
X_vec.start=c(3,4,6:13,18,21:35,38:44,46:48)
split_type="region" #randomly splits the data into test/training sets according to ko/ki
#split_type="random" #splits the data according to 'region', which is the study name

#######res objects contain:
#res[[1]] - gbm model object
#res[[2]] - mean RMSE of model
#res[[3]] - best model hyperparameters

#source functions
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/DropUnimportantVars.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/Optimize_Alpha_Lasso.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/run_Lasso_function.R"))

#make cutoff set for sigma models
cutoff_sl=pigsums[pigsums$count>150,]

#Run models
res.sl=Run.GBM.Model(pigsums,"sl_",X_vec.start,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

#Export CV stats from first run
write.csv(res.sl[[2]],"Outputs/05APR23_Runs/Region/sl_meanRMSE.csv")
write.csv(res.sl[[3]],"Outputs/05APR23_Runs/Region/sl_bestmodelparams.csv")
write.csv(res.sl[[4]],"Outputs/05APR23_Runs/Region/sl_meanR2.csv")
saveRDS(res.sl[[5]],"Outputs/05APR23_Runs/Region/sl_preds.rds")
saveRDS(res.sl[[6]],"Outputs/05APR23_Runs/Region/sl_testobs.rds")

res.sigma.sl=Run.GBM.Model(cutoff_sl,"sigma_sl",X_vec.start,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)

write.csv(res.sigma.sl[[2]],"Outputs/05APR23_Runs/Region/sigma.sl_meanRMSE.csv")
write.csv(res.sigma.sl[[3]],"Outputs/05APR23_Runs/Region/sigma.sl_bestmodelparams.csv")
write.csv(res.sigma.sl[[4]],"Outputs/05APR23_Runs/Region/sigma.sl_meanR2.csv")
saveRDS(res.sigma.sl[[5]],"Outputs/05APR23_Runs/Region/sigma.sl_preds.rds")
saveRDS(res.sigma.sl[[6]],"Outputs/05APR23_Runs/Region/sigma.sl_testobs.rds")

res.disp=Run.GBM.Model(pigsums,"displacement",X_vec.start,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

write.csv(res.disp[[2]],"Outputs/05APR23_Runs/Region/disp_meanRMSE.csv")
write.csv(res.disp[[3]],"Outputs/05APR23_Runs/Region/disp_bestmodelparams.csv")
write.csv(res.disp[[4]],"Outputs/05APR23_Runs/Region/disp_meanR2.csv")
saveRDS(res.disp[[5]],"Outputs/05APR23_Runs/Region/disp_preds.rds")
saveRDS(res.disp[[6]],"Outputs/05APR23_Runs/Region/disp_testobs.rds")

res.sigma.disp=Run.GBM.Model(cutoff_sl,"sigma_disp",X_vec.start,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)

write.csv(res.sigma.disp[[2]],"Outputs/05APR23_Runs/Region/sigma.disp_meanRMSE.csv")
write.csv(res.sigma.disp[[3]],"Outputs/05APR23_Runs/Region/sigma.disp_bestmodelparams.csv")
write.csv(res.sigma.disp[[4]],"Outputs/05APR23_Runs/Region/sigma.disp_meanR2.csv")
saveRDS(res.sigma.disp[[5]],"Outputs/05APR23_Runs/Region/sigma.disp_preds.rds")
saveRDS(res.sigma.disp[[6]],"Outputs/05APR23_Runs/Region/sigma.disp_testobs.rds")

res.tenavg=Run.GBM.Model(pigsums,"tenavg",X_vec.start,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

write.csv(res.tenavg[[2]],"Outputs/05APR23_Runs/Region/tenavg_meanRMSE.csv")
write.csv(res.tenavg[[3]],"Outputs/05APR23_Runs/Region/tenavg_bestmodelparams.csv")
write.csv(res.tenavg[[4]],"Outputs/05APR23_Runs/Region/tenavg_meanR2.csv")
saveRDS(res.tenavg[[5]],"Outputs/05APR23_Runs/Region/tenavg_preds.rds")
saveRDS(res.tenavg[[6]],"Outputs/05APR23_Runs/Region/tenavg_testobs.rds")

})

#Get x_vecs without unimportant variables
#Remove variables that had less than 0.1 influence
X_vec.2.sl=drop.unimportant(res.sl[[1]],X_vec.start,0.1)
X_vec.2.sigma.sl=drop.unimportant(res.sigma.sl[[1]],X_vec.start,0.1)
X_vec.2.disp=drop.unimportant(res.disp[[1]],X_vec.start,0.1)
X_vec.2.sigma.disp=drop.unimportant(res.sigma.disp[[1]],X_vec.start,0.1)
X_vec.2.tenavg=drop.unimportant(res.tenavg[[1]],X_vec.start,0.1)

#Save new vectors to file
write.csv(X_vec.2.sl,"Outputs/05APR23_Runs/Region/X_vec_01Drop/Xsl.csv")
write.csv(X_vec.2.sigma.sl,"Outputs/05APR23_Runs/Region/X_vec_01Drop/Xsigma.sl.csv")
write.csv(X_vec.2.disp,"Outputs/05APR23_Runs/Region/X_vec_01Drop/Xdisp.csv")
write.csv(X_vec.2.sigma.disp,"Outputs/05APR23_Runs/Region/X_vec_01Drop/Xsigma.disp.csv")
write.csv(X_vec.2.tenavg,"Outputs/05APR23_Runs/Region/X_vec_01Drop/Xtenavg.csv")

#Run models without dropped variables
job::job({
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/DropUnimportantVars.R"))
  cutoff_sl=pigsums[pigsums$count>150,]
  
res.sl=Run.GBM.Model(pigsums,"sl_",X_vec.2.sl,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_vec.2.sigma.sl,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
res.disp=Run.GBM.Model(pigsums,"displacement",X_vec.2.disp,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_vec.2.sigma.disp,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
res.tenavg=Run.GBM.Model(pigsums,"tenavg",X_vec.2.tenavg,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

#Export CV stats
write.csv(res.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sl_meanRMSE.csv")
write.csv(res.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sl_bestmodelparams.csv")
write.csv(res.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sl_meanR2.csv")
saveRDS(res.sl[[5]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sl_preds.rds")
saveRDS(res.sl[[6]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sl_testobs.rds")

write.csv(res.sigma.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.sl_meanRMSE.csv")
write.csv(res.sigma.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.sl_bestmodelparams.csv")
write.csv(res.sigma.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.sl_meanR2.csv")
saveRDS(res.sigma.sl[[5]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.sl_preds.rds")
saveRDS(res.sigma.sl[[6]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.sl_testobs.rds")

write.csv(res.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_01Drop/disp_meanRMSE.csv")
write.csv(res.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_01Drop/disp_bestmodelparams.csv")
write.csv(res.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_01Drop/disp_meanR2.csv")
saveRDS(res.disp[[5]],"Outputs/05APR23_Runs/Region/GBM_01Drop/disp_preds.rds")
saveRDS(res.disp[[6]],"Outputs/05APR23_Runs/Region/GBM_01Drop/disp_testobs.rds")

write.csv(res.sigma.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.disp_meanRMSE.csv")
write.csv(res.sigma.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.disp_bestmodelparams.csv")
write.csv(res.sigma.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.disp_meanR2.csv")
saveRDS(res.sigma.disp[[5]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.disp_preds.rds")
saveRDS(res.sigma.disp[[6]],"Outputs/05APR23_Runs/Region/GBM_01Drop/sigma.disp_testobs.rds")

write.csv(res.tenavg[[2]],"Outputs/05APR23_Runs/Region/GBM_01Drop/tenavg_meanRMSE.csv")
write.csv(res.tenavg[[3]],"Outputs/05APR23_Runs/Region/GBM_01Drop/tenavg_bestmodelparams.csv")
write.csv(res.tenavg[[4]],"Outputs/05APR23_Runs/Region/GBM_01Drop/tenavg_meanR2.csv")
saveRDS(res.tenavg[[5]],"Outputs/05APR23_Runs/Region/GBM_01Drop/tenavg_preds.rds")
saveRDS(res.tenavg[[6]],"Outputs/05APR23_Runs/Region/GBM_01Drop/tenavg_testobs.rds")

})

#Run lasso function for each model
#get list of remaining important variables to use
nrep=100
typemeasure="mae"

#some tidying of x, cant have characters
x=pigsums[,X_vec.start]
x$sexf=NA
x[x$sex=="Male",]$sexf=1
x[x$sex=="Female",]$sexf=0
x=x[,-2]
#class(x$sexf)
#make alpha search grid
grid=seq(0.01,1,by=0.01)

job::job({
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/run_Lasso_function.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/Optimize_Alpha_Lasso.R"))
  
#Run lasso on each model
keepVars.sl=runLassofunc(x,pigsums$sl_,"poisson",grid,typemeasure,nrep)
X_vec.sl.lasso=which(colnames(pigsums)%in%names(keepVars.sl))

keepVars.sigma.sl=runLassofunc(x,pigsums$sigma_sl,"poisson",grid,typemeasure,nrep)
X_vec.sigma.sl.lasso=which(colnames(pigsums)%in%names(keepVars.sigma.sl))

keepVars.disp=runLassofunc(x,pigsums$displacement,"poisson",grid,typemeasure,nrep)
X_vec.disp.lasso=which(colnames(pigsums)%in%names(keepVars.disp))

keepVars.sigma.disp=runLassofunc(x,pigsums$sigma_disp,"poisson",grid,typemeasure,nrep)
X_vec.sigma.disp.lasso=which(colnames(pigsums)%in%names(keepVars.sigma.disp))

keepVars.tenavg=runLassofunc(x,pigsums$tenavg,"poisson",grid,typemeasure,nrep)
X_vec.tenavg.lasso=which(colnames(pigsums)%in%names(keepVars.tenavg))

#Write out the post lasso X vecs
write.csv(X_vec.sl.lasso,"Outputs/05APR23_Runs/Region/X_vec_postlasso/Xsl_lasso.csv")
write.csv(X_vec.sigma.sl.lasso,"Outputs/05APR23_Runs/Region/X_vec_postlasso/Xsigmasl_lasso.csv")
write.csv(X_vec.disp.lasso,"Outputs/05APR23_Runs/Region/X_vec_postlasso/Xdisp_lasso.csv")
write.csv(X_vec.sigma.disp.lasso,"Outputs/05APR23_Runs/Region/X_vec_postlasso/Xsigmadisp_lasso.csv")
write.csv(X_vec.tenavg.lasso,"Outputs/05APR23_Runs/Region/X_vec_postlasso/Xtenavg_lasso.csv")

})


#Run lasso models
job::job({
split_type="region"

  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
  
  res.sl=Run.GBM.Model(pigsums,"sl_",X_vec.sl.lasso,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sl_meanRMSE.csv")
  write.csv(res.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sl_bestmodelparams.csv")
  write.csv(res.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sl_meanR2.csv")
  saveRDS(res.sl[[5]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sl_preds.rds")
  saveRDS(res.sl[[6]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sl_testobs.rds")
  
  res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_vec.sigma.sl.lasso,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.sl_meanRMSE.csv")
  write.csv(res.sigma.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.sl_bestmodelparams.csv")
  write.csv(res.sigma.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.sl_meanR2.csv")
  saveRDS(res.sigma.sl[[5]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.sl_preds.rds")
  saveRDS(res.sigma.sl[[6]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.sl_testobs.rds")
  
  res.disp=Run.GBM.Model(pigsums,"displacement",X_vec.disp.lasso,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_Lasso/disp_meanRMSE.csv")
  write.csv(res.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_Lasso/disp_bestmodelparams.csv")
  write.csv(res.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_Lasso/disp_meanR2.csv")
  saveRDS(res.disp[[5]],"Outputs/05APR23_Runs/Region/GBM_Lasso/disp_preds.rds")
  saveRDS(res.disp[[6]],"Outputs/05APR23_Runs/Region/GBM_Lasso/disp_testobs.rds")
  
  res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_vec.sigma.disp.lasso,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.disp_meanRMSE.csv")
  write.csv(res.sigma.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.disp_bestmodelparams.csv")
  write.csv(res.sigma.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.disp_meanR2.csv")
  saveRDS(res.sigma.disp[[5]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.disp_preds.rds")
  saveRDS(res.sigma.disp[[6]],"Outputs/05APR23_Runs/Region/GBM_Lasso/sigma.disp_testobs.rds")
  
  res.tenavg=Run.GBM.Model(pigsums,"tenavg",X_vec.tenavg.lasso,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.tenavg[[2]],"Outputs/05APR23_Runs/Region/GBM_Lasso/tenavg_meanRMSE.csv")
  write.csv(res.tenavg[[3]],"Outputs/05APR23_Runs/Region/GBM_Lasso/tenavg_bestmodelparams.csv")
  write.csv(res.tenavg[[4]],"Outputs/05APR23_Runs/Region/GBM_Lasso/tenavg_meanR2.csv")
  saveRDS(res.tenavg[[5]],"Outputs/05APR23_Runs/Region/GBM_Lasso/tenavg_preds.rds")
  saveRDS(res.tenavg[[6]],"Outputs/05APR23_Runs/Region/GBM_Lasso/tenavg_testobs.rds")
  
})


#Run null models
#gbm.null=rep(1,nrow(pigsums))
#igsums$gbm.null=rep(1,nrow(pigsums))
#X_null=ncol(pigsums)
#####figure out above, check make sure is good

job::job({
  pigsums$gbm.null=rep(1,nrow(pigsums))
  X_null=ncol(pigsums)
  split_type="region"
  
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
  
  res.sl=Run.GBM.Model(pigsums,"sl_",X_null,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_Null/sl_meanRMSE.csv")
  write.csv(res.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_Null/sl_bestmodelparams.csv")
  write.csv(res.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_Null/sl_meanR2.csv")
  
  res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_null,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.sl[[2]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.sl_meanRMSE.csv")
  write.csv(res.sigma.sl[[3]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.sl_bestmodelparams.csv")
  write.csv(res.sigma.sl[[4]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.sl_meanR2.csv")
  
  res.disp=Run.GBM.Model(pigsums,"displacement",X_null,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_Null/disp_meanRMSE.csv")
  write.csv(res.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_Null/disp_bestmodelparams.csv")
  write.csv(res.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_Null/disp_meanR2.csv")
  
  res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_null,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.disp[[2]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.disp_meanRMSE.csv")
  write.csv(res.sigma.disp[[3]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.disp_bestmodelparams.csv")
  write.csv(res.sigma.disp[[4]],"Outputs/05APR23_Runs/Region/GBM_Null/sigma.disp_meanR2.csv")
  
  res.tenavg=Run.GBM.Model(pigsums,"tenavg",X_null,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.tenavg[[2]],"Outputs/05APR23_Runs/Region/GBM_Null/tenavg_meanRMSE.csv")
  write.csv(res.tenavg[[3]],"Outputs/05APR23_Runs/Region/GBM_Null/tenavg_bestmodelparams.csv")
  write.csv(res.tenavg[[4]],"Outputs/05APR23_Runs/Region/GBM_Null/tenavg_meanR2.csv")
  
})

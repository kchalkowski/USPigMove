library(raster)
library(amt)
library(sf)
library(tidyverse)
library(udpipe)
library(lubridate)
library(sp)
library(dismo)
library(survival)
        
 # Start with clean data set and in the same coordinate system and same data format

MSPigsStep2 <- st_as_sf(MSPigsSteps, coords = c("x1_", "y1_"), crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
MSPigsStep3<-st_transform(MSPigsStep2, crs="EPSG:4326")
MSPigsStep3$lon.x <- sapply(MSPigsStep3$geometry, "[", 1)
MSPigsStep3$lat.x <- sapply(MSPigsStep3$geometry, "[", 2)
MSPigsStep3$date<-as.POSIXct(MSPigsStep3$t1_, format= "%Y-%m-%d %H:%M")
table(MSPigsStep3$date)
MSPigsStepOct <- MSPigsStep3[, c("sl_","lon.x","lat.x","id","period","sex","region","State","date")]

# Remove Judas and Canada as outside this study
PigsDataGam2<-subset(PigsDataGam,is.na(PigsDataGam$lon.x)==FALSE)
PigsDataGam2<-subset(PigsDataGam2,PigsDataGam2$region!="Judas")
PigsDataGam2<-subset(PigsDataGam2,PigsDataGam2$State!="Canada")
table(PigsDataGam2$region)
range(PigsDataGam$date)

#Adding in last data addition - MS Pigs
PigsTotal2 <- st_as_sf(PigsDataGam2, coords = c("lon.x", "lat.x"), crs="EPSG:4326")
PigsTotal2$lon.x <- sapply(PigsTotal2$geometry, "[", 1)
PigsTotal2$lat.x <- sapply(PigsTotal2$geometry, "[", 2)
PigsTotal2$date<-as.POSIXct(PigsTotal2$t1_, format= "%m/%d/%Y %H:%M")
PigsTotal3 <- PigsTotal2[, c("sl_","lon.x","lat.x","id","period","sex","region","State","date")]
PigsTotalOct<-rbind(PigsTotal3,MSPigsStepOct)
write.csv(PigsTotalOct,"PigsTotalOct14.csv")
table(PigsTotalOct$State)

#Setting Up Data for use with AMT package
PigsTotalOct2<-st_transform(PigsTotalOct, crs="+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
st_crs(PigsTotalOct2)
points<-PigsTotalOct2
points$x_ <- sapply(points$geometry, "[", 1)
points$y_ <- sapply(points$geometry, "[", 2)
points$t_<-as.POSIXct(points$date)
table(is.na(points$t_))
points$animalid<-points$id
points$hour<-hour(points$t_)
table(points$hour)
Steps <- points[, c("x_","y_","t_","animalid","hour")]

#iSSF setting up steps by every 6 hours

Alln <- Steps %>%
  group_by(animalid) %>%
  summarise(
    n=n(),
  ) %>%
  as.data.frame()
ids <- Alln[, c("animalid")]

Steps <- points[, c("x_","y_","t_","animalid","hour")]

full<-data.frame()
Data.Path <- "X"
setwd(Data.Path)
for(i in 1:length(ids)) { 
  dat_1<-subset(Steps, animalid==ids[[i]][1])
  dat_1$t_<-as.POSIXct(dat_1$t_, format="%Y-%m-%d %H:%M:%S")
  dat_2<-subset(dat_1,is.na(dat_1$t_)=="FALSE")
  Dat2 <- dat_2[, c("x_","y_","t_","hour")]
  tr1<-amt::make_track(Dat2,x_,y_,t_, hour=hour,crs=crs)
  tr2 <- tr1 %>% track_resample(rate = hours(6), tolerance = minutes(60))
  trp1<-subset(tr2,hour=="23"|hour=="22"|hour=="21"|hour=="20"|hour=="19"|hour=="18")
  trp2<-subset(tr2,hour=="17"|hour=="16"|hour=="15"|hour=="14"|hour=="13"|hour=="12")
  trp3<-subset(tr2,hour=="11"|hour=="10"|hour=="9"|hour=="8"|hour=="7"|hour=="6")
  trp4<-subset(tr2,hour=="0"|hour=="1"|hour=="2"|hour=="3"|hour=="4"|hour=="5")
  #round1
  if(nrow(trp1)>1) {
    steps1<-trp1%>%steps(,diff_time_units="hours")
        trp1s<-steps1
        trp1s<-trp1s%>%mutate(id=ids[i])
        trp1s<-trp1s%>%mutate(period=1)
        full<-rbind(full,trp1s)
      }
  else {}
  #Round2
  if(nrow(trp2)>1) {
    steps1<-trp2%>%steps(,diff_time_units="hours")
        trp2s<-steps1
        trp2s<-trp2s%>%mutate(id=ids[i])
        trp2s<-trp2s%>%mutate(period=2)
        full<-rbind(full,trp2s)
      }
  else {}
  #round 3
  if(nrow(trp3)>1) {
    steps1<-trp3%>%steps(,diff_time_units="hours")
        trp3s<-steps1
        trp3s<-trp3s%>%mutate(id=ids[i])
        trp3s<-trp3s%>%mutate(period=3)
        full<-rbind(full,trp3s)
      }
  else {}
  #Round 4
  if(nrow(trp4)>1) {
    steps1<-trp4%>%steps(,diff_time_units="hours")
        trp4s<-steps1
        trp4s<-trp4s%>%mutate(id=ids[i])
        trp4s<-trp4s%>%mutate(period=4)
        full<-rbind(full,trp4s)
      }
  else {}
}
range(full$sl_)

#Subset hour intervals between 21-27 hours
FullSteps<-subset(full, full$dt_<27.0001 & full$dt_>20.9999)
range(FullSteps$sl_)
FullSteps2<-subset(FullSteps, FullSteps$sl_<25000)
write.csv(FullSteps2, "FullSteps2ISSF.csv")
hist(as.numeric(FullSteps2$dt_))
#hist(FullSteps2$sl_)
range(FullSteps2$sl_)

#Check pig days per ID
AllFull <- FullSteps %>%
  group_by(id) %>%
  summarise(
    n=n(),
  ) %>%
  as.data.frame()


#ADD 15 random steps
LCD <- raster("C:/Downloads/united_states_2015_v2/USA_NALCMS_2015_v2_land_cover_30m/USA_NALCMS_2015_v2_land_cover_30m.tif")

FullSteps3<-FullSteps2%>%random_steps(n=15)%>%extract_covariates(LCD)

#Rename
FullSteps3$nlcd<-FullSteps3$USA_NALCMS_2015_v2_land_cover_30m

FullSteps3$set<-unique_identifier(FullSteps3, fields = c("id","step_id_","period"), start_from = 1)
FullSteps3$nlcd<-as.factor(FullSteps3$nlcd)
FullSteps3$year<-year(FullSteps3$t1_)
write.csv(FullSteps3, "FullStepsOctISSF.csv")


#Add state info with info file

PigsInfo  <-subset(PigsTotalOct,select=c(id,sex,region,State))
PigsInfo <- st_set_geometry(PigsInfo, NULL)
PigsInfo2<-unique(PigsInfo)

ISSFSteps<-merge(FullSteps3,PigsInfo2, by="id") #Make sure each pig has unique id incase multiple studies used the same numbers
write.csv(ISSFSteps, "ISSFStepsOCt17.csv")
ISSFSteps<-ISSFStepsOCt17

#Categorize to combine temperate and subtropical and water and wetland
ISSFStepsc<-ISSFSteps %>%
  mutate(nlcd2 = case_when(nlcd < 3 ~ 'needleleaf',
                             nlcd < 6 ~ 'broadleaf',
                           nlcd == 6 ~ 'mixed',
                           nlcd < 9 ~ 'shrub',
                           nlcd < 11 ~ 'grass',
                           nlcd == 14 ~ 'wetland',
                           nlcd == 15 ~ 'cropland',
                           nlcd == 16 ~ 'barren',
                           nlcd == 17 ~ 'urban',
                           nlcd == 18 ~ 'wetland'))
table(ISSFStepsc$nlcd,ISSFStepsc$nlcd2)

Ids2b<-as.data.frame(table(ISSFStepsc$id,ISSFStepsc$nlcd2,ISSFStepsc$case_))

#Remove 1 case that is misscoded
ISSFStepsG_9<-subset(ISSFStepsG,ISSFStepsG$group=="9")
ISSFStepsG_N9<-subset(ISSFStepsG,ISSFStepsG$group!="9")
ISSFStepsG_9<-subset(ISSFStepsG_9,ISSFStepsG_9$nlcd2!="broadleaf")
ISSFStepsG<-rbind(ISSFStepsG_N9,ISSFStepsG_9)
write.csv(ISSFStepsG,"ISSFStepsG.csv")

#Just group summary for table
groupsIDs <- ISSFStepsG[, c("id","State","region","sex")]
groupsIDs2<-unique(groupsIDs)
write.csv(groupsIDs2,"AllPigsIDs.csv")

#Start here if loading in ISSFStepsG.csv data set

groups <- ISSFStepsG[, c("group")]
groups2<-unique(groups)

#NO Step Length or Boosted
Data.Path <- "X"
setwd(Data.Path)
ISSFStepsG$nlcd2<-as.factor(ISSFStepsG$nlcd2)


for(i in 1:length(groups2)) { 
  Dat_T<-subset(ISSFStepsG,group==groups2[[i]][1])
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
    T1<-survival::clogit(case_ ~ nlcd2 + strata(set), data=Dat_T2)
    summary(T1)
    preds.clog=predict(T1,Dat_T2,type="expected")
    Dat_T2=cbind(Dat_T2,preds.clog)
    Dat_T2$pred2<-round(10*Dat_T2$preds.clog)
    TPred<-table(Dat_T2$case_,Dat_T2$pred2)
    sumT<-summary(T1)$coefficients
    AIC<-AIC(T1)
    SumName <- paste0("Sum_",i,".csv")
    LandName<- paste0("LT_",i,".csv")
    IDName<- paste0("IDT_",i,".csv")
    PredName<- paste0("Pred_",i,".csv")
    AICName<- paste0("AIC_",i,".csv")
    write.csv(sumT, file=SumName)
    write.csv(NLCDT, file=LandName)
    write.csv(IDT, file=IDName)
    write.csv(AIC, file=AICName)
    write.csv(TPred, file=PredName)}
  else {}
}

# Comparison with SL
Data.Path <- "X"
setwd(Data.Path)
ISSFStepsGb<-subset(ISSFStepsG,ISSFStepsG$group!=2) # only has access to wetlands
groups <- ISSFStepsGb[, c("group")]
groups2<-unique(groups)


for(i in 1:length(groups2)) { 
  Dat_T<-subset(ISSFStepsG,group==groups2[[i]][1])
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
    T1<-survival::clogit(case_ ~ nlcd2*sl_ + strata(set), data=Dat_T2)
    summary(T1)
    preds.clog=predict(T1,Dat_T2,type="expected")
    Dat_T2=cbind(Dat_T2,preds.clog)
    Dat_T2$pred2<-round(10*Dat_T2$preds.clog)
    TPred<-table(Dat_T2$case_,Dat_T2$pred2)
    sumT<-summary(T1)$coefficients
    AIC<-AIC(T1)
    SumName <- paste0("Sum_",i,".csv")
    LandName<- paste0("LT_",i,".csv")
    IDName<- paste0("IDT_",i,".csv")
    PredName<- paste0("Pred_",i,".csv")
    AICName<- paste0("AIC_",i,".csv")
    write.csv(sumT, file=SumName)
    write.csv(NLCDT, file=LandName)
    write.csv(IDT, file=IDName)
    write.csv(AIC, file=AICName)
    write.csv(TPred, file=PredName)}
  else {}
}


#Without SL - verification
Data.Path <- "X"
setwd(Data.Path)
# Remove groups made of 1 study 
table(ISSFStepsG$group,ISSFStepsG$region)
ISSFStepsGb<-subset(ISSFStepsG,ISSFStepsG$group!=2)
ISSFStepsGb<-subset(ISSFStepsGb,ISSFStepsGb$group!=5)

# Check groups
table(ISSFStepsGb$group)
table(ISSFStepsGb$nlcd,ISSFStepsGb$group, ISSFStepsGb$case_)

#Loop removing one study at a time for all groups
regions <- ISSFStepsGb[, c("region")]
region2<-unique(regions)
groups <- ISSFStepsGb[, c("group")]
groups2<-unique(groups)

for(j in 1:length(region2)) {
  Dat_Ta<-subset(ISSFStepsGb,region!=region2[[j]][1])
  for(i in 1:length(groups2)) {
    Dat_T<-subset(Dat_Ta,group==groups2[[i]][1])
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
      T1<-survival::clogit(case_ ~ nlcd2 + strata(set), data=Dat_T2)
      summary(T1)
      preds.clog=predict(T1,Dat_T2,type="expected")
      Dat_T2=cbind(Dat_T2,preds.clog)
      Dat_T2$pred2<-round(10*Dat_T2$preds.clog)
      TPred<-table(Dat_T2$case_,Dat_T2$pred2)
      sumT<-summary(T1)$coefficients
      SumName <- paste0("Sum_",j,"_",i,".csv")
      LandName<- paste0("LT_",j,"_",i,".csv")
      IDName<- paste0("IDT_",j,"_",i,".csv")
      PredName<- paste0("Pred_",j,"_",i,".csv")
      write.csv(sumT, file=SumName)
      write.csv(NLCDT, file=LandName)
      write.csv(IDT, file=IDName)
      write.csv(TPred, file=PredName)}
    else {}
  }
}

Table1<-table(ISSFStepsG$region,ISSFStepsG$group)
Table2<-table(ISSFStepsG$nlcd2,ISSFStepsG$group)
write.csv(Table1,"Group_studysite.csv")
write.csv(Table2,"Group_NLCD.csv")

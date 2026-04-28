library(MARSS)
library(ggplot2)
library(forecast)
library(dplyr)
library(lubridate)
library(mgcv)
library(tidyr)
library(corrplot)
library(reshape2)
library(corrplot)
likelihoods  <- read.csv("likelihoodsSens.csv")%>%
  rename(Oceanographic=oceanographic_index2,SMURF=base,'No Index'=no_smurf, 
         "RREAS North"=RREASN, "RREAS Coastwide"=RREAS)
recdevs<- read.csv("RecruitmentDeviationsSens.csv")%>%
  select(-X)%>%
  rename(Oceanographic=oceanographic_index2,SMURF=base, "No Index"=no_smurf, 
         "RREAS North"=RREASN, "RREAS Coastwide"=RREAS)

OCNMS<-read.csv('Data/JuvenileIndexDatasets/Estimated-YOY-trend-coast.csv')%>%
  rename(Index=grand.mean)%>%
  select(year, Index)%>%
  mutate(dataset="OCNMS Dive Survey")

RREAS_North<-read.csv('Data/JuvenileIndexDatasets/YOYGroundfishCPUEPerYr_North.csv')%>%
  rename(year=YEAR, Index=est)%>%
  select(year, Index)%>%
  mutate(dataset="Northern YOY Rockfish")

Oceanographic<-read.csv('Data/JuvenileIndexDatasets/OceanographicIndexV1.csv')%>%
  rename(Index=fit)%>%
  select(year, Index)%>%
  mutate(dataset="Oceanographic Index")

SMURF<-read.csv('Data/JuvenileIndexDatasets/index_forSS.csv')%>%
  rename(Index=obs)%>%
  select(year, Index)%>%
  mutate(dataset="SMURF Juvenile Survey")

RREAS_Yellowtail_Coast<-read.csv('Data/JuvenileIndexDatasets/ytail_coastwide_indices.csv')%>%
  rename(Index=est, year=YEAR)%>%
  select(year, Index)%>%
  mutate(dataset="Yellowtail RREAS Coastwide")

dat_long1<-bind_rows(OCNMS,RREAS_North)%>%
  bind_rows(Oceanographic)%>%
  bind_rows(SMURF)%>%
  bind_rows(RREAS_Yellowtail_Coast)

dat_long<-dat_long1%>%filter(year>2013)%>%
  group_by(dataset)%>%
  mutate(scaled_index=scale(Index))%>%
  rename(Dataset=dataset)

YOYIndex<-ggplot(data=dat_long%>%filter(year>2013)%>%
                   group_by(Dataset)%>%
                   mutate(scaled_index=scale(Index)),
                 aes(group=Dataset,col=Dataset,x=year,y=scaled_index))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  ylab("Standardized Index")+
  xlab("Year")+
  theme_classic()

YOYno2016<-ggplot(data=dat_long%>%filter(year>2016)%>%
                    group_by(Dataset)%>%
                    mutate(scaled_index=scale(Index)),
                  aes(group=Dataset,col=Dataset,x=year,y=scaled_index))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  ylab("Standardized Index (omitting 2016)")+
  xlab("Year")+
  theme_classic()

png("YOYindices.png",width=6,height=3,units="in",res=1200)
YOYIndex
dev.off()

png("YOYindicesNo2016.png",width=6,height=3,units="in",res=1200)
YOYno2016
dev.off()


##### Table functions ####

# convert likelihoods to difference from base (assumed to be in column 2)
tab<-likelihoods
like_rows <- grep("_like", tab$Label)
for (icol in ncol(tab):2) {
  tab[like_rows, icol] <- tab[like_rows, icol] - tab[like_rows, 2]
}

# round length values to 1 decimal place
tab[grep("L_at", tab$Label), -1] <-
  round(tab[grep("L_at", tab$Label), -1], 1)
tab$Label <- gsub("_like", "", tab$Label )


write.csv(tab, "like_table.csv")


ggplot(SMURF) +
  geom_point(aes(y=obs, x=year))+
  geom_line(aes(y=obs, x=year))+
  ylab("Relative Index")+
  geom_errorbar(aes(ymin=obs-se_log,ymax=obs+se_log, x=year))+
  ylab("Year")+
  theme_bw()

colnames(dat_long)<-c("year","Index","Dataset","scaled_index")
dat_wide<-data.frame(dat_long%>%pivot_wider(c(year),names_from = Dataset, values_from = scaled_index))

colnames(dat_wide)<-c("year","OCNMS", "RREAS North", "Oceanographic", "SMURF", "RREAS Coastwide")%>%
  left_join()
corr<-dat_wide%>%filter(year!=2014&year!=2015&year!=2020)%>%select(-year)

M = cor(corr)
corrplot(M, method = 'number')

rec<-recdevs%>%filter(Yr<2018)%>%select(Yr, 'No Index')%>%rename(year=Yr)%>%
  left_join(dat_wide)%>%
  select(-OCNMS)

corr2<-rec%>%filter(year>2013)%>%select(-year)                
N<-cor(corr2)

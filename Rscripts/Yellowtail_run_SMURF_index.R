
#########################################################################
### Run the combined SMURF and ODFW oceanography data to get an index of abundance
### Yellowtail rockfish assessment 2025
### Original script - Melissa Monk (SWFSC), edited by Alison Whitman (ODFW)
#########################################################################
# updated by A. Whitman 02/20/2025
# to use new oceanography variables from Megan Feddern (NMFS)

# updated on 2/20 to finalize candidate models 
# updated on 2/28 to include final oceanography variables (16d) and run final tables 
# updated on 4/24/2026 to re-run a model we missed!

rm(list = ls(all = TRUE))
graphics.off()

library(sdmTMB)
library(tmbstan)
library(ggeffects)
library(MuMIn)
library(here)
library(glue)
library(tidyr)
library(dplyr)
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(bayesplot)
library(grid)
library(devtools)
library(ggeffects)
library(tidybayes)
library(gridExtra)
library(fitdistrplus)

#species and area identifiers - eventually put in function
pacfinSpecies <- "YTRF"
speciesName <- "Yellowtail Rockfish"
modelArea = "oregon"
indexName <-  "SMURF"
modelName <- "full"

# loading helper functions 

dir<-file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/Raw Index Scripts")
setwd(dir)
list.files()
#source("helper_functions.R")
source("diagnostics.R")
source("do_diagnostics.R")
source("format_hkl_data.R")
source("format_index.R")
source("get_index.R")
source("match.f.R")
source("plot_betas.R")
source("plot_index.R")
source("refactor.R")

# load data
dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS")
setwd(dir)

#load("data_for_GLM.RData")
dat<-read.csv("combined_settlement_updatedocean_v2.csv")

# subset to species of interest  - unnecessary here
#dat<-dat[dat$Common_Name==speciesName,]

# set dir for full model
#dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/yellowtail_oregon_SMURF_full")
#setwd(dir)

# explore the data 

summary(dat)
names(dat)
str(dat)

#any(is.na(dat$temp_c))

#dat$CPUE<-dat$Counts/dat$Total_Drift_Effort

## note using their calculated settlement rate rather calculating my own ##

dat$CPUE<-dat$SFLA_settlement_rate

# add some categorical factors to explore 
dat$region<-ifelse(dat$Site %in% c("HH","RR"),"South","North")
dat$treatment<-ifelse(dat$Site %in% c("HH","CF"),"Comparison","Reserve")
dat$Date<-as.POSIXct(dat$Date, format = "%Y-%m-%d")
dat$month<-as.numeric(format(dat$Date,"%m"))
table(dat$month)
summary(dat$julian)
dat$season<-ifelse(dat$julian<196,"early","late") # used July 15 as the divider 
table(dat$season)

#Look at the data
pos <- dat[dat$SFLA_count>0,]

# 255 positives out of 1376 = 18% not bad! 

with(pos, table(Site)) # all pretty even across sites
with(pos, table(treatment)) 
with(pos, table(region)) 
with(pos, table(season))
with(pos, table(year)) # should exclude 2013, potentially 2015 or 2022 too... 

# notes from meeting with Meghan on 1/17 - try a model with complete cases to see if code works 
# for the dredge function (this could then be fed into sdmTMB, as normal)
# also need to try some standard filtering to address the singularity error I'm getting in sdmTMB
# (e.g. there are some years with low positives, there are some months too, so things to try)

ggplot(pos, aes(Site, CPUE)) + # 
  geom_boxplot()

ggplot(pos, aes(treatment, CPUE)) + 
  geom_boxplot()

ggplot(pos, aes(region, CPUE)) + # oh wow, higher in south...
  geom_boxplot()

ggplot(pos, aes(season, CPUE)) + # more in the early season 
  geom_boxplot()

# looking at new temp data covariates! 

ggplot(pos, aes(temp_index, CPUE)) + # primarily in the -1 to 0 on the standardized index, could bin these
  geom_point(alpha = 0.5)

ggplot(pos, aes(temp_c_mid, CPUE)) + # 8 - 10 degree range
  geom_point(alpha = 0.5)

ggplot(pos, aes(cdd_8, CPUE)) + # hmmm
  geom_point(alpha = 0.5)

ggplot(pos, aes(rolling8d, CPUE)) + # hmmm
  geom_point(alpha = 0.5)


with(pos, table(year, month)) # thin on samples in August and Sept
with(pos, table(month))
with(dat, table(month))

# fix or add anything 
dat <- dat %>%
  #rename(Year = year) %>%
  rename(Effort = Sampling.Interval) %>%
  mutate(logEffort = log(Effort)) %>% 
  #create temp bins for drill
  #mutate(temp_bin = cut(temp_c, breaks=c(6,7,8,9,10,11,12)))
  # remove samples - not going to remove anything at this time
  filter(year > 2013) %>% # taking out 2011 - 2013 - low sample sizes
  filter(month %in% c(5:7)) # taking out shoulder months with few data

length(which(dat$CPUE>0))
   
# # subset to MR or CA only 
# dat <- dat %>%
#   filter(Treatment == "CA")

# reset dir if not running full model
#dir <- file.path("C:/Users/daubleal/Desktop/2023 Assessment Cycle/Index_MARINE RESERVES/black_oregon_marres_hnl_CAonly")
#setwd(dir)

# define my covars
covars <- c("year", "region", "treatment","month") # full model, starting with no temp data
#covars <- c("year", "season" ,"Site","depth_bin")

#Ensure columns named appropriately and covariates are factors
dat <- dat %>%
  #filter(Treatment == "MR") %>% # testing one with just the MR sites
  mutate_at(covars, as.factor) 

# model selection 

model.full <- MASS::glm.nb(
  SFLA_count ~ year + region + treatment + month + offset(logEffort),
  data = dat,
  na.action = "na.fail") # I'm getting errors using "na.fail" here (I also tried changing it but then it wouldn't work in the dredge function below) 
summary(model.full)
anova(model.full)
#use ggpredict to get an estimate of the logEffort for sdmTMB predictions
#MuMIn will fit all models and then rank them by AICc
model.suite <- MuMIn::dredge(model.full,
                             rank = "AICc", 
                             fixed= c("offset(logEffort)", "year"))

#Create model selection dataframe for the document
Model_selection <- as.data.frame(model.suite) %>%
  dplyr::select(-weight)
Model_selection

# sdmTMB model 

#set the grid
grid <- expand.grid(
  year = unique(dat$year),
  region = levels(dat$region)[1],
  #treatment = levels(dat$treatment)[1]
  month = levels(dat$month)[1]
  #season = levels(dat$season)[1],
  #Site = levels(dat$Site)[1],
  #temp_bin = levels(dat$temp_bin)[1]
  #rock_bin = levels(dat$rock_bin)[1]
)

fit.nb <- sdmTMB(
  SFLA_count ~ year + region + month,
  data = dat,
  offset = dat$logEffort,
  time = "year",
  spatial="off",
  spatiotemporal = "off",
  family = nbinom2(link = "log"),
  control = sdmTMBcontrol(newton_loops = 1)) #documentation states sometimes aids convergence?

#}

#Get diagnostics and index for SS
do_diagnostics(
  dir = file.path(dir), 
  fit = fit.nb,
  plot_resids = F)

calc_index(
  dir = file.path(dir), 
  fit = fit.nb,
  grid = grid)

#-------------------------------------------------------------------------------
#Format data filtering table and the model selection table for document

# will need to modify my filter script to get the dataFilters to work
#View(dataFilters)

dataFilters <- dataFilters %>%
  rowwise() %>%
  filter(!all(is.na(across((everything()))))) %>%
  ungroup() %>%
  rename(`Positive Samples` = Positive_Samples)
dataFilters <- data.frame(lapply(dataFilters, as.character), stringsasFactors = FALSE)
#write.csv(dataFilters, file = file.path(dir, "data_filters.csv"), row.names = FALSE)

#View(Model_selection)
#format table for the document
out <- Model_selection %>%
  dplyr::select(-`(Intercept)`) %>%
  mutate_at(vars(covars,"year","offset(logEffort)"), as.character) %>%
  mutate(across(c("logLik","AICc","delta"), round, 1)) %>%
  # replace_na(list(district = "Excluded",                      # fix these later
  #                 targetSpecies = "Excluded", month = "Excluded")) %>% # fix later
  mutate_at(c(covars,"year","offset(logEffort)"), 
            funs(stringr::str_replace(.,"\\+","Included"))) %>%
  rename(`Effort offset` = `offset(logEffort)`, 
         `log-likelihood` = logLik) %>%
  rename_with(stringr::str_to_title,-AICc)
View(out)
#write.csv(out, file = file.path(dir,  "model_selection.csv"), row.names = FALSE)

#summary of trips and  percent pos per year
summaries <- dat %>%
  group_by(year) %>%
  summarise(tripsWithTarget = sum(SFLA_count>0),
            tripsWOTarget = sum(SFLA_count==0)) %>%
  mutate(totalTrips = tripsWithTarget+tripsWOTarget,
         percentpos = round(tripsWithTarget/(tripsWithTarget+tripsWOTarget),2)) 
View(summaries)
#write.csv(summaries, file.path(dir,  "percent_pos.csv"), row.names=FALSE)




######### testing out models with temp covariates ########

# so there are 31 NAs for temp - these will need to be removed and right now, i can only run factor variables?
# so attempting to bin some of Megan's new variables

# first reload data and start from scratch

# load data and add some stuff
dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS")
setwd(dir)
#load("data_for_GLM.RData")
dat<-read.csv("combined_settlement_updatedocean_v2.csv")

dat$region<-ifelse(dat$Site %in% c("HH","RR"),"South","North")
dat$treatment<-ifelse(dat$Site %in% c("HH","CF"),"Comparison","Reserve")
dat$Date<-as.POSIXct(dat$Date, format = "%Y-%m-%d")
dat$month<-as.numeric(format(dat$Date,"%m"))

# reset to a different directory 
dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/yellowtail_oregon_SMURF_addrolling_16_nomonth")
setwd(dir)

summary(dat$rolling16d)

# fix or add anything 
dat <- dat %>%
  #rename(Year = year) %>%
  rename(Effort = Sampling.Interval) %>%
  mutate(logEffort = log(Effort)) %>% 
  #create temp bins for the standardized index
  #mutate(temp_bin = cut(temp_index, breaks=c(-2,-1,0,1,2,3))) %>%
  #mutate(temp_bin = cut(temp_c_mid, breaks=c(6,7,8,9,10,11,15))) %>%
  #mutate(temp_bin = cut(cdd_8, breaks=c(50,60,70,80,90,120))) %>%
  #mutate(temp_bin = cut(cdd_16, breaks=c(110,125,140,155,170,215))) %>%
  #mutate(temp_bin = cut(rolling8d, breaks=c(-1.6,-1,-0.5,0,0.5,2))) %>%
  mutate(temp_bin = cut(rolling16d, breaks=c(-1.4,-1,-0.5,0,0.5,2.1))) %>%
  # remove samples
  filter(year > 2013) %>% # taking out 2011 - 2013 - low sample sizes
  filter(month %in% c(5:7)) %>% # taking out shoulder months with few data/positives
  filter(!is.na(temp_bin))

# define my covars
#covars <- c("year", "region", "treatment","month") # full model, starting with no temp data
covars <- c("year", "region","temp_bin")

#Ensure columns named appropriately and covariates are factors
dat <- dat %>%
  #filter(Treatment == "MR") %>% # testing one with just the MR sites
  mutate_at(covars, as.factor) 

# some boxplots for fun 

# ggplot(dat, aes(temp_bin, SFLA_settlement_rate)) +  
#   geom_boxplot() + 
#   labs(x = "rolling8d") 
# ggsave("rolling8d(binned).png",width = 10,height = 5)

# model selection 

model.full <- MASS::glm.nb(
  SFLA_count ~ year + region + temp_bin + offset(logEffort), 
  data = dat,
  na.action = "na.fail") 
summary(model.full)
anova(model.full)
#use ggpredict to get an estimate of the logEffort for sdmTMB predictions
#MuMIn will fit all models and then rank them by AICc
model.suite <- MuMIn::dredge(model.full,
                             rank = "AICc", 
                             fixed= c("offset(logEffort)", "year"))

#Create model selection dataframe for the document
Model_selection <- as.data.frame(model.suite) %>%
  dplyr::select(-weight)
Model_selection

# sdmTMB model 

#set the grid
grid <- expand.grid(
  year = unique(dat$year),
  region = levels(dat$region)[1],
  #treatment = levels(dat$treatment)[1]
  #month = levels(dat$month)[1],
  #season = levels(dat$season)[1],
  #Site = levels(dat$Site)[1],
  temp_bin = levels(dat$temp_bin)[1]
  #rock_bin = levels(dat$rock_bin)[1]
)

fit.nb <- sdmTMB(
  SFLA_count ~ year + region + temp_bin,
  data = dat,
  offset = dat$logEffort,
  time = "year",
  spatial="off",
  spatiotemporal = "off",
  family = nbinom2(link = "log"),
  control = sdmTMBcontrol(newton_loops = 1)) #documentation states sometimes aids convergence?

#}

#Get diagnostics and index for SS
do_diagnostics(
  dir = file.path(dir), 
  fit = fit.nb,
  plot_resids = F)

calc_index(
  dir = file.path(dir), 
  fit = fit.nb,
  grid = grid)

##### SUPPLEMENTARY TABLES #### 

#View(Model_selection)
#format table for the document
out <- Model_selection %>%
  dplyr::select(-`(Intercept)`) %>%
  mutate_at(vars(covars,"year","offset(logEffort)"), as.character) %>%
  mutate(across(c("logLik","AICc","delta"), round, 1)) %>%
  # replace_na(list(district = "Excluded",                      # fix these later
  #                 targetSpecies = "Excluded", month = "Excluded")) %>% # fix later
  mutate_at(c(covars,"year","offset(logEffort)"), 
            funs(stringr::str_replace(.,"\\+","Included"))) %>%
  rename(`Effort offset` = `offset(logEffort)`, 
         `log-likelihood` = logLik) %>%
  rename_with(stringr::str_to_title,-AICc)
View(out)
#write.csv(out, file = file.path(dir,  "model_selection.csv"), row.names = FALSE)

#summary of trips and  percent pos per year
summaries <- dat %>%
  group_by(year) %>%
  summarise(tripsWithTarget = sum(SFLA_count>0),
            tripsWOTarget = sum(SFLA_count==0)) %>%
  mutate(totalTrips = tripsWithTarget+tripsWOTarget,
         percentpos = round(tripsWithTarget/(tripsWithTarget+tripsWOTarget),2)) 
View(summaries)
#write.csv(summaries, file.path(dir,  "percent_pos.csv"), row.names=FALSE)




######## preliminary comparison with all temp models (all binned!) ###### 

# results notes - had to bin all of them!  didn't go for the mid-depths but could add that later 
# the full model for comparison is year + region + month (treatment was not significant)
# all temp_bins were added to the model
# the temp index was the strongest candidate but the model wouldn't converge when month was included? 
# however, the other two did work with month included 
# (some were and weren't significant! but I forgot to take notes on this)
# rolling8d was not significant,
# all models provided similar indices with a big peak in 2021 but some variation in the magnitude

# more detailed notes are in the google drive doc 

# reset wd
dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/temp bin comparison")
setwd(dir)

# pull in all three
full_notemp<-read.csv("index_forSS_full_notemp.csv")
temp_index<-read.csv("index_forSS_tempindex.csv")
temp_mid<-read.csv("index_forSS_mid.csv")
cdd8<-read.csv("index_forSS_cdd8.csv")
rolling<-read.csv("index_forSS_rolling.csv")
cdd8_month<-read.csv("index_forSS_cdd8_nomonth.csv")
rolling_month<-read.csv("index_forSS_rolling_nomonth.csv")

# standardize 
full_notemp$std_full<-full_notemp$obs/mean(full_notemp$obs)
temp_index$std_index<-temp_index$obs/mean(temp_index$obs)
temp_mid$std_index<-temp_mid$obs/mean(temp_mid$obs)
cdd8$std_cdd8<-cdd8$obs/mean(cdd8$obs)
rolling$std_cdd8<-rolling$obs/mean(rolling$obs)
cdd8_month$std_cdd8<-cdd8_month$obs/mean(cdd8_month$obs)
rolling_month$std_cdd8<-rolling_month$obs/mean(rolling_month$obs)

# combine and plot
df1 <- data.frame(Year = full_notemp$year, Variable = full_notemp$std_full)
df2 <- data.frame(Year = temp_index$year, Variable = temp_index$std_index)
df3 <- data.frame(Year = temp_mid$year, Variable = temp_mid$std_index)
df4 <- data.frame(Year = cdd8$year, Variable = cdd8$std_cdd8)
df5 <- data.frame(Year = rolling$year, Variable = rolling$std_cdd8)
df6 <- data.frame(Year = cdd8_month$year, Variable = cdd8_month$std_cdd8)
df7 <- data.frame(Year = rolling_month$year, Variable = rolling_month$std_cdd8)


df8 <- df1 %>%  mutate(Index = 'full (no temp)') %>%
  bind_rows(df2 %>%
              mutate(Index = 'temp_index (no month)')) %>%
  bind_rows(df3 %>%
              mutate(Index = 'temp_c_mid (no region/month)')) %>%
  bind_rows(df4 %>% 
              mutate(Index = 'cdd8')) %>%
  bind_rows(df5 %>%
              mutate(Index = 'rolling8')) %>%
  bind_rows(df6 %>% 
              mutate(Index = 'cdd8 (no month)')) %>%
  bind_rows(df7 %>%
              mutate(Index = 'rolling8 (no month)'))

df8$Year<-as.factor(df8$Year)

ggplot(df8,aes(y = Variable,x = Year,color = Index)) + 
  geom_line(aes(group = Index))+
  geom_point(size = 1.5) +
  labs(y = "Standardized Index Value")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x = element_blank())
#ggsave("temperature (binned) index comparison.png",width = 10,height = 5)






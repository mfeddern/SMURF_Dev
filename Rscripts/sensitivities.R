library(here)
library(r4ss)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(NatParksPalettes)
library(ggpubr)
# option 1: install via {pak}
#install.packages("pak")
#pak::pkg_install("r4ss/r4ss")

colpal<-natparks.pals("DeathValley",n=6)
colpalLONG<-natparks.pals("DeathValley",n=12)
model_directory <- 'model_runs'
base_model_name <- '5.09_no_extra_se' #base model is the assessment base model == contains SMURF!
exe_loc <- here('model_runs/ss3.exe')
base_model <- SS_read(file.path(model_directory, base_model_name))
base_out <- SS_output(
  file.path(model_directory, base_model_name),
  # SpawnOutputLabel = "Spawning output (trillions of eggs)",
  printstats = FALSE,
  verbose = FALSE,
  covar =FALSE
)


# Write sensitivities -----------------------------------------------------


## No SMURF ------------------------------------------------------------------

sensi_mod <- base_model

sensi_mod$dat$CPUE <- sensi_mod$dat$CPUE |>
  filter(index != 7)

sensi_mod$ctl$Q_options <- sensi_mod$ctl$Q_options[-grep('SMURF', rownames(sensi_mod$ctl$Q_options)),]
sensi_mod$ctl$Q_parms <- sensi_mod$ctl$Q_parms[-grep('SMURF', rownames(sensi_mod$ctl$Q_parms)),]

SS_write(sensi_mod, file.path(model_directory, 'sensitivities', 'no_SMURF'), overwrite = TRUE)

no_SMURF <- SSgetoutput(dirvec = c(file.path(model_directory, base_model_name),
                              file.path(model_directory, 'sensitivities', 'oceanographic_index')))

## Setting base model to no SMURF 
no_smurf<- SS_read(file.path(model_directory, 'sensitivities', 'no_SMURF'))

 
## Oceanographic index2 -----------------------------------------------------
sensi_mod <-no_smurf

flt <- sensi_mod$dat$Nfleets + 1

ocean <- read.csv('Data/raw_not_confidential/OceanographicIndex/OceanographicIndexV1.csv') |>
  mutate(month = 7, 
         index = ifelse(year >= 2020, flt, -flt) # include 5 years of index
  ) |>
  select(year, month, index, obs = fit, se_log = se.p)

# data file updates
sensi_mod$dat$Nfleets <- flt
sensi_mod$dat$fleetnames[flt] <- 'ocean'
sensi_mod$dat$fleetinfo[flt,] <- c(3,1,1,2,0,'ocean')
sensi_mod$dat$CPUEinfo[flt,] <- c(flt,36,-1,0)
sensi_mod$dat$len_info[flt,] <- sensi_mod$dat$len_info[flt-1,]
sensi_mod$dat$age_info[flt,] <- sensi_mod$dat$age_info[flt-1,]
sensi_mod$dat$fleetinfo1$ocean <- sensi_mod$dat$fleetinfo1$WCGBTS
sensi_mod$dat$fleetinfo2$ocean <- sensi_mod$dat$fleetinfo2$WCGBTS
sensi_mod$dat$CPUE <- bind_rows(sensi_mod$dat$CPUE,
                                ocean)

# control file updates
sensi_mod$ctl$size_selex_types[flt,] <- rep(0, 4)
sensi_mod$ctl$age_selex_types[flt,] <- sensi_mod$ctl$age_selex_types[flt-1,]
sensi_mod$ctl$Q_options <- rbind(sensi_mod$ctl$Q_options,
                                 ocean = c(flt,1,0,0,0,0))
sensi_mod$ctl$Q_parms <- bind_rows(sensi_mod$ctl$Q_parms,
                                   sensi_mod$ctl$Q_parms[1,])

sensi_mod$ctl$Q_parms[nrow(sensi_mod$ctl$Q_parms), c('INIT', 'PHASE')] <- c(1, -99)

SS_write(sensi_mod, file.path(model_directory, 'sensitivities', 'oceanographic_index'), overwrite = TRUE)

## RREAS Coastwide-------------------------------------------------------------------

sensi_mod <- no_smurf

flt <- base_model$dat$Nfleets + 1

rreas <- read.csv('Data/raw_not_confidential/RREAS/ytail_coastwide_indices.csv') |>
  mutate(index = flt, month = 7) |>
  rename(se_log = logse, year = YEAR, obs = est)

sensi_mod$dat$Nfleets <- flt
sensi_mod$dat$fleetnames[flt] <- 'RREAS'
sensi_mod$dat$fleetinfo[flt,] <- c(3,1,1,2,0,'RREAS')
sensi_mod$dat$CPUEinfo[flt,] <- c(flt,33,0,0)
sensi_mod$dat$len_info[flt,] <- sensi_mod$dat$len_info[flt-1,]
sensi_mod$dat$age_info[flt,] <- sensi_mod$dat$age_info[flt-1,]
sensi_mod$dat$fleetinfo1$SMURF <- sensi_mod$dat$fleetinfo1$WCGBTS
sensi_mod$dat$fleetinfo2$SMURF <- sensi_mod$dat$fleetinfo2$WCGBTS
sensi_mod$dat$CPUE <- bind_rows(sensi_mod$dat$CPUE,
                                rreas)

# control file updates
sensi_mod$ctl$size_selex_types[flt,] <- rep(0, 4)
sensi_mod$ctl$age_selex_types[flt,] <- sensi_mod$ctl$age_selex_types[flt-1,]
sensi_mod$ctl$Q_options <- rbind(sensi_mod$ctl$Q_options,
                                 RREAS = c(flt,1,0,0,0,0))
sensi_mod$ctl$Q_parms <- bind_rows(sensi_mod$ctl$Q_parms,
                                   sensi_mod$ctl$Q_parms[1,])

SS_write(sensi_mod, file.path(model_directory, 'sensitivities', 'RREAS'),
         overwrite = TRUE)

## RREAS North-------------------------------------------------------------------

sensi_mod <- no_smurf

flt <- base_model$dat$Nfleets + 1

rreasn <- read.csv('Data/raw_not_confidential/RREAS/YOYNorth.csv') |>
  mutate(index = flt, month = 7) |>
  rename(se_log = logse, year = YEAR, obs = est)

sensi_mod$dat$Nfleets <- flt
sensi_mod$dat$fleetnames[flt] <- 'RREASN'
sensi_mod$dat$fleetinfo[flt,] <- c(3,1,1,2,0,'RREASN')
sensi_mod$dat$CPUEinfo[flt,] <- c(flt,33,0,0)
sensi_mod$dat$len_info[flt,] <- sensi_mod$dat$len_info[flt-1,]
sensi_mod$dat$age_info[flt,] <- sensi_mod$dat$age_info[flt-1,]
sensi_mod$dat$fleetinfo1$RREASN <- sensi_mod$dat$fleetinfo1$WCGBTS
sensi_mod$dat$fleetinfo2$RREASN <- sensi_mod$dat$fleetinfo2$WCGBTS
sensi_mod$dat$CPUE <- bind_rows(sensi_mod$dat$CPUE,
                                rreasn)

# control file updates
sensi_mod$ctl$size_selex_types[flt,] <- rep(0, 4)
sensi_mod$ctl$age_selex_types[flt,] <- sensi_mod$ctl$age_selex_types[flt-1,]
sensi_mod$ctl$Q_options <- rbind(sensi_mod$ctl$Q_options,
                                 RREASN = c(flt,1,0,0,0,0))
sensi_mod$ctl$Q_parms <- bind_rows(sensi_mod$ctl$Q_parms,
                                   sensi_mod$ctl$Q_parms[1,])

SS_write(sensi_mod, file.path(model_directory, 'sensitivities', 'RREASN'),
         overwrite = TRUE)

## OCNMS-------------------------------------------------------------------

sensi_mod <- no_smurf

flt <- base_model$dat$Nfleets + 1

OCNMS <- read.csv('Data/raw_not_confidential/Estimated-YOY-trend-coast.csv')|>
  rename(Index=grand.mean)|>
  select(year, Index, SE)|>
  mutate(index = flt, month = 7) |>
  rename(se_log = SE,  obs = Index)
  
sensi_mod$dat$Nfleets <- flt
sensi_mod$dat$fleetnames[flt] <- 'OCNMS'
sensi_mod$dat$fleetinfo[flt,] <- c(3,1,1,2,0,'OCNMS')
sensi_mod$dat$CPUEinfo[flt,] <- c(flt,33,0,0)
sensi_mod$dat$len_info[flt,] <- sensi_mod$dat$len_info[flt-1,]
sensi_mod$dat$age_info[flt,] <- sensi_mod$dat$age_info[flt-1,]
sensi_mod$dat$fleetinfo1$SMURF <- sensi_mod$dat$fleetinfo1$WCGBTS
sensi_mod$dat$fleetinfo2$SMURF <- sensi_mod$dat$fleetinfo2$WCGBTS
sensi_mod$dat$CPUE <- bind_rows(sensi_mod$dat$CPUE,
                                OCNMS)

# control file updates
sensi_mod$ctl$size_selex_types[flt,] <- rep(0, 4)
sensi_mod$ctl$age_selex_types[flt,] <- sensi_mod$ctl$age_selex_types[flt-1,]
sensi_mod$ctl$Q_options <- rbind(sensi_mod$ctl$Q_options,
                                 OCNMS = c(flt,1,0,0,0,0))
sensi_mod$ctl$Q_parms <- bind_rows(sensi_mod$ctl$Q_parms,
                                   sensi_mod$ctl$Q_parms[1,])

SS_write(sensi_mod, file.path(model_directory, 'sensitivities', 'OCNMS'),
         overwrite = TRUE)
# Run stuff ---------------------------------------------------------------

sensi_dirs <- list.files(file.path(model_directory, 'sensitivities'))

tuning_mods <- grep('weighting', sensi_mod)

future::plan(future::multisession(workers = parallelly::availableCores(omit = 1)))


furrr::future_map(c("OCNMS"), \(x) 
                  run(file.path(model_directory, 'sensitivities', x), 
                      exe = exe_loc, extras = '-nohess', skipfinished = FALSE)
)

future::plan(future::sequential)
# Plot stuff --------------------------------------------------------------
out <- SSgetoutput(dirvec = c(file.path(model_directory, base_model_name),
                              file.path(model_directory, 'sensitivities', 'no_SMURF'),
                              file.path(model_directory, 'sensitivities', 'OCNMS'),
                              file.path(model_directory, 'sensitivities', 'RREAS'),
                              file.path(model_directory, 'sensitivities', 'RREASN'),
                              file.path(model_directory, 'sensitivities', 'oceanographic_index')))
SSsummarize(out) |>
  SSplotComparisons(subplots = c(1,3), endyrvec = 2037, new = FALSE)

## function ----------------------------------------------------------------

make_detailed_sensitivites <- function(biglist, mods, 
                                       outdir, grp_name) {
  
  shortlist <- big_sensitivity_output[c('base', mods$dir)] |>
    r4ss::SSsummarize() 
  
  r4ss::SSplotComparisons(shortlist,
                          subplots = c(2,4, 18), 
                          print = TRUE,  
                          plot = FALSE,
                          plotdir = outdir, 
                          filenameprefix = grp_name,
                          legendlabels = c('Base', mods$pretty), 
                          endyrvec = 2036)

  r4ss::plot_twopanel_comparison(big_sensitivity_output[c('base', mods$dir)],
                          dir = outdir, 
                          filename = paste0(grp_name, '_comparison.png'),
                          legendlabels = c('Base', mods$pretty), 
                          legendloc = 'bottomleft',
                          endyrvec = 2036)
  
  SStableComparisons(shortlist, 
                     modelnames = c('Base', mods$pretty),
                     names =c("Recr_Virgin", "R0", "NatM", "L_at_Amax", "VonBert_K", "SmryBio_unfished", "SSB_Virg",
                              "SSB_2025", "Bratio_2025", "SPRratio_2024", "LnQ_base_WCGBTS"),
                     likenames = c(
                              "TOTAL", "Survey", "Length_comp", "Age_comp",
                              "Discard", "Mean_body_wt", "Recruitment", "priors"
                            )
                          ) |>
    # dplyr::filter(!(Label %in% c('NatM_break_1_Fem_GP_1',
    #                              'NatM_break_1_Mal_GP_1', 'NatM_break_2_Mal_GP_1')),
    #               Label != 'NatM_uniform_Mal_GP_1' | any(grep('break', Label)),
    #               Label != 'SR_BH_steep' | any(grep('break', Label))) |>
    # dplyr::mutate(dplyr::across(-Label, ~ sapply(., format, digits = 3, scientific = FALSE) |>
    #                               stringr::str_replace('NA', ''))) |>
    `names<-`(c('Label', 'Base', mods$pretty)) |>
    write.csv(file.path(outdir, paste0(grp_name, '_table.csv')), 
              row.names = FALSE, )
  
}



## grouped plots -----------------------------------------------------------


indices <- data.frame(dir = c('no_smurf',
                              'OCNMS',
                              "RREASN",
                              'oceanographic_index',
                              'RREAS'),
                      pretty = c('- No Indices',
                                 '+ OCNMS',
                                 '+ RREAS North',
                                 '+ Oceanographic',
                                 '+ RREAS Coastwide')
)


sens_names <- bind_rows(indices)

big_sensitivity_output <- SSgetoutput(
  dirvec = file.path(
    model_directory,
    c(
      base_model_name,
      glue::glue("sensitivities/{subdir}", subdir = sens_names$dir)
    )
  ),
 # SpawnOutputLabel = "Spawning output (trillions of eggs)"
) |>
  `names<-`(c('base', sens_names$dir))

saveRDS(object = big_sensitivity_output, file = "data/processed/index_sensitivities_forplots.rds")

# test to make sure they all read correctly:
which(sapply(big_sensitivity_output, length) < 180) # all lengths should be >180

sens_names_ls <- list(
                      indices = indices)

outdir <- 'figures/sensitivities'

purrr::imap(sens_names_ls, \(sens_df, grp_name) 
            make_detailed_sensitivites(biglist = big_sensitivity_output, 
                                       mods = sens_df, 
                                       outdir = outdir, 
                                       grp_name = grp_name))


## spawning output -------------------------------------------------------------

sensitivity_output <- SSsummarize(big_sensitivity_output) 

SpawningBiomass <- sensitivity_output$SpawnBio
colnames(SpawningBiomass)<- c( "NoIndices","SMURF", "OCNMS", "RREASN", "Oceanographic", "RREAS", "Label", "Year")
SpawningData<-SpawningBiomass%>%pivot_longer(c( NoIndices,SMURF, OCNMS, RREASN, Oceanographic, RREAS))%>%
  rename(SpawnBio=value)%>%
  mutate(Index=ifelse(name=="NoIndices", "No Indices", 
                      ifelse(name=="SMURF", "SMURF",
                             ifelse(name=="OCNMS", "OCNMS",
                                    ifelse(name=="RREASN","RREAS North",
                                           ifelse(name=="Oceanographic", "Oceanographic","RREAS Coastwide"))))))
  

spawnlong<-ggplot(SpawningData) +
  geom_rect(xmin=2025, xmax=2050,ymin=2.25, ymax=16, alpha=0.1,fill='#F0F0F0',col='#F0F0F0')+
  geom_point(aes(x = Year, y = SpawnBio, col = Index, pch = Index))+
  geom_line(aes(x = Year, y = SpawnBio, col = Index, pch = Index))+
  ylab("Spawning Output")+
  ylim(c(3, 15))+
  xlim(c(1900, 2036))+
  scale_color_manual(values=colpal)+
  geom_text(x = 2032, y=14, label="Forecast \n Period")+
  scale_x_continuous(limits = c(1900, 2036), oob = oob_keep)+
  theme_bw()

spawnshort<-ggplot(SpawningData%>%filter(Year>=2000)) +
  geom_rect(xmin=2025, xmax=2050,ymin=2.25, ymax=16, alpha=0.1,fill='#F0F0F0',col='#F0F0F0')+
  geom_point(aes(x = Year, y = SpawnBio, col = Index, pch = Index))+
  geom_line(aes(x = Year, y = SpawnBio, col = Index))+
  ylab("Spawning Output")+
  ylim(c(3, 15))+
  scale_color_manual(values=colpal)+
  scale_x_continuous(limits = c(1900, 2036), oob = oob_keep)+
  geom_text(x = 2032, y=14, label="Forecast \n Period")+
  xlim(c(2000, 2036))+
  theme_bw()

## fraction unfished -------------------------------------------------------------

unfished <- sensitivity_output$Bratio
colnames(unfished)<- c("SMURF", "NoIndices", "OCNMS", "RREASN", "Oceanographic", "RREAS", "Label", "Year")
FractionData<-unfished%>%pivot_longer(c(SMURF, NoIndices, OCNMS, RREASN, Oceanographic, RREAS))%>%
  rename(Bratio=value)%>%
  mutate(Index=ifelse(name=="SMURF", "SMURF", 
                      ifelse(name=="NoIndices", "No Indices",
                             ifelse(name=="OCNMS", "OCNMS",
                                    ifelse(name=="RREASN","RREAS North",
                                           ifelse(name=="Oceanographic", "Oceanographic","RREAS Coastwide"))))))


fractionshort<-ggplot(FractionData) +
  geom_rect(xmin=2025, xmax=2045,ymin=-2.25, ymax=16, alpha=0.1,fill='#F0F0F0',col='#F0F0F0')+
  geom_point(aes(x = Year, y = Bratio, col = Index, pch = Index))+
  geom_line(aes(x = Year, y = Bratio, col = Index, pch = Index))+
  ylab("Fraction Unfished")+
  ylim(c(0, 1.05))+
  geom_text(x = 2032, y=1, label="Forecast \n Period")+
  scale_color_manual(values=colpal)+
  geom_hline(yintercept=0.40, lty=2, col =colpalLONG[1])+
  geom_hline(yintercept=0.25, lty=2, col =colpalLONG[1])+
  geom_text(x = 2006, y=0.45, label="Management target")+
  geom_text(x = 2006, y=0.3, label="Minimum stock size threshold")+
  geom_hline(yintercept=1, lty=2)+
  scale_x_continuous(limits = c(2000, 2036), oob = oob_keep)+
  geom_text(x = 2032, y=14, label="Forecast Period")+
  theme_bw()

fractionlong<-ggplot(FractionData) +
  geom_rect(xmin=2025, xmax=2045,ymin=-2.25, ymax=16, alpha=0.1,fill='#F0F0F0',col='#F0F0F0')+
  geom_point(aes(x = Year, y = Bratio, col = Index, pch = Index))+
  geom_line(aes(x = Year, y = Bratio, col = Index, pch = Index))+
  ylab("Fraction Unfished")+
  ylim(c(0, 1.05))+
  geom_text(x = 2032, y=1, label="Forecast \n Period")+
  scale_color_manual(values=colpal)+
  geom_hline(yintercept=0.40, lty=2, col =colpalLONG[1])+
  geom_hline(yintercept=0.25, lty=2, col =colpalLONG[1])+
  geom_text(x = 1910, y=0.45, label="Management target")+
  geom_text(x = 1915, y=0.3, label="Minimum stock size threshold")+
  geom_hline(yintercept=1, lty=2)+
  scale_x_continuous(limits = c(1900, 2036), oob = oob_keep)+
  geom_text(x = 2032, y=14, label="Forecast Period")+
  theme_bw()

## full plot ----------------------------------------------------------------

#all_ts<- ggarrange(ps_ts, sb_ts,yt_ts, hk_ts,ncol = 1, nrow = 4)
#all_ts
pdf(file = "Figures/Manuscript/MainText/pdf/spawnplot.pdf", width =8, height =6)
spawnshort
dev.off()

ggsave("Figures/SMURF_Manuscript/png/spawnplot.png",  dpi = 300,  
       width = 7, height = 5, units = "in")
spawnshort
dev.off() 
## big plot ----------------------------------------------------------------

current.year <- 2025
CI <- 0.95



lapply(big_sensitivity_output, function(.)
  .$warnings[grep('gradient', .$warnings)]) # check gradients

dev.quants.SD <- c(
  sensitivity_output$quantsSD[sensitivity_output$quantsSD$Label == "SSB_Initial", 1],
  (sensitivity_output$quantsSD[sensitivity_output$quantsSD$Label == paste0("SSB_", current.year), 1]),
  sensitivity_output$quantsSD[sensitivity_output$quantsSD$Label == paste0("Bratio_", current.year), 1],
  sensitivity_output$quantsSD[sensitivity_output$quantsSD$Label == "Dead_Catch_SPR", 1],
  sensitivity_output$quantsSD[sensitivity_output$quantsSD$Label == "annF_SPR", 1]
)

dev.quants <- rbind(
  sensitivity_output$quants[sensitivity_output$quants$Label == "SSB_Initial", 
                            1:(dim(sensitivity_output$quants)[2] - 2)],
  sensitivity_output$quants[sensitivity_output$quants$Label == paste0("SSB_", current.year), 
                            1:(dim(sensitivity_output$quants)[2] - 2)],
  sensitivity_output$quants[sensitivity_output$quants$Label == paste0("Bratio_", current.year), 
                            1:(dim(sensitivity_output$quants)[2] - 2)],
  sensitivity_output$quants[sensitivity_output$quants$Label == "Dead_Catch_SPR", 
                            1:(dim(sensitivity_output$quants)[2] - 2)],
  sensitivity_output$quants[sensitivity_output$quants$Label == "annF_SPR", 
                            1:(dim(sensitivity_output$quants)[2] - 2)]
) |>
  cbind(baseSD = dev.quants.SD) |>
  dplyr::mutate(Metric = c("SB0", paste0("SSB_", current.year), paste0("Bratio_", current.year), "MSY_SPR", "F_SPR")) |>
  tidyr::pivot_longer(-c(base, Metric, baseSD), names_to = 'Model', values_to = 'Est') |>
  dplyr::mutate(relErr = (Est - base)/base,
                logRelErr = log(Est/base),
                mod_num = rep(1:nrow(sens_names), 5))

metric.labs <- c(
  SB0 = expression(SB[0]),
  SSB_2023 = as.expression(bquote("SB"[.(current.year)])),
  Bratio_2023 = bquote(frac(SB[.(current.year)], SB[0])),
  MSY_SPR = expression(Yield['SPR=0.50']),
  F_SPR = expression(F['SPR=0.50'])
)

CI.quants <- dev.quants |>
  dplyr::filter(Model == unique(dev.quants$Model)[1]) |>
  dplyr::select(base, baseSD, Metric) |>
  dplyr::mutate(CI = qnorm((1-CI)/2, 0, baseSD)/base)

ggplot(dev.quants, aes(x = relErr, y = mod_num, col = Metric, pch = Metric)) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_point() +
  geom_segment(aes(x = CI, xend = abs(CI), col = Metric,
                   y = nrow(sens_names) + 1.5 + seq(-0.5, 0.5, length.out = length(metric.labs)),
                   yend = nrow(sens_names) + 1.5 + seq(-0.5, 0.5, length.out = length(metric.labs))), 
               data = CI.quants, linewidth = 2, show.legend = FALSE, lineend = 'round') +
  theme_bw() +
  scale_shape_manual(
    values = c(15:18, 12),
    # name = "",
    labels = metric.labs
  ) +
  # scale_color_discrete(labels = metric.labs) +
  scale_y_continuous(breaks = 1:nrow(sens_names), name = '', labels = sens_names$pretty, 
                     limits = c(1, nrow(sens_names) + 2), minor_breaks = NULL) +
  xlab("Relative change") +
  viridis::scale_color_viridis(discrete = TRUE, labels = metric.labs)
ggsave(file.path(outdir, 'sens_summary.png'),  dpi = 300,  
       width = 6, height = 6.5, units = "in")


## comparison plot ----------------------------------------------------------------

SMURF <- read.csv('Data/raw_not_confidential/SMURF index/index_forSS.csv') |>
  select(year, month, obs, se_log = logse)|>
  mutate(index=7)
no_index<- read.csv('Data/raw_not_confidential/SMURF index/ghost.csv')|>
  mutate(index=7)

dat_long<-bind_rows(OCNMS%>%mutate(Index="OCNMS"),rreasn%>%mutate(Index="RREAS North"))%>%
  bind_rows(ocean%>%mutate(Index="Oceanographic"))%>%
  bind_rows(SMURF%>%mutate(Index="SMURF"))%>%
  bind_rows(rreas%>%mutate(Index="RREAS Coastwide"))%>%
  bind_rows(no_index%>%mutate(Index="No Index"))

dat_stand<-dat_long%>%group_by(Index)%>%mutate(Stand=scale(obs))%>%ungroup()
dat_stand2<-dat_long%>%filter(year>2015)%>%group_by(Index)%>%mutate(Stand=scale(obs))%>%ungroup()%>%
  bind_rows(no_index%>%mutate(Index="No Index"))

stand<-ggplot(dat_stand%>%filter(index>0)) +
  geom_point(aes(x = year, y = Stand, col = Index, pch = Index))+
  geom_line(aes(x = year, y = Stand, col = Index, pch = Index))+
  ylab("Standardized Index")+
 xlim(c(2000, 2025))+
  xlab("")+
 # geom_text(x = 2032, y=1, label="Forecast \n Period")+
  scale_color_manual(values=colpal)+
#  geom_hline(yintercept=0.40, lty=2, col =colpalLONG[1])+
 # geom_hline(yintercept=0.25, lty=2, col =colpalLONG[1])+
 # geom_text(x = 1910, y=0.45, label="Management target")+
#  geom_text(x = 1915, y=0.3, label="Minimum stock size threshold")+
#  geom_hline(yintercept=1, lty=2)+
#  scale_x_continuous(limits = c(1900, 2036), oob = oob_keep)+
 # geom_text(x = 2032, y=14, label="Forecast Period")+
  theme_bw()+
  theme(legend.position = "none")

stand2016<-ggplot(dat_stand2%>%filter(index>0)) +
  geom_point(aes(x = year, y = Stand, col = Index, pch = Index))+
  geom_line(aes(x = year, y = Stand, col = Index, pch = Index))+
  ylab("")+
  xlab("")+
  xlim(c(2016, 2024))+
  # ylim(c(0, 1.05))+
  # geom_text(x = 2032, y=1, label="Forecast \n Period")+
  scale_color_manual(values=colpal)+
  #  geom_hline(yintercept=0.40, lty=2, col =colpalLONG[1])+
  # geom_hline(yintercept=0.25, lty=2, col =colpalLONG[1])+
  # geom_text(x = 1910, y=0.45, label="Management target")+
  #  geom_text(x = 1915, y=0.3, label="Minimum stock size threshold")+
  #  geom_hline(yintercept=1, lty=2)+
  #  scale_x_continuous(limits = c(1900, 2036), oob = oob_keep)+
  # geom_text(x = 2032, y=14, label="Forecast Period")+
  theme_bw()

indices<-annotate_figure(ggarrange(stand, stand2016,ncol = 2, nrow = 1, widths = c(1,1.5),labels=c("A.", "B.")),
                bottom = "Year")

pdf(file = "Figures/Manuscript/yoyindices.pdf", width =9, height =6)
indices
dev.off()

ggsave("Figures/Manuscript/yoyindices.png",  dpi = 300,  
       width = 9, height = 5, units = "in", bg="white")
indices
dev.off() 


#### csvs ####
likelihoods<- sensitivity_output[["likelihoods"]]
write.csv(likelihoods, "data/processed/likelihoodsSens.csv")

recdevs<-sensitivity_output[["recdevs"]]
write.csv(recdevs, "data/processed/RecruitmentDeviationsSens.csv")

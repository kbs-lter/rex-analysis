# TITLE:          REX: 2022 Solidago leaves analysis
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_VOC L0 folder
# DATA OUTPUT:    Calculating total biomass of leaves measured for VOCs in each REX treatment
# PROJECT:        REX
# DATE:           July 2022


### Data importing & checking ###
# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Read in data
length_weight <- read.csv(file.path(dir, "T7_warmx_VOC/L0/REX_T7_VOC_2022_Solidago_length_weight_L0.csv"))
voc_leaves <- read.csv(file.path(dir, "T7_warmx_VOC/L0/REX_T7_VOC_2022_Solidago_leaves_L0.csv"))

# checking spreadsheets
unique(length_weight$Treatment)
unique(voc_leaves$Treatment)

# removing plants from voc_leaves dataframe
# removing these because I lost their VOC data, so I need to make sure I am not calculating more biomass than what was measured for VOCs
voc_leaves <- voc_leaves %>%
  filter(!(Rep == 1 & Treatment == "Warmed" & Plant_Number == 1)) %>%
  filter(!(Rep == 2 & Treatment == "Ambient" & Plant_Number == 1)) %>%
  filter(!(Rep == 4 & Treatment == "Warmed" & Plant_Number == 1)) %>%
  filter(!(Rep == 4 & Treatment == "Warmed_Drought" & Plant_Number == 1))


### Getting regression btwn leaf length and weight for each treatment ###
# making dataframes for each treatment in regression dataframe
warmed <- length_weight %>%
  filter(Treatment == "Warmed")
ambient <- length_weight %>%
  filter(Treatment == "Ambient")
irr <- length_weight %>%
  filter(Treatment == "Irrigated")
drought <- length_weight %>%
  filter(Treatment == "Drought")
warm_drought <- length_weight %>%
  filter(Treatment == "Warmed_Drought")

# regression btwn length and weight for leaves in each treatment
w_mod <- lm(Weight_g ~ Length_cm, data = warmed)
summary(w_mod)$coef
plot(Weight_g ~ Length_cm, data=warmed)
abline(w_mod)
a_mod <- lm(Weight_g ~ Length_cm, data = ambient)
summary(a_mod)$coef
plot(Weight_g ~ Length_cm, data=ambient)
abline(a_mod)
i_mod <- lm(Weight_g ~ Length_cm, data = irr)
summary(i_mod)$coef
plot(Weight_g ~ Length_cm, data=irr)
abline(i_mod)
d_mod <- lm(Weight_g ~ Length_cm, data = drought)
summary(d_mod)$coef
plot(Weight_g ~ Length_cm, data=drought)
abline(d_mod)
wd_mod <- lm(Weight_g ~ Length_cm, data = warm_drought)
summary(wd_mod)$coef
plot(Weight_g ~ Length_cm, data=warm_drought)
abline(wd_mod)

# making dataframe for each treatment in voc_leaves dataframe
warmed_voc <- voc_leaves %>%
  filter(Treatment == "Warmed") %>%
  select(-Treatment, -Rep, -Plant_Number)
ambient_voc <- voc_leaves %>%
  filter(Treatment == "Ambient") %>%
  select(-Treatment, -Rep, -Plant_Number)
irr_voc <- voc_leaves %>%
  filter(Treatment == "Irrigated") %>%
  select(-Treatment, -Rep, -Plant_Number)
drought_voc <- voc_leaves %>%
  filter(Treatment == "Drought") %>%
  select(-Treatment, -Rep, -Plant_Number)
warm_drought_voc <- voc_leaves %>%
  filter(Treatment == "Warmed_Drought") %>%
  select(-Treatment, -Rep, -Plant_Number)

# predicting biomass from leaf lengths
sum(predict(w_mod, warmed_voc)) # 65.61699 g sampled
sum(predict(a_mod, ambient_voc)) # 60.93225 g sampled
sum(predict(i_mod, irr_voc)) # 65.71546 g sampled
sum(predict(d_mod, drought_voc)) # 57.64733 g sampled
sum(predict(wd_mod, warm_drought_voc)) # 51.75015 g sampled



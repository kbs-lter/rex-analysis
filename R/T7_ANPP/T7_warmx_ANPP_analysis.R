# TITLE:          REX: Zarnetske plots ANPP analyses
# AUTHORS:        Kara Dobson, Moriah Young
# COLLABORATORS:  Phoebe Zarnetske, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_ANPP L1 folder
# DATA OUTPUT:    Analyses
# PROJECT:        REX
# DATE:           April 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(lmerTest)
library(car)
library(bbmle)
library(sjPlot)
library(emmeans)
library(stats)
library(plotrix)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Read in data
anpp <- read.csv(file.path(dir, "T7_ANPP/L1/T7_warmx_ANPP_L1.csv")) # 2019-2023 ANPP data

# filter out 2019 data and non 0.2 data
anpp1 <- anpp %>%
        filter(!(Year=='2019')) %>%
        filter(Scale_meter_square==0.2)

# calculate plot level data???


###### Data exploration #######
# The first thing we do is check the distribution of the raw data
# We're checking for a normal distribution; the distribution of the raw data
# doesn't matter as much as the distribution of a model's residuals (which we'll also check),
# but its a good way to get a glimpse of how the data look
descdist(anpp1$plant_biomass_gm2, discrete = FALSE) # checks what type of distribution your data matches
hist(anpp1$plant_biomass_gm2)
qqnorm(anpp1$plant_biomass_gm2) # should be a straight diagonal line if normal
shapiro.test(anpp1$plant_biomass_gm2) # p > 0.05 if normal
# making a model to test distribution of model residuals
raw.data.test <- lm(plant_biomass_gm2 ~ Subplot_Descriptions * Year, data=anpp1) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))

# log transformation
anpp1$log_biomass <- log(anpp1$plant_biomass_gm2)
descdist(anpp1$log_biomass, discrete = FALSE)
hist(anpp1$log_biomass)
qqnorm(anpp1$log_biomass) # a little weird
shapiro.test(anpp1$log_biomass) # still not normal
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(log_biomass ~ Subplot_Descriptions * Year, data=anpp1) # testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # not normal
# better - going with log transformation

# beta transformation


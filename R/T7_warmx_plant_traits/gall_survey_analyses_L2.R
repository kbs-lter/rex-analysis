# TITLE:          REX: Gall density analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Results from analyses
# PROJECT:        REX
# DATE:           July 2021

# Clear all existing data
rm(list=ls())

# Load packages
library(bbmle)
library(lmerTest)
library(fitdistrplus)
library(sjPlot)
library(tidyverse)
library(car)
library(emmeans)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
gall <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_gall_survey_L1.csv"))


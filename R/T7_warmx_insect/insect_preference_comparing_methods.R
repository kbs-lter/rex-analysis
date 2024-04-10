# TITLE:          REX: comparing 2021 & 2022 insect preference trial methods
# AUTHORS:        Moriah Young
# COLLABORATORS:  
# DATA INPUT:     
# DATA OUTPUT:    Analyses
# PROJECT:        REX
# DATE:           April 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(plotrix)
library(lmerTest)
library(fitdistrplus)
library(emmeans)

#set directory
dir<-Sys.getenv("DATA_DIR")

#read in data
insects21 <- read.csv(file.path(dir, "T7_warmx_insect/L1/T7_warmx_insect_preference_2021_L1.csv"))
insects22 <- read.csv(file.path(dir, "T7_warmx_insect/L0/T7_warmx_insect_preference_2022_L0.csv"))
meta <- read.csv(file.path(dir, "REX_warmx_metadata.csv"))





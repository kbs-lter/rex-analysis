# TITLE:          REX: Abiotic analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Mark Hammond, Moriah Young, Emily Parker
# DATA INPUT:     Data imported as csv files from shared REX Google drive sensor data L1 folder
# DATA OUTPUT:    Analyses of abiotic data
# PROJECT:        REX
# DATE:           July 2023


# Clear all existing data
rm(list=ls())


# Set working directory
dir <- Sys.getenv("DATA_DIR")
list.files(dir)


# Load packages
library(tidyverse)
library(fitdistrplus)
library(lmerTest)
library(glmmTMB)
library(scales)
library(emmeans)


# Read in data
hobo_data <- read.csv(file.path(dir, "sensors/OTC Footprints/L1/T7_warmx_HOBO_L1.csv"), stringsAsFactors=T)
soil_data <- read.csv(file.path(dir, "sensors/OTC Footprints/L1/T7_warmx_soil_sensors_L1.csv"))
str(hobo_data)
str(soil_data)
# re-making the data column a date
hobo_data$Date_Time <- as.POSIXct(hobo_data$Date_Time, format = "%Y-%m-%d %H:%M:%S")
soil_data$Date_Time <- as.POSIXct(soil_data$sample_datetime, format = "%Y-%m-%d %H:%M:%S")


# adding month, day, year, and hour columns
hobo_data$month <- format(hobo_data$Date_Time,format="%m")
hobo_data$year <- format(hobo_data$Date_Time,format="%Y")
hobo_data$hour <- format(hobo_data$Date_Time, format="%H")
hobo_data$day <- format(hobo_data$Date_Time, format="%d")
soil_data$month <- format(soil_data$Date_Time,format="%m")
soil_data$year <- format(soil_data$Date_Time,format="%Y")
soil_data$hour <- format(soil_data$Date_Time, format="%H")
soil_data$day <- format(soil_data$Date_Time, format="%d")


# taking 2022 temp average for july 11 - july 15 (sampling period)
hobo_data$date <- paste0(hobo_data$month,"",hobo_data$day)
hobo_data$date <- as.numeric(hobo_data$date)
# selecting 2022 and june 1 - july 15
hobo_sampling <- hobo_data %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
hobo_sampling <- hobo_sampling %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)


# taking 2022 soil average for july 11 - july 15 (sampling period)
soil_data$date <- paste0(soil_data$month,"",soil_data$day)
soil_data$date <- as.numeric(soil_data$date)
# selecting 2022 and june 1 - july 15
soil_sampling <- soil_data %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
soil_sampling <- soil_sampling %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)


# HOBO model exploration
descdist(hobo_sampling$Temperature_C, discrete = FALSE)
hist(hobo_sampling$Temperature_C)
qqnorm(hobo_sampling$Temperature_C)
shapiro.test(hobo_sampling$Temperature_C)
# kinda right skewed, lets see once its in a model
hobo_sampling <- within(hobo_sampling, Treatment <- relevel(factor(Treatment), ref = "Warmed_Drought"))
m.hobo.test <- lmer(Temperature_C ~ Treatment + (1|Rep), data = hobo_sampling, REML=FALSE)
anova(m.hobo.test)
summary(m.hobo.test)
hist(resid(m.hobo.test))
shapiro.test(resid(m.hobo.test)) # still not normal, but going with this model
descdist(resid(m.hobo.test), discrete = FALSE)
emmeans(m.hobo.test, list(pairwise ~ Treatment), adjust = "tukey")





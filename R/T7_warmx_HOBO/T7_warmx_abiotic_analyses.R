# TITLE:          REX: Abiotic analyses
# AUTHORS:        Kara Dobson, Moriah Young
# COLLABORATORS:  Phoebe Zarnetske, Mark Hammond, Emily Parker
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

hobo_data$date <- paste0(hobo_data$month,"",hobo_data$day)
hobo_data$date <- as.numeric(hobo_data$date)
soil_data$date <- paste0(soil_data$month,"",soil_data$day)
soil_data$date <- as.numeric(soil_data$date)

###############################################################################
# Below is for Kara's VOC sampling in 2022
# taking 2022 temp average for july 11 - july 15 (VOC sampling period)
# selecting 2022 and june 1 - july 15
hobo_sampling_VOC <- hobo_data %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
hobo_sampling_VOC <- hobo_sampling_VOC %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)


# taking 2022 soil average for july 11 - july 15 (VOC sampling period)
# selecting 2022 and june 1 - july 15
soil_sampling_VOC <- soil_data %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
soil_sampling_VOC <- soil_sampling_VOC %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)


# HOBO model exploration for air temperature
descdist(hobo_sampling_VOC$Temperature_C, discrete = FALSE)
hist(hobo_sampling_VOC$Temperature_C)
qqnorm(hobo_sampling_VOC$Temperature_C)
shapiro.test(hobo_sampling_VOC$Temperature_C)
# kinda right skewed, lets see once its in a model
#hobo_sampling <- within(hobo_sampling, Treatment <- relevel(factor(Treatment), ref = "Warmed_Drought"))
m.hobo.test <- lmer(Temperature_C ~ Treatment + (1|Rep), data = hobo_sampling_VOC, REML=FALSE)
anova(m.hobo.test)
summary(m.hobo.test)
hist(resid(m.hobo.test))
shapiro.test(resid(m.hobo.test)) # still not normal, but going with this model
descdist(resid(m.hobo.test), discrete = FALSE)
emmeans(m.hobo.test, list(pairwise ~ Treatment), adjust = "tukey")
################################################################################

# select just for 2022 data
hobo_sampling_2022 <- hobo_data %>%
        filter(year == 2022)
soil_sampling_VOC <- soil_data %>%
        filter(year == 2022)

# HOBO model exploration for air temperature
descdist(hobo_sampling_2022$Temperature_C, discrete = FALSE)
hist(hobo_sampling_2022$Temperature_C)
qqnorm(hobo_sampling_2022$Temperature_C)
shapiro.test(hobo_sampling_2022$Temperature_C)
m.hobo.test <- lmer(Temperature_C ~ Treatment + (1|Rep), data = hobo_sampling_2022, REML=FALSE)
anova(m.hobo.test)
summary(m.hobo.test)
hist(resid(m.hobo.test))
shapiro.test(resid(m.hobo.test))
descdist(resid(m.hobo.test), discrete = FALSE)
emmeans(m.hobo.test, list(pairwise ~ Treatment), adjust = "tukey")



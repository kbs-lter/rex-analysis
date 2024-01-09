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
# air
hobo_data$month <- format(hobo_data$Date_Time,format="%m")
hobo_data$year <- format(hobo_data$Date_Time,format="%Y")
hobo_data$hour <- format(hobo_data$Date_Time, format="%H")
hobo_data$day <- format(hobo_data$Date_Time, format="%d")
# soil
soil_data$month <- format(soil_data$Date_Time,format="%m")
soil_data$year <- format(soil_data$Date_Time,format="%Y")
soil_data$hour <- format(soil_data$Date_Time, format="%H")
soil_data$day <- format(soil_data$Date_Time, format="%d")

# create date column that combines month and day together, i.e. 09 (month) and 28 (day) = 928 (sept 28th)
hobo_data$date <- paste0(hobo_data$month,"",hobo_data$day)
hobo_data$date <- as.numeric(hobo_data$date)
soil_data$date <- paste0(soil_data$month,"",soil_data$day)
soil_data$date <- as.numeric(soil_data$date)

###############################################################################
# Below is for Kara's VOC sampling in 2022
# taking 2022 temps for july 11 - july 15 (VOC sampling period)
hobo_sampling_VOC <- hobo_data %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
hobo_sampling_VOC <- hobo_sampling_VOC %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)


# taking 2022 soil for july 11 - july 15 (VOC sampling period)
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
# check normality
hist(resid(m.hobo.test))
shapiro.test(resid(m.hobo.test)) # still not normal, but going with this model
descdist(resid(m.hobo.test), discrete = FALSE)
# results
anova(m.hobo.test)
summary(m.hobo.test)
contrast(emmeans(m.hobo.test, ~Treatment), "pairwise", simple = "each", combine = F, adjust = "mvt")


# soil model exploration for soil temperature
descdist(soil_sampling_VOC$temperature, discrete = FALSE)
hist(soil_sampling_VOC$temperature)
qqnorm(soil_sampling_VOC$temperature)
shapiro.test(soil_sampling_VOC$temperature)
# kinda right skewed, lets see once its in a model
#soil_sampling <- within(soil_sampling, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "Warmed_Drought"))
m.soil.test <- lmer(temperature ~ Subplot_Descriptions + (1|Rep), data = soil_sampling_VOC, REML=FALSE)
# assumption checking
hist(resid(m.soil.test))
shapiro.test(resid(m.soil.test)) # still not normal, but going with this model
descdist(resid(m.soil.test), discrete = FALSE)
# results
anova(m.soil.test)
summary(m.soil.test)
contrast(emmeans(m.soil.test, ~Subplot_Descriptions), "pairwise", simple = "each", combine = F, adjust = "mvt")


# soil model exploration for soil moisture
descdist(soil_sampling_VOC$vwc, discrete = FALSE)
hist(soil_sampling_VOC$vwc)
qqnorm(soil_sampling_VOC$vwc)
shapiro.test(soil_sampling_VOC$vwc)
# kinda right skewed, lets see once its in a model
#soil_sampling <- within(soil_sampling, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "Warmed_Drought"))
m.soil.moist.test <- lmer(vwc ~ Subplot_Descriptions + (1|Rep), data = soil_sampling_VOC, REML=FALSE)
# assumption checking
hist(resid(m.soil.moist.test))
shapiro.test(resid(m.soil.moist.test)) # still not normal, but going with this model
descdist(resid(m.soil.moist.test), discrete = FALSE)
# results
anova(m.soil.moist.test)
summary(m.soil.moist.test)
contrast(emmeans(m.soil.moist.test, ~Subplot_Descriptions), "pairwise", simple = "each", combine = F, adjust = "mvt")


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

################################################################################

# HOBO model exploration for air temperature across all year (all HOBO data available)
descdist(hobo_data$Temperature_C, discrete = FALSE) # this doesn't work, error
hist(hobo_data$Temperature_C)
qqnorm(hobo_data$Temperature_C)
shapiro.test(hobo_data$Temperature_C)
m.hobo.test.all <- lmer(Temperature_C ~ Treatment + year + (1|Rep), data = hobo_data, REML=FALSE)
anova(m.hobo.test.all)
summary(m.hobo.test.all)
hist(resid(m.hobo.test.all))
shapiro.test(resid(m.hobo.test.all))
descdist(resid(m.hobo.test.all), discrete = FALSE)
emmeans(m.hobo.test.all, list(pairwise ~ Treatment), adjust = "tukey")
#$`pairwise differences of Treatment`
#1                        estimate     SE  df z.ratio p.value
#Ambient - Drought          -0.054 0.0569 Inf  -0.948  0.7789
#Ambient - Warmed           -1.338 0.0561 Inf -23.825  <.0001
#Ambient - Warmed_Drought   -0.996 0.0570 Inf -17.481  <.0001
#Drought - Warmed           -1.284 0.0576 Inf -22.301  <.0001
#Drought - Warmed_Drought   -0.942 0.0577 Inf -16.317  <.0001
#Warmed - Warmed_Drought     0.342 0.0576 Inf   5.935  <.0001

emmeans(m.hobo.test.all, list(pairwise ~ year), adjust = "tukey")

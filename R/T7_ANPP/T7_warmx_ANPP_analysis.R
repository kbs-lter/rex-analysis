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
raw.data.test <- lm(plant_biomass_gm2 ~ Drought*Warming*Insecticide*Year, data=anpp1) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))

# log transformation
anpp1$log_biomass <- log(anpp1$plant_biomass_gm2)
descdist(anpp1$log_biomass, discrete = FALSE)
hist(anpp1$log_biomass)
qqnorm(anpp1$log_biomass) # a little weird
shapiro.test(anpp1$log_biomass) # still not normal
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(plant_biomass_gm2 ~ Drought*Warming*Insecticide*Year, data=anpp1)# testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # not normal
# better - going with log transformation


######### Assumption checking ##########
# in the assumption checking, we're making sure that our full model meets the assumptions of the model.
# the model below is the model structure we can use for all response variables; its testing to see
# if there is 1. an effect of climate treatment on height, 2. an effect of galling status on height, and 3. does the effect
# of climate on height depend on galling status. Subplot nested within footprint nested within rep is used as our random effect
# to account for variation between plots. Year is also included as a random effect to account for variation between years.
m1 <- lmer(log_biomass ~ Drought*Warming*Insecticide*as.factor(Year) + (1|Rep/Footprint/Subplot), data = anpp1, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Biomass")
# Homogeneity of variance looks a bit off (increasing variance in resids does increase with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ anpp1$Subplot_Descriptions) # met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Biomass")
hist(residuals(m1), main = "Biomass")
shapiro.test(resid(m1))
outlierTest(m1) # checking for outliers - none

###### Checking model results ########
anova(m1)
summary(m1)

contrast(emmeans(m1, ~Subplot_Descriptions + as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")

################################################################################

# filter for 2023 data and 1.0m2 data
anpp23 <- anpp %>%
        filter(Year=='2023') %>%
        filter(Scale_meter_square==1)

# re-leveling the dataframe
# first make into factors
anpp23$Drought <-as.factor(anpp23$Drought)
anpp23$Warming <-as.factor(anpp23$Warming)
anpp23$Insecticide <-as.factor(anpp23$Insecticide)

anpp23 <- within(anpp23, Drought <- relevel(factor(Drought), ref = "drought"))
anpp23 <- within(anpp23, Warming <- relevel(factor(Warming), ref = "warmed"))
anpp23 <- within(anpp23, Insecticide <- relevel(factor(Insecticide), ref = "insecticide"))

###### Data exploration #######
descdist(anpp23$plant_biomass_gm2, discrete = FALSE) # checks what type of distribution your data matches
hist(anpp23$plant_biomass_gm2)
qqnorm(anpp23$plant_biomass_gm2) # should be a straight diagonal line if normal
shapiro.test(anpp23$plant_biomass_gm2) # p > 0.05 if normal
# making a model to test distribution of model residuals
raw.data.test <- lm(plant_biomass_gm2 ~ Drought, data=anpp23) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))

# log transformation
anpp23$log_biomass <- log(anpp23$plant_biomass_gm2)
descdist(anpp23$log_biomass, discrete = FALSE)
hist(anpp23$log_biomass)
qqnorm(anpp23$log_biomass) # a little weird
shapiro.test(anpp23$log_biomass) # still not normal
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(log_biomass ~ Subplot_Descriptions * Year, data=anpp23) # testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # not normal
# better - going with log transformation


######### Assumption checking ##########
# in the assumption checking, we're making sure that our full model meets the assumptions of the model.
# the model below is the model structure we can use for all response variables; its testing to see
# if there is 1. an effect of climate treatment on height, 2. an effect of galling status on height, and 3. does the effect
# of climate on height depend on galling status. Subplot nested within footprint nested within rep is used as our random effect
# to account for variation between plots. Year is also included as a random effect to account for variation between years.
m23_full <- lmer(log_biomass ~ Warming*Drought*Insecticide + (1|Rep/Footprint), data = anpp23, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m23_full, main = "Biomass")
# Homogeneity of variance looks a bit off (increasing variance in resids does increase with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m23_full) ~ anpp23$Drought) # met
leveneTest(residuals(m23_full) ~ anpp23$Warming) # met
leveneTest(residuals(m23_full) ~ anpp23$Insecticide) # met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m23), main = "Biomass")
hist(residuals(m23), main = "Biomass")
shapiro.test(resid(m23))
outlierTest(m23) # checking for outliers - none

###### Checking model results ########
anova(m23_full)
summary(m23_full)

contrast(emmeans(m1, ~Subplot_Descriptions + as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")
# TITLE:          REX: Soca plant height analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    analyses
# PROJECT:        REX
# DATE:           Jan. 2022; updated Oct. 2023

# Clear all existing data
rm(list=ls())

# Load packages
library(lmerTest)
library(fitdistrplus)
library(sjPlot)
library(tidyverse)
library(car)
library(emmeans)
library(bbmle)
library(multcomp)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
height <- read.csv(file.path(dir, "/T7_warmx_plant_traits/L1/T7_warmx_soca_height_harvest_L1.csv"))


###### Data exploration #######
# The first thing we do is check the distribution of the raw data
# We're checking for a normal distribution; the distribution of the raw data
# doesn't matter as much as the distribution of a model's residuals (which we'll also check),
# but its a good way to get a glimpse of how the data look
descdist(height$Height_cm, discrete = FALSE) # checks what type of distribution your data matches
hist(height$Height_cm) # should be a bell curve if normal
qqnorm(height$Height_cm) # should be a straight diagonal line if normal
shapiro.test(height$Height_cm) # p > 0.05 if normal
# making a model to test distribution of model residuals
raw.data.test <- lm(Height_cm ~ Climate_Treatment * Galling_Status, data=height) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))
# looks pretty good, going to try some transformations but might not need them

# log transformation
height$log_height <- log(height$Height_cm)
descdist(height$log_height, discrete = FALSE)
hist(height$log_height)
qqnorm(height$log_height)
shapiro.test(height$log_height) # looks good
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(log_height ~ Climate_Treatment * Galling_Status, data=height) # testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log))
# better - going with log transformation



######### Assumption checking ##########
# in the assumption checking, we're making sure that our full model meets the assumptions of the model.
# the model below is the model structure we can use for all response variables; its testing to see
# if there is 1. an effect of climate treatment on height, 2. an effect of galling status on height, and 3. does the effect
# of climate on height depend on galling status. Subplot nested within footprint nested within rep is used as our random effect
# to account for variation between plots. Year is also included as a random effect to account for variation between years.
m1 <- lmer(log_height ~ Climate_Treatment * Galling_Status + (1|Rep/Footprint/Subplot) + (1|Year), data = height, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Plant height")
# Homogeneity of variance looks a bit off (increasing variance in resids does increase with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ height$Climate_Treatment) # met
leveneTest(residuals(m1) ~ height$Galling_Status) # not met, but pretty close
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Plant height")
hist(residuals(m1), main = "Plant height")
shapiro.test(resid(m1))
outlierTest(m1) # checking for outliers - none



###### Checking model results ########
anova(m1)
emmip(m1, Climate_Treatment~Galling_Status)
# this outcome shows us that climate has an effect but not galling. there is also no interaction btwn climate and galling
# model w/o interaction term, since it was not significant
height <- within(height, Climate_Treatment <- relevel(factor(Climate_Treatment), ref = "Warm Drought")) # re-set the reference level for diff. comparisons
m2 <- lmer(log(Height_cm) ~ Climate_Treatment + Galling_Status + (1|Rep/Footprint/Subplot) + (1|Year), data = height, REML=F)
anova(m2)
summary(m2)
# back-transforming - W vs. A and W vs. I
exp(4.43673)-exp(4.43673-0.22449) # warmed plants 17 cm taller than ambient
exp(4.43673)-exp(4.43673-0.19599) # warmed plants 15 cm taller than irrigated
# back-transforming - WD vs. A and WD vs. D (set the reference level to warm drought instead of warm)
exp(4.42690)-exp(4.42690-0.21466) # warmed drought plants 16 cm taller than ambient
exp(4.42690)-exp(4.42690-0.18240) # warmed drought plants 14 cm taller than drought
contrast(emmeans(m2, ~Climate_Treatment), "pairwise", simple = "each", combine = F, adjust = "mvt")
(0.799 - 1) * 100 # 20% decrease from ambient to warmed

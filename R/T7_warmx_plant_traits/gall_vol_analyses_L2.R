# TITLE:          REX: Gall volume analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    analyses
# PROJECT:        REX
# DATE:           July 2021

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(lmerTest)
library(fitdistrplus)
library(sjPlot)
library(tidyverse)
library(car)
library(emmeans)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
galls <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_gall_vol_L1.csv"))

# Data exploration
descdist(galls$sphere_vol, discrete = FALSE)
hist(galls$sphere_vol)
qqnorm(galls$sphere_vol)
shapiro.test(galls$sphere_vol)
# pretty right skewed, going to try a few transformations

# square root transformation
galls$sqrt_vol <- sqrt(galls$sphere_vol)
descdist(galls$sqrt_vol, discrete = FALSE)
hist(galls$sqrt_vol)
qqnorm(galls$sqrt_vol)
shapiro.test(galls$sqrt_vol)
# pretty good

# Assumption checking
m1 <- lmer(sqrt_vol ~ treatment + (1|rep), data = galls, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Gall volume")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ galls$treatment)
# Assumption met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Gall volume")
hist(residuals(m1), main = "Gall volume")
shapiro.test(resid(m1)) # Good

# comparison with other models
m2 <- lmer(sqrt_vol ~ treatment + (1|rep/footprint), data = galls, REML=FALSE)
AICctab(m1, m2, weights=T) # going with model 2 because of nested design

summary(m2)
emmeans(m2, list(pairwise ~ treatment), adjust = "tukey")

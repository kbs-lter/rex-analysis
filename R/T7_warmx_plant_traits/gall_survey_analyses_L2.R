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

# Gall count data #
# Data exploration
descdist(gall$total_rosette_quadrat, discrete = FALSE)
hist(gall$total_rosette_quadrat)
qqnorm(gall$total_rosette_quadrat)
shapiro.test(gall$total_rosette_quadrat)
# pretty right skewed, going to try a few transformations

# square root transformation
gall$sqrt_total_rose <- sqrt(gall$total_rosette_quadrat)
descdist(gall$sqrt_total_rose, discrete = FALSE)
hist(gall$sqrt_total_rose)
qqnorm(gall$sqrt_total_rose)
shapiro.test(gall$sqrt_total_rose)
# better, not perfect

# cubed root transformation
gall$cubed_total_rose <- (gall$total_rosette_quadrat)^(1/3)
descdist(gall$cubed_total_rose, discrete = FALSE)
hist(gall$cubed_total_rose)
qqnorm(gall$cubed_total_rose)
shapiro.test(gall$cubed_total_rose)
# pretty good

# Assumption checking
m1 <- lmer(cubed_total_rose ~ treatment + (1|rep), data = gall, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Total rosette count")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ gall$treatment)
# Assumption met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Total rosette count")
hist(residuals(m1), main = "Total rosette count")
shapiro.test(resid(m1)) # Good

# comparison with other models
m2 <- lmer(cubed_total_rose ~ treatment + (1|rep/footprint_number), data = gall, REML=FALSE)
AICctab(m1, m2, weights=T)

summary(m2)
emmeans(m2, list(pairwise ~ treatment), adjust = "tukey")



# Galling proportion data #
gall <- gall %>%
  mutate(prop_gall = rosette_galls/total_stems)
# Data exploration
descdist(gall$prop_gall, discrete = FALSE)
hist(gall$prop_gall)
qqnorm(gall$prop_gall)
shapiro.test(gall$prop_gall)
# right skewed, going to try a few transformations

# square root transformation
gall$sqrt_prop <- sqrt(gall$prop_gall)
descdist(gall$sqrt_prop, discrete = FALSE)
hist(gall$sqrt_prop)
qqnorm(gall$sqrt_prop)
shapiro.test(gall$sqrt_prop)
# better, not perfect

# cubed root transformation
gall$cubed_prop <- (gall$prop_gall)^(1/3)
descdist(gall$cubed_prop, discrete = FALSE)
hist(gall$cubed_prop)
qqnorm(gall$cubed_prop)
shapiro.test(gall$cubed_prop)
# sqrt was better

# Assumption checking
m3 <- lmer(sqrt_prop ~ treatment + (1|rep), data = gall, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m3, main = "Proportion galled")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m3) ~ gall$treatment)
# Error, go back to this another time
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m3), main = "Proportion galled")
hist(residuals(m3), main = "Proportion galled")
shapiro.test(resid(m3)) # Good enough

# comparison with other models
m4 <- lmer(sqrt_prop ~ treatment + (1|rep/footprint_number), data = gall, REML=FALSE)
AICctab(m3, m4, weights=T)

summary(m4)
emmeans(m4, list(pairwise ~ treatment), adjust = "tukey")

# emmeans is a lot different than summary (?) so re-leveling to get pairwise comparisons
gall <- within(gall, treatment <- relevel(factor(treatment), ref = "warmed"))
summary(m4)
gall <- within(gall, treatment <- relevel(factor(treatment), ref = "drought"))
summary(m4)
gall <- within(gall, treatment <- relevel(factor(treatment), ref = "irrigated_control"))
summary(m4)


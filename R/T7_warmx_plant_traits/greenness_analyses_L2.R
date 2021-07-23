# TITLE:          REX: Greenness analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Mark Hammond
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
green <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_greenness_L1.csv"))

# Take subplot average
green2 <- green %>%
  group_by(rep, footprint, treatment,gall_present) %>%
  summarize(greenness = mean(greenness, na.rm = TRUE))

# Data exploration
descdist(green2$greenness, discrete = FALSE)
hist(green2$greenness)
qqnorm(green2$greenness)
shapiro.test(green2$greenness)
# Normally distributed

# Assumption checking
m1 <- lmer(greenness ~ treatment + gall_present + (1|rep), data = green2, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Greenness")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ green2$treatment)
# Assumption met
leveneTest(residuals(m1) ~ green2$gall_present)
# Assumption met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Greenness")
hist(residuals(m1), main = "Greenness")
shapiro.test(resid(m1)) # Normal

# Model comparisons
m2 <- lm(greenness ~ treatment, data=green2)
m3 <- lm(greenness ~ gall_present, data=green2)
m4 <- lmer(greenness ~ treatment + (1|rep), data=green2, REML=F)
m5 <- lmer(greenness ~ treatment * gall_present + (1|rep), data=green2, REML=F)
AICctab(m1, m2, m3, m4, m5, weights=T)
# Model 1 fits the best
summary(m1)

# Post hoc test to compare different levels
emmeans(m1, list(pairwise ~ treatment + gall_present), adjust = "tukey")

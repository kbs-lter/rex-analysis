# TITLE:          REX: Soca plant height analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    analyses
# PROJECT:        REX
# DATE:           Jan 2022

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
height <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_Soca_plant_height_L1.csv"))

# Data exploration
descdist(height$plant_height_cm, discrete = FALSE)
hist(height$plant_height_cm)
qqnorm(height$plant_height_cm)
shapiro.test(height$plant_height_cm)
# looks pretty normal, going to try some transformations but might not need them

# square root transformation
height$sqrt_height <- sqrt(height$plant_height_cm)
descdist(height$sqrt_height, discrete = FALSE)
hist(height$sqrt_height)
qqnorm(height$sqrt_height)
shapiro.test(height$sqrt_height) # looks good

# log transformation
height$log_height <- log(height$plant_height_cm)
descdist(height$log_height, discrete = FALSE)
hist(height$log_height)
qqnorm(height$log_height)
shapiro.test(height$log_height) # slightly better than sqrt

# Assumption checking - log transformation
m1 <- lmer(log_height ~ treatment + drought_period + gall_present + (1|rep), data = height, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Plant height")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ height$treatment)
leveneTest(residuals(m1) ~ height$drought_period)
leveneTest(residuals(m1) ~ height$gall_present)
# Assumption not met - ignoring for now
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Plant height")
hist(residuals(m1), main = "Plant height")
shapiro.test(resid(m1))
outlierTest(m1)

# comparison with other models
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
m3 <- lmer(log_height ~ treatment + (1|drought_period) + (1|rep), data = height, REML=F)
m4 <- lmer(log_height ~ treatment * drought_period + (1|rep), data = height, REML=F)
m5 <- lmer(log_height ~ treatment + drought_period + (1|rep), data = height, REML=F)
m6 <- lmer(log_height ~ treatment * gall_present + drought_period + (1|rep), data = height, REML=F)
m7 <- lmer(log_height ~ treatment * gall_present + drought_period + (1|rep/footprint), data = height, REML=F)
m8 <- lmer(log_height ~ treatment + gall_present + drought_period + (1|rep/footprint), data = height, REML=F)
AICctab(m1, m2, m3, m4, m5,m6,m7,m8,weights=T) # model 7
leveneTest(residuals(m7) ~ height$treatment)
leveneTest(residuals(m7) ~ height$drought_period)
summary(m7)
emmeans(m7, list(pairwise ~ treatment*gall_present), adjust = "tukey")
png("height_galling_lm.png", units="in", width=7, height=6, res=300)
emmip(m7, treatment ~ gall_present)
dev.off()

## note: these comparisons below need to be changed for m7
# re-leveling the data & checking summary
height <- within(height, treatment <- relevel(factor(treatment), ref = "Ambient Drought"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)
height <- within(height, treatment <- relevel(factor(treatment), ref = "Irrigated Control"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)
height <- within(height, treatment <- relevel(factor(treatment), ref = "Warm"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)
height <- within(height, treatment <- relevel(factor(treatment), ref = "Warm Drought"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)
height <- within(height, treatment <- relevel(factor(treatment), ref = "Ambient"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)

height <- within(height, drought_period <- relevel(factor(drought_period), ref = "Post-Drought"))
m2 <- lmer(log_height ~ treatment + drought_period + (1|rep/footprint), data = height, REML=F)
summary(m2)

## note: these estimates below are slightly different now for m7 - didn't update after adding in interaction term
# calculating effect size accounting for log
exp(4.173e+00 + 1.954e-01*0) # 64.90989 - average for ambient
exp(4.173e+00 + 1.954e-01*1) # 78.91726 - average for warmed
# effect:
78.91726 - 64.90989 # 14.00737cm taller plants in warmed compared to ambient

# calculating effect size accounting for log
exp(4.173e+00 + 1.899e-01*0) # 64.90989 - average for ambient
exp(4.173e+00 + 1.899e-01*1) # 78.48441 - average for warmed drought
# effect:
78.48441 - 64.90989 # 13.57452cm taller plants in warmed drought compared to ambient

# calculating effect size accounting for log
exp(4.166e+00 + 2.023e-01*0) # 64.45711 - average for drought
exp(4.166e+00 + 2.023e-01*1) # 78.90937 - average for warmed
# effect:
78.90937 - 64.45711  # 14.45226cm taller plants in warmed compared to drought

# calculating effect size accounting for log
exp(4.166e+00 + 1.968e-01*0) # 64.45711 - average for drought
exp(4.166e+00 + 1.968e-01*1) # 78.47656 - average for warmed drought
# effect:
78.47656 - 64.45711  # 14.01945cm taller plants in warmed drought compared to drought

# calculating effect size accounting for log
exp(4.12117 + 0.24760*0) # 61.63131 - average for irr control
exp(4.12117 + 0.24760*1) # 78.94647 - average for warmed
# effect:
78.94647 - 61.63131  # 17.31516cm taller plants in warmed compared to irr control

# calculating effect size accounting for log
exp(4.12117 + 0.24208 *0) # 61.63131 - average for irr control
exp(4.12117 + 0.24208 *1) # 78.51188 - average for warmed
# effect:
78.51188 - 61.63131  # 16.88057cm taller plants in warmed drought compared to irr control

# this is updated for m7 values
# calculating effect size accounting for log
exp(4.14017 + 0.13394 *0) # 62.8135 - average for gall
exp(4.14017 + 0.13394 *1) # 71.81619 - average for no gall
# effect:
71.81619 - 62.8135  # 9.00269cm taller non galled plants compared to galled plants


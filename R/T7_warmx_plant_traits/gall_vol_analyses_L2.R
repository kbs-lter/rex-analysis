# TITLE:          REX: Gall volume analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    analyses
# PROJECT:        REX
# DATE:           July 2021; updated Jan 2022

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
galls <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_Soca_gall_vol_L1.csv"))

# Data exploration
descdist(galls$sphere_vol_cm3, discrete = FALSE)
hist(galls$sphere_vol_cm3)
qqnorm(galls$sphere_vol_cm3)
shapiro.test(galls$sphere_vol_cm3)
# pretty right skewed, going to try a few transformations & gamma distribution

# square root transformation
galls$sqrt_vol <- sqrt(galls$sphere_vol_cm3)
descdist(galls$sqrt_vol, discrete = FALSE)
hist(galls$sqrt_vol)
qqnorm(galls$sqrt_vol)
shapiro.test(galls$sqrt_vol)

# log transformation
galls$log_vol <- log(galls$sphere_vol_cm3)
descdist(galls$log_vol, discrete = FALSE)
hist(galls$log_vol)
qqnorm(galls$log_vol)
shapiro.test(galls$log_vol) # slightly better than sqrt

# gamma distribution?
fit.gamma <- fitdist(galls$sphere_vol_cm3, "gamma")
plot(fit.gamma) # this looks pretty good

# Assumption checking - log transformation
m1 <- lmer(log_vol ~ treatment + drought_period + (1|rep), data = galls, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Gall volume")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ galls$treatment)
leveneTest(residuals(m1) ~ galls$drought_period)
# Assumption met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Gall volume")
hist(residuals(m1), main = "Gall volume")
shapiro.test(resid(m1))
outlierTest(m1)

# Assumption checking - sqrt transformation
m2 <- lmer(sqrt_vol ~ treatment + drought_period + (1|rep), data = galls, REML=FALSE)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m2, main = "Gall volume")
# Homogeneity of variance is ok here (increasing variance in resids is not increasing with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m2) ~ galls$treatment)
# Assumption met
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m2), main = "Gall volume")
hist(residuals(m2), main = "Gall volume")
shapiro.test(resid(m2))
outlierTest(m2)

# Assumption checking - gamma distribution (do these apply to non-normal distributions?)
m3 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep), data = galls, family = Gamma)
plot(m3, main = "Gall volume")
leveneTest(residuals(m3) ~ galls$treatment)
leveneTest(residuals(m3) ~ galls$drought_period)
qqPlot(resid(m3), main = "Gall volume")
hist(residuals(m3), main = "Gall volume")
shapiro.test(resid(m3)) # gonna go with gamma distribution

# comparison with other models
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
m5 <- glmer(sphere_vol_cm3 ~ treatment + (1|drought_period) + (1|rep/footprint), data = galls, family = Gamma)
m6 <- glmer(sphere_vol_cm3 ~ treatment * drought_period + (1|rep/footprint), data = galls, family = Gamma)
AICctab(m3, m4, m5, m6, weights=T) # model 4

summary(m4)
emmeans(m4, list(pairwise ~ treatment), adjust = "tukey")
summary(glht(m4, mcp(treatment="Tukey")), test = adjusted(type = "bonferroni"))

# re-leveling the data & checking summary since post-hoc gives diff results (?)
galls <- within(galls, treatment <- relevel(factor(treatment), ref = "Ambient Drought"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)
galls <- within(galls, treatment <- relevel(factor(treatment), ref = "Irrigated Control"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)
galls <- within(galls, treatment <- relevel(factor(treatment), ref = "Warm"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)
galls <- within(galls, treatment <- relevel(factor(treatment), ref = "Warm Drought"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)
galls <- within(galls, treatment <- relevel(factor(treatment), ref = "Ambient"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)

galls <- within(galls, drought_period <- relevel(factor(drought_period), ref = "Pre-Drought"))
m4 <- glmer(sphere_vol_cm3 ~ treatment + drought_period + (1|rep/footprint), data = galls, family = Gamma)
summary(m4)

# calculating effect size accounting for inverse link
1/(0.025608 + 0.018283*0) # 39.0503 - average for irrigated control
1/(0.025608 + 0.018283*1) # 22.78371 - average for warmed
# effect:
22.78371 - 39.0503 # 16.26659 cm smaller galls in warmed compared to irrigated

# calculating effect size accounting for inverse link
1/(0.043890  + -0.016302*0) # 22.78423 - average for warmed
1/(0.043890 + -0.016302*1) # 36.24764 - average for drought
# effect:
36.24764 - 22.78423 # 13.46341 cm smaller galls in warmed compared to drought

# calculating effect size accounting for inverse link
1/(0.043890  + -0.013871*0) # 22.78423 - average for warmed
1/(0.043890 + -0.013871*1) # 33.31224 - average for warmed drought
# effect:
33.31224 - 22.78423 # 10.52801cm smaller galls in warmed compared to warmed drought

# calculating effect size accounting for inverse link
1/(0.043890  + 0.004594*0) # 22.78423 - average for drought
1/(0.043890 + 0.004594*1) # 20.62536 - average for post-drought
# effect:
20.62536 - 22.78423 # 2.15887cm smaller galls in post-drought compared to drought

# calculating effect size accounting for inverse link
1/(0.043890  + 0.014814*0) # 22.78423 - average for drought
1/(0.043890 + 0.014814*1) # 17.03461 - average for pre-drought
# effect:
17.03461 - 22.78423 # 5.74962cm smaller galls in pre-drought compared to drought

# calculating effect size accounting for inverse link
1/(0.058703  + -0.010220*0) # 17.0349 - average for pre-drought
1/(0.058703 + -0.010220*1) # 20.62579 - average for post-drought
# effect:
20.62579 - 17.0349 # 3.59089cm smaller galls in pre-drought compared to post-drought

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
library(fitdistrplus)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Read in data
anpp <- read.csv(file.path(dir, "T7_ANPP/L1/T7_warmx_ANPP_L1.csv")) # 2019-2023 ANPP data

# paste
anpp$footprint_unique <- paste(anpp$Replicate, anpp$FP_location)
anpp$subplot_unique <- paste(anpp$Replicate, anpp$FP_location, anpp$Subplot_location)


# Convert relevant columns to a factor
anpp[12:15] <- lapply(anpp[12:15], as.factor)
anpp[1] <- lapply(anpp[1], as.factor)

# filter out 2019 data and non 0.2 data
anpp1 <- anpp %>%
        filter(!(Year=='2019')) %>%
        filter(Scale_meter_square==0.2)

# filter out insecticide plots
# non insecticide, no 2019 data, and just 0.2 scale
anpp_noinsect <- anpp1 %>%
        filter(!(Insecticide =='insecticide'))

###############################################################################
# calculate 0.2 data to 1m2
# Multiply all values in the column by 5
anpp2 <- anpp1 %>%
        filter(!(Year==2023))
anpp2$scaled_biomass <- anpp2$plant_biomass_gm2 * 5
anpp2 <- anpp2[,-c(16,17)]
names(anpp2)[names(anpp2)=="scaled_biomass"] <- "plant_biomass_gm2" 
anpp2$Scale_meter_square <- 1

anpp19 <- anpp %>% filter(Year==2019)

anpp23_1 <- anpp %>%
        filter(Year==2023) %>% 
        filter(Scale_meter_square==1)

anpp_merge <- full_join(anpp19, anpp23_1)
anpp_scaled <- full_join(anpp_merge, anpp2)

###############################################################################
# Using the warming, drought, and insecticide columns rather than the subplot_descriptions column

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
raw.data.test.log <- lm(log_biomass ~ Drought*Warming*Insecticide*Year, data=anpp1)# testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # still not normal but plot looks a lot better
# better - going with log transformation

######### Assumption checking ##########
# in the assumption checking, we're making sure that our full model meets the assumptions of the model.
# the model below is the model structure we can use for all response variables; its testing to see
# if there is 1. an effect of climate treatment on height, 2. an effect of galling status on height, and 3. does the effect
# of climate on height depend on galling status. Subplot nested within footprint nested within rep is used as our random effect
# to account for variation between plots. Year is also included as a random effect to account for variation between years.
m1 <- lmer(log_biomass ~ Drought*Warming*Insecticide*as.factor(Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m1, main = "Biomass")
# Homogeneity of variance looks a bit off (increasing variance in resids does increase with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m1) ~ anpp1$Drought) # met
leveneTest(residuals(m1) ~ anpp1$Warming) 
leveneTest(residuals(m1) ~ anpp1$Insecticide) 
leveneTest(residuals(m1) ~ anpp1$Year) 
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m1), main = "Biomass")
hist(residuals(m1), main = "Biomass")
shapiro.test(resid(m1))
outlierTest(m1) # checking for outliers - none

###### Checking model results ########
anova(m1)
summary(m1)

contrast(emmeans(m1, ~ Drought + as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")
contrast(emmeans(m1, ~ Warming + as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")

m1 <- lmer(log_biomass ~ Drought*Warming*Insecticide*as.factor(Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
m2 <- lmer(log_biomass ~ Drought*Warming*as.factor(Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
m3 <- lmer(log_biomass ~ Drought*Warming + (1|Rep/footprint_unique/subplot_unique) + (1|Year) , data = anpp1, REML=F)
m4 <- lmer(log_biomass ~ Drought + Warming + Insecticide + Year + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)

anova(m1)
summary(m1)
anova(m2)
summary(m2)
anova(m3)
summary(m3)
anova(m4)
summary(m4)

################################################################################
# Using the "Subplot_Descriptions" column

###### Data exploration #######
descdist(anpp1$plant_biomass_gm2, discrete = FALSE)
hist(anpp1$plant_biomass_gm2)
qqnorm(anpp1$plant_biomass_gm2) # should be a straight diagonal line if normal
shapiro.test(anpp1$plant_biomass_gm2) # p > 0.05 if normal
# making a model to test distribution of model residuals
raw.data.test <- lm(plant_biomass_gm2 ~ Subplot_Descriptions*as.factor(Year), data=anpp1) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))

# log transformation
descdist(anpp1$log_biomass, discrete = FALSE)
hist(anpp1$log_biomass)
qqnorm(anpp1$log_biomass) # a little weird
shapiro.test(anpp1$log_biomass) # still not normal
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(log_biomass ~ Subplot_Descriptions*as.factor(Year), data=anpp1)# testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # still not normal but plot looks a lot better
# better - going with log transformation

m5 <- lmer(log_biomass ~ Subplot_Descriptions*as.factor(Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
anova(m5)
summary(m5)
contrast(emmeans(m5, ~ Subplot_Descriptions*as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")

# using this model 
m6 <- lmer(log_biomass ~ Subplot_Descriptions + as.factor(Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
anova(m6)
summary(m6)
contrast(emmeans(m6, ~as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")
# back-transforming - 2022 vs 2023
exp(1.1680)-exp(4.43673-0.22449) 
exp(1.1680)-exp(4.43673-0.22449) 

# Extract fixed effects
fixed_effects <- fixef(m6)
# Extract random effects
random_effects <- ranef(m6)$`Rep %in% footprint_unique`$`subplot_unique`
# Calculate squared random effects
squared_random_effects <- random_effects^2
# Calculate squared fixed effects
squared_fixed_effects <- fixed_effects^2

# Calculate sum of squared fixed effects
sum_squared_fixed_effects <- sum(squared_fixed_effects)

# Calculate sum of squared random effects
sum_squared_random_effects <- sum(squared_random_effects)

# Calculate residual variance
residual_variance <- sigma(m6)^2

# Calculate partial eta-squared
partial_eta_squared <- sum_squared_fixed_effects / (sum_squared_fixed_effects + sum_squared_random_effects + residual_variance)


m7 <- lmer(log_biomass ~ Subplot_Descriptions + (1|Year) + (1|Rep/footprint_unique/subplot_unique), data = anpp1, REML=F)
anova(m7)
summary(m7)

################################################################################
# 2023 DATA

# filter for 2023 data and 1.0m2 data
anpp23 <- anpp %>%
        filter(Year=='2023') %>%
        filter(Scale_meter_square==1)

anpp23_2 <- anpp %>%
        filter(Year=='2023') %>%
        filter(!(Scale_meter_square==1))

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
qqPlot(resid(m23_full), main = "Biomass")
hist(residuals(m23_full), main = "Biomass")
shapiro.test(resid(m23_full))
outlierTest(m23_full) # checking for outliers - none

###### Checking model results ########
anova(m23_full)
summary(m23_full)

################################################################################
# Using all years data with the 2021 and 2022 data scaled to 1m2 (done at the beginning and multiplied 0.2 
# values by 5)

###### Data exploration #######
descdist(anpp_scaled$plant_biomass_gm2, discrete = FALSE) # checks what type of distribution your data matches
hist(anpp_scaled$plant_biomass_gm2)
qqnorm(anpp_scaled$plant_biomass_gm2) # should be a straight diagonal line if normal
shapiro.test(anpp_scaled$plant_biomass_gm2) # p > 0.05 if normal
# making a model to test distribution of model residuals
raw.data.test <- lm(plant_biomass_gm2 ~ Subplot_Descriptions*Year, data=anpp_scaled) # testing model
hist(resid(raw.data.test)) # checking model residuals
shapiro.test(resid(raw.data.test))

# log transformation
anpp_scaled$log_biomass <- log(anpp_scaled$plant_biomass_gm2)
anpp_scaled <- anpp_scaled[is.finite(anpp_scaled$log_biomass), ] # Remove rows with Inf in column 'x'
descdist(anpp_scaled$log_biomass, discrete = FALSE)
hist(anpp_scaled$log_biomass)
qqnorm(anpp_scaled$log_biomass) # a little weird
shapiro.test(anpp_scaled$log_biomass) # still not normal
# making a model to test distribution of model residuals w/ log transformation
raw.data.test.log <- lm(log_biomass ~ Subplot_Descriptions*Year, data=anpp_scaled)# testing model
hist(resid(raw.data.test.log)) # checking model residuals
shapiro.test(resid(raw.data.test.log)) # still not normal but plot looks a lot better
# better - going with log transformation

######### Assumption checking ##########
# in the assumption checking, we're making sure that our full model meets the assumptions of the model.
# the model below is the model structure we can use for all response variables; its testing to see
# if there is 1. an effect of climate treatment on height, 2. an effect of galling status on height, and 3. does the effect
# of climate on height depend on galling status. Subplot nested within footprint nested within rep is used as our random effect
# to account for variation between plots. Year is also included as a random effect to account for variation between years.
m8 <- lmer(log_biomass ~ Subplot_Descriptions*Year + (1|Rep/footprint_unique/subplot_unique), data = anpp_scaled, REML=F)
# Check Assumptions:
# (1) Linearity: if covariates are not categorical
# (2) Homogeneity: Need to Check by plotting residuals vs predicted values.
plot(m8, main = "Biomass")
# Homogeneity of variance looks a bit off (increasing variance in resids does increase with fitted values)
# Check for homogeneity of variances (true if p>0.05). If the result is not significant, the assumption of equal variances (homoscedasticity) is met (no significant difference between the group variances).
leveneTest(residuals(m8) ~ anpp_scaled$Subplot_Descriptions)
# (3) Normality of error term: need to check by histogram, QQplot of residuals, could do Kolmogorov-Smirnov test.
# Check for normal residuals
qqPlot(resid(m8), main = "Biomass")
hist(residuals(m8), main = "Biomass")
shapiro.test(resid(m8))
outlierTest(m8) # checking for outliers - none

anova(m8)
summary(m8)

m9 <- lmer(log_biomass ~ Subplot_Descriptions + Year + (1|Rep/footprint_unique/subplot_unique), data = anpp_scaled, REML=F)
anova(m9)
summary(m9)

contrast(emmeans(m1, ~Subplot_Descriptions + as.factor(Year)), "pairwise", simple = "each", combine = F, adjust = "mvt")
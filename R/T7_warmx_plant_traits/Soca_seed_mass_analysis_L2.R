# TITLE:          REX: Soca seed mass analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Emily Parker, Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Feb 2024


# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(pscl)
library(performance)
library(fitdistrplus)
library(MASS)
library(lmtest)
library(glmmTMB)
library(emmeans)

# get data
dir<-Sys.getenv("DATA_DIR")
list.files(dir)
seed <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_soca_seeds_mass_L1.csv"))
biomass <- read.csv(file.path(dir, "/T7_warmx_plant_traits/L1/T7_warmx_soca_biomass_L1.csv"))
str(seed)
# storing year as a factor
seed$Year <- as.factor(seed$Year)
biomass$Year <- as.factor(biomass$Year)
# adding unique plant ID for each subplot
seed <- seed %>%
  group_by(Treatment, Rep, Footprint, Subplot) %>%
  mutate(plant_num = 1:n())

# merge seed and biomass dataframes
# note to self: 267 is duplicated in 2021 - fix this
seed <- left_join(seed,biomass,by=c("Treatment","Rep","Footprint","Subplot","Climate_Treatment","Galling_Status","Year","Unique_ID"))
#removing rep 4
seed <- seed %>%
  filter(!(Rep == 4 & Climate_Treatment == "Warm Drought")) %>%
  filter(!(Rep == 4 & Climate_Treatment == "Ambient Drought"))


##################################### distribution check #####################################
# first checking basic histogram
hist(seed$Seeds_Mass) # very right-skewed

# how much of the data is zeros?
100*sum(seed$Seeds_Mass == 0)/nrow(seed) # 35.8%

### determining distribution ###
descdist(seed$Seeds_Mass, discrete = FALSE)
# log transform?
seed$seed_log <- log(seed$Seeds_Mass+1)
hist(seed$seed_log)
fit <- lm(seed_log~1, data = seed)
hist(resid(fit))
# mean centering?
seed$seed_scaled <- seed$seed_log - mean(seed$seed_log)
hist(seed$seed_scaled)
fit_scaled <- lm(seed_scaled~1, data = seed)
hist(resid(fit_scaled))
# square root?
seed$seed_sqrt <- sqrt(seed$Seeds_Mass)
hist(seed$seed_sqrt)
fit_sqrt <- lm(seed_sqrt~1, data = seed)
hist(resid(fit_sqrt))
shapiro.test(seed$seed_sqrt)

# transformations are a no-go
# gamma distribution?
fit.gamma <- fitdist(seed$Seeds_Mass, distr = "gamma", method = "mme")
plot(fit.gamma)
# gamma seems to fit well - we can do a zero-inflated gamma model using glmmTMB



##### full model #####
seed <- within(seed, Climate_Treatment <- relevel(factor(Climate_Treatment), ref = "Ambient Drought"))
seed <- within(seed, Galling_Status <- relevel(factor(Galling_Status), ref = "Non-Galled"))
full.model <- glmmTMB(Seeds_Mass ~ Climate_Treatment * Galling_Status + (1|Year:Rep:Footprint) + (1|Biomass),
                      data=seed,
                      family=ziGamma(link="log"),
                      zi=~.)
summary(full.model)
car::Anova(full.model) # slight significance of interaction term, going to check all pairwise comparisons
# pairwise comparisons #
# note: in zero-inflated model, we're testing the probability of being zero
# negative estimate (<0) means fewer 0's, positive estimate (>0) means more 0's
# positive estimate means it is more likely to be zero (or, there is a reduced probability of of having a seed)
# negative estimates means it is less likely to be zero (or, there is an increased probability of having a seed)
emmeans(full.model, ~Climate_Treatment*Galling_Status, type = "response") # seed weight
emmeans(full.model, ~Climate_Treatment*Galling_Status, component = "zi", type = "response") # probability

contrast(emmeans(full.model, ~Climate_Treatment*Galling_Status), "pairwise", simple="each", combine = F, adjust="mvt")
contrast(emmeans(full.model, ~Climate_Treatment*Galling_Status, component = "zi",), "pairwise", simple="each", combine = F, adjust="mvt")


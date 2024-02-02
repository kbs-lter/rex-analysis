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

# get data
dir<-Sys.getenv("DATA_DIR")
list.files(dir)
seed <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_soca_seeds_mass_L1.csv"))
str(seed)
# storing year as a factor
seed$Year <- as.factor(seed$Year)
# adding unique plant ID for each subplot
seed <- seed %>%
  group_by(Treatment, Rep, Footprint, Subplot) %>%
  mutate(plant_num = 1:n())



##################################### distribution check #####################################
# first checking basic histogram
hist(seed$Seeds_Mass) # very right-skewed

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

# how much of the data is zeros?
100*sum(seed$Seeds_Mass == 0)/nrow(seed) # 35.8%

# checking mean and variance of non-zero counts
seed %>%
  dplyr::filter(Seeds_Mass != "0") %>%
  dplyr::summarize(mean_seed = mean(Seeds_Mass, na.rm=T), var_seed = var(Seeds_Mass, na.rm=T))
# mean is > than variance. poisson distributions assume mean = variance.
# negative binomial might fit better

# comparing poisson vs. neg bin
pois <- glm(Seeds_Mass ~ Climate_Treatment + Galling_Status, family = poisson(), data = seed)
nb <- glm.nb(Seeds_Mass ~ Climate_Treatment + Galling_Status, data = seed)
lrtest(pois,nb) # negative binomial outperforms poisson


# checking zero-inflation
m.test <- glm(Seeds_Mass ~ Climate_Treatment + Galling_Status, family = poisson, data = seed)
check_zeroinflation(m.test)



##### full model #####
full.model <- glmmTMB(Seeds_Mass ~ Climate_Treatment * Galling_Status + (1|Year/Treatment/Rep/Footprint/Subplot/plant_num),
                      data=seed,
                      family=ziGamma(link="log"),
                      zi=~.)
summary(full.model) # * used in paper * #







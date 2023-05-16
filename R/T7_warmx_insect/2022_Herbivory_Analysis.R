# TITLE:          2022 Herbivory Analysis
# AUTHORS:        Emily Parker
# COLLABORATORS:  Phoebe Zarnetske, Kara Dobson, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv https://drive.google.com/drive/folders/1mFnCu6bHhdMzkecyh6PDWCI6WSObz3or
# DATA OUTPUT:    plots, stats
# PROJECT:        REX
# DATE:           2 May 2023
# NOTES:          

#  Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(plotrix)
library(lmerTest)
library(fitdistrplus)
library(emmeans)

#set directory
dir <- setwd("/Users/emilyparker/Documents/R")
dir<-Sys.getenv("DATA_DIR") # Kara's

#read in data
herbvar <- read.csv(file.path(dir,"T7_REX_goldenrod_herbvar_2022.csv"))
zarnetske <- read.csv(file.path(dir,"T7_REX_goldenrod_herbivory_2022.csv"))

herbvar <- read.csv(file.path(dir,"T7_warmx_insect/L0/T7_REX_goldenrod_herbvar_2022.csv")) # Kara's
zarnetske <- read.csv(file.path(dir,"T7_warmx_insect/L0/T7_REX_goldenrod_herbivory_2022.csv")) # Kara's


# merge w/ meta data
# this should be done in a cleaning script, but doing it here for now
meta_warmx <- read.csv(file.path(dir, "REX_warmx_metadata.csv"))
# renaming meta data columns to match data & fixing lowercase subplot names to match
colnames(meta_warmx)[colnames(meta_warmx) == "Rep"] ="rep"
colnames(meta_warmx)[colnames(meta_warmx) == "Footprint_Location"] ="footprint"
colnames(meta_warmx)[colnames(meta_warmx) == "Subplot_Location"] ="quad"
meta_warmx$quad = toupper(meta_warmx$quad)
# merge
zarnetske <- left_join(zarnetske,meta_warmx,by=c("rep","footprint","quad"))
herbvar <- left_join(herbvar,meta_warmx,by=c("rep","footprint","quad"))
# remove insecticide plots from herbvar (for now) since the zarnetske data does not have them
herbvar <- herbvar[- grep("insecticide", herbvar$Subplot_Descriptions),]



### correlation btwn previous (zarnetske) method & new herbvar method ###
#define variables
plant_num <- 1:294

#zarnetske - average of leaves sampled
zarnetske_leaf_avg <- zarnetske %>%
  group_by(plant_num) %>%
  summarize(zarnetske_leaf_avg = mean(herbiv_percent, na.rm = TRUE),
            zarnetske_leaf_se = std.error(herbiv_percent, na.rm = TRUE))

#herbvar - average of leaves sampled
herbvar_leaf_avg <- herbvar %>%
  group_by(plant_num) %>%
  summarize(herbvar_leaf_avg = mean(leaf_avg, na.rm = TRUE),
            herbvar_leaf_se = std.error(leaf_avg, na.rm = TRUE))

#herbvar - average of total plant herbivory
herbvar_plant_avg <- herbvar %>%
  group_by(plant_num) %>%
  summarize(herbvar_plant_avg = mean(total_plant_herb, na.rm = TRUE),
            herbvar_plant_se = std.error(total_plant_herb, na.rm = TRUE))

## would this be helpful? not sure how to do it - would need help
# herbvar - estimated herbivory vs theoretical (from leaf average and number of leaves)
# theoretical = (leaves_herb * leaf_avg) / total_leaves

#create dataframe w/ only plants with herbvar and zarnetske avgs
herbivory <- data.frame(plant_num)

herbivory <- zarnetske_leaf_avg %>% inner_join(herbivory, by =c('plant_num' = 'plant_num'))
herbivory <- herbvar_leaf_avg %>% inner_join(herbivory, by =c('plant_num' = 'plant_num'))
herbivory <- herbvar_plant_avg %>% inner_join(herbivory, by =c('plant_num' = 'plant_num'))

# regression - zarnetske vs herbvar leaf avg
cor.test(herbivory$herbvar_leaf_avg,herbivory$zarnetske_leaf_avg, method = "pearson")
lm1 <- lm(herbvar_leaf_avg ~ zarnetske_leaf_avg, data = herbivory)
plot(herbvar_leaf_avg ~ zarnetske_leaf_avg, data = herbivory)
abline(lm1)
summary(lm1)

# regression - herbvar total plant avg vs leaves avg
cor.test(herbivory$herbvar_leaf_avg,herbivory$herbvar_plant_avg, method = "pearson")
lm2 <- lm(herbvar_leaf_avg ~ herbvar_plant_avg, data = herbivory)
plot(herbvar_leaf_avg ~ herbvar_plant_avg, data = herbivory)
abline(lm2)
summary(lm2)

# regression - herbvar total plant avg vs zarnetske
cor.test(herbivory$herbvar_plant_avg,herbivory$zarnetske_leaf_avg, method = "pearson")
lm3 <- lm(herbvar_plant_avg ~ zarnetske_leaf_avg, data = herbivory)
plot(herbvar_plant_avg ~ zarnetske_leaf_avg, data = herbivory)
abline(lm3)
summary(lm3)



#### herbivory distribution check ####
# previous herbivory data had an excess of zeros (>50% of data as zeros), so here I'm checking
# what percent of the REX herb. data is zeros for the two herbivory methods
100*sum(zarnetske$herbiv_percent == 0)/nrow(zarnetske)
100*sum(herbvar$leaf_avg == 0)/nrow(herbvar)
# low # of 0's, so a zero-inflated model is not needed

# determining distribution - zarnetske data
descdist(zarnetske$herbiv_percent, discrete = FALSE)
hist(zarnetske$herbiv_percent)
qqnorm(zarnetske$herbiv_percent)
shapiro.test(zarnetske$herbiv_percent)
# right skewed, trying transformations
# square root
zarnetske$herbiv_percent_sqrt <- sqrt(zarnetske$herbiv_percent)
descdist(zarnetske$herbiv_percent_sqrt, discrete = FALSE)
hist(zarnetske$herbiv_percent_sqrt)
qqnorm(zarnetske$herbiv_percent_sqrt)
shapiro.test(zarnetske$herbiv_percent_sqrt)

# determining distribution - herbvar data
herbvar_na_rem <- herbvar %>% # remove rows w/ NA in leaf_avg column
  drop_na(leaf_avg)
descdist(herbvar_na_rem$leaf_avg, discrete = FALSE)
hist(herbvar_na_rem$leaf_avg)
qqnorm(herbvar_na_rem$leaf_avg)
shapiro.test(herbvar_na_rem$leaf_avg)
# right skewed, trying transformations
# log transformation
herbvar_na_rem$leaf_avg_log <- log(herbvar_na_rem$leaf_avg+1)
descdist(herbvar_na_rem$leaf_avg_log, discrete = FALSE)
hist(herbvar_na_rem$leaf_avg_log)
qqnorm(herbvar_na_rem$leaf_avg_log)
shapiro.test(herbvar_na_rem$leaf_avg_log)



### herbivory models ###
# zarnetske data
zarn_m1 <- lmer(herbiv_percent_sqrt ~ Subplot_Descriptions + (1|rep/footprint/plant_num), zarnetske, REML=FALSE)
anova(zarn_m1)
# comparisons
emmeans(zarn_m1, list(pairwise ~ Subplot_Descriptions), adjust = "tukey")
zarnm1_emm <- emmeans(zarn_m1, ~ Subplot_Descriptions)
contrast(zarnm1_emm, "pairwise", simple = "each", combine = F, adjust = "mvt")

# herbvar data
herbvar_m1 <- lmer(leaf_avg_log ~ Subplot_Descriptions + (1|rep/footprint), herbvar_na_rem, REML=FALSE)
anova(herbvar_m1)
# comparisons
emmeans(herbvar_m1, list(pairwise ~ Subplot_Descriptions), adjust = "tukey")
herbvarm1_emm <- emmeans(herbvar_m1, ~ Subplot_Descriptions)
contrast(herbvarm1_emm, "pairwise", simple = "each", combine = F, adjust = "mvt")


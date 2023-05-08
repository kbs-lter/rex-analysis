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

#set directory
dir <- setwd("/Users/emilyparker/Documents/R")

#read in data
herbvar <- read.csv(file.path(dir,"T7_REX_goldenrod_herbvar_2022.csv"))
zarnetske <- read.csv(file.path(dir,"T7_REX_goldenrod_herbivory_2022.csv"))


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


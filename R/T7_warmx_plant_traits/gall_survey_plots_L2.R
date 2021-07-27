# TITLE:          REX: Gall density plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Figures
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

# removing outlier from irrigated control before plotting
gall <- gall[!grepl(40, gall$total_rosette_quadrat),]

# Figure for gall counts between each treatment
png("gall_count.png", units="in", width=6, height=5, res=300)
ggplot(gall, aes(x = treatment, y = total_rosette_quadrat)) +
  geom_boxplot(color = "black", fill = "olivedrab") +
  labs(x = "Treatment", y = "Total # of Galls") +
  scale_x_discrete(limits = c("irrigated_control", "ambient_control", "drought", "warmed", "warmed_drought"),
                   labels=c("ambient_control" = "Control",
                            "drought" = "Drought",
                            "irrigated_control" = "Irrigated Control",
                            "warmed" = "Warmed",
                            "warmed_drought" = "Warmed & Drought"),
                   guide = guide_axis(n.dodge=2)) +
  #geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6) +
  theme_classic()
dev.off()

# getting proportion of galled stems in each treatment
gall <- gall %>%
  mutate(prop_gall = rosette_galls/total_stems)

# Figure for proportion of galled stems between each treatment
png("gall_prop.png", units="in", width=6, height=5, res=300)
ggplot(gall, aes(x = treatment, y = prop_gall)) +
  geom_boxplot(color = "black", fill = "olivedrab") +
  labs(x = "Treatment", y = "Proportion of Galled Stems") +
  scale_x_discrete(limits = c("irrigated_control", "ambient_control", "drought", "warmed", "warmed_drought"),
                   labels=c("ambient_control" = "Control",
                            "drought" = "Drought",
                            "irrigated_control" = "Irrigated Control",
                            "warmed" = "Warmed",
                            "warmed_drought" = "Warmed & Drought"),
                   guide = guide_axis(n.dodge=2)) +
  #geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6) +
  theme_classic()
dev.off()

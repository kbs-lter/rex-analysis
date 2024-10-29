# TITLE:          REX: Gall biomass plots
# AUTHORS:        Emily Parker
# COLLABORATORS:  Kara Dobson, Phoebe Zarnetske, Moriah Young, Kristin Wolford, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Figures
# PROJECT:        REX
# DATE:           October 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(plotrix)

# Set ggplot2 plotting
# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(14))
theme_update(axis.text.x = element_text(size = 12),
             axis.text.y = element_text(size = 16),
             axis.title = element_text(size=16,face="bold"))

# Set working directory
dir <- setwd("C:/Users/Emily/Documents/R/Goldenrod Project")

# Read in data
galls <- read.csv(file.path(dir, "L1/T7_warmx_soca_galls_L1.csv"))

# Take subplot average dried of galls
gall_avg <- galls %>%
  group_by(Climate_Treatment) %>%
  summarize(avg_gall = mean(Dried_Weight, na.rm = TRUE),
            se = std.error(Dried_Weight, na.rm = TRUE))

#pointrange
## color isn't adding correctly, but rest is correct
png("gall_mass_overall_point.png", units="in", width=9, height=6, res=300)
ggplot(gall_avg, aes(x = Climate_Treatment, y = avg_gall)) +
  geom_pointrange(aes(ymin = avg_gall - se, ymax = avg_gall + se),pch=21,size=1,position=position_dodge(0.2), fill = "purple4") +
  labs(x = NULL, y = "Gall dried mass (g)", title=NULL) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient", "Warm" = "Warmed",
                            "Ambient Drought" = "Drought", "Irrigated Control" = "Irrigated\nControl",
                            "Warm Drought" = "Warmed\nDrought"))

dev.off()


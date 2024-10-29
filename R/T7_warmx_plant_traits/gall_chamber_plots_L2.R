# TITLE:          REX: gall chamber plots
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
chmbr <- read.csv(file.path(dir, "L1/T7_warmx_soca_chmbr_chmb_vol_L1.csv"))
count <- read.csv(file.path(dir, "L1/T7_warmx_Soca_chmbr_chmb_count_L1.csv"))

# Take subplot average vol of chmbr
chmbr_avg <- chmbr %>%
  group_by(treatment) %>%
  summarize(avg_chmbr = mean(chamber_volume_mm3, na.rm = TRUE),
            se = std.error(chamber_volume_mm3, na.rm = TRUE))

#avg count
count_avg <- count %>%
  group_by(treatment) %>%
  summarize(avg_count = mean(num_of_chambers, na.rm = TRUE),
            se = std.error(num_of_chambers, na.rm = TRUE))

#pointrange - vol
png("chmbr_vol_overall_point.png", units="in", width=9, height=6, res=300)
ggplot(chmbr_avg, aes(x = treatment, y = avg_chmbr)) +
  geom_pointrange(aes(ymin = avg_chmbr - se, ymax = avg_chmbr + se),pch=21,size=1,position=position_dodge(0.2), fill = "purple4") +
  labs(x = NULL, y = "Chamber Volume (mm^3)", title=NULL) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient", "Warm" = "Warmed",
                            "Ambient Drought" = "Drought", "Irrigated Control" = "Irrigated\nControl",
                            "Warm Drought" = "Warmed\nDrought"))

dev.off()

#pointrange - count
png("chmbr_count_overall_point.png", units="in", width=9, height=6, res=300)
ggplot(count_avg, aes(x = treatment, y = avg_count)) +
  geom_pointrange(aes(ymin = avg_count - se, ymax = avg_count + se),pch=21,size=1,position=position_dodge(0.2), fill = "purple4") +
  labs(x = NULL, y = "Number of Larval Chambers", title=NULL) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient", "Warm" = "Warmed",
                            "Ambient Drought" = "Drought", "Irrigated Control" = "Irrigated\nControl",
                            "Warm Drought" = "Warmed\nDrought"))

dev.off()


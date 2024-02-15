# TITLE:          REX: Soca plant height plots for harvest time point
# AUTHORS:        Emily Parker, Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Feb 2024

# Clear all existing data
rm(list=ls())

# Load packages
# install.packages("tidyverse") - Emily
# install.packages("plotrix") - Emily
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
dir <- setwd("/Users/emilyparker/Documents/R/Goldenrod Project 2022")

# Read in data
height <- read.csv(file.path(dir, "L1/T7_warmx_soca_height_harvest_L1.csv"))

# Take averages
height_avg <- height %>%
  group_by(Climate_Treatment) %>%
  summarize(avg_height = mean(Height_cm, na.rm = TRUE),
            se = std.error(Height_cm, na.rm = TRUE))
height_avg2 <- height %>%
  group_by(Galling_Status) %>%
  summarize(avg_height = mean(Height_cm, na.rm = TRUE),
            se = std.error(Height_cm, na.rm = TRUE))

height_avg3 <- height %>%
  group_by(Climate_Treatment, Galling_Status) %>%
  summarize(avg_height = mean(Height_cm, na.rm = TRUE),
            se = std.error(Height_cm, na.rm = TRUE))


# Barplot of temp treatment
png("height_harv_temp_bar.png", units="in", width=10, height=6, res=300)
ggplot(height_avg, aes(x = Climate_Treatment, y = avg_height)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Height (g)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=2)) +
  theme(legend.position = "none")
dev.off()


#barplot of galling X mass
png("height_harv_galling_bar.png", units="in", width=10, height=6, res=300)
ggplot(height_avg2, aes(x = Galling_Status, y = avg_height)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = "identity") +
  labs(x = "Galling Status", y = "Height (g)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Non-Galled", "Galled"),
                   labels=c("Galled" = "Galled",
                            "Non-Galled" = "Non-Galled"),
                   guide = guide_axis(n.dodge=2)) +
  theme(legend.position = "none")
dev.off()

# barplot of treatment X height X galling
png("height_harv_overall_bar.png", units="in", width=9, height=6, res=300)
ggplot(height_avg3, aes(x = Climate_Treatment, y = avg_height, fill=Galling_Status)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = position_dodge(0.85)) +
  labs(x = "Treatment", y = "Height (g)") +
  scale_fill_manual(values = c("olivedrab4", "darkseagreen2"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "right")
dev.off()


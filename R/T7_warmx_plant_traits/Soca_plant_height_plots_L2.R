# TITLE:          REX: Soca plant height plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Oct 2023

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
# Emily doesn't need this line
dir<-Sys.getenv("DATA_DIR")

# Read in data
height <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_soca_height_L1.csv"))
height2 <- height[-c(437,225,422,449),] # removing same outliers from analyses
# height <- read.csv("downloads/T7_warmx_Soca_plant_height_L1.csv") - Emily

# Take averages
height_avg <- height2 %>%
  group_by(Climate_Treatment, Galling_Status) %>%
  summarize(avg_height = mean(Height_cm, na.rm = TRUE),
            se = std.error(Height_cm, na.rm = TRUE))

# plot of means for climate treatment x galling status
png("height_galling_point.png", units="in", width=7, height=5, res=300)
ggplot(height_avg, aes(x = Climate_Treatment, y = avg_height, fill=Galling_Status)) +
  geom_pointrange(aes(ymin = avg_height - se, ymax = avg_height + se),pch=21,size=1,position=position_dodge(0.2)) +
  labs(x = NULL, y = "Plant height (cm)") +
  scale_fill_manual(values = c("purple4", "plum2"),
                    name="Galling status",
                    labels=c("Galled","Non-galled")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme_bw(14) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"))
dev.off()

# barplot of climate treatment + galling status
png("height_galling_bar.png", units="in", width=9, height=6, res=300)
ggplot(height_avg, aes(x = Climate_Treatment, y = avg_height, fill=Galling_Status)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = position_dodge(0.85)) +
  labs(x = "Treatment", y = "Plant Height (cm)") +
  scale_fill_manual(values = c("olivedrab4", "darkseagreen2"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "right")
dev.off()

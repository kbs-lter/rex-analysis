# TITLE:          REX: Soca infl mass plots
# AUTHORS:        Kara Dobson, Emily Parker
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Dec 2023

# Clear all existing data
rm(list=ls())

# Load packages
#install.packages("tidyverse") #Emily
#install.packages("plotrix") #Emily
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
mass <- read.csv(file.path(dir, "L1/T7_warmx_soca_infl_mass_L1.csv"))


# Take averages
seeds_avg <- mass %>%
  group_by(Climate_Treatment) %>%
  summarize(avg_mass = mean(Seeds_Mass, na.rm = TRUE),
            se = std.error(Seeds_Mass, na.rm = TRUE))
seeds_avg2 <- mass %>%
 group_by(Galling_Status) %>%
  summarize(avg_mass = mean(Seeds_Mass, na.rm = TRUE),
         se = std.error(Seeds_Mass, na.rm = TRUE))

seeds_avg3 <- mass %>%
  group_by(Climate_Treatment, Galling_Status) %>%
  summarize(avg_mass = mean(Seeds_Mass, na.rm = TRUE),
            se = std.error(Seeds_Mass, na.rm = TRUE))


# Barplot of temp treatment
png("mass_temp_bar.png", units="in", width=10, height=6, res=300)
ggplot(seeds_avg, aes(x = Climate_Treatment, y = avg_mass)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_mass - se, ymax = avg_mass + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Seeds Mass (g)") +
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
png("mass_galling_bar.png", units="in", width=10, height=6, res=300)
ggplot(seeds_avg2, aes(x = Galling_Status, y = avg_mass)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_mass - se, ymax = avg_mass + se), width = 0.2,
                position = "identity") +
  labs(x = "Galling Status", y = "Seeds Mass (g)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Non-Galled", "Galled"),
                   labels=c("Galled" = "Galled",
                            "Non-Galled" = "Non-Galled"),
                   guide = guide_axis(n.dodge=2)) +
  theme(legend.position = "none")
dev.off()

# barplot of treatment X height X galling
png("mass_overall_bar.png", units="in", width=9, height=6, res=300)
ggplot(seeds_avg3, aes(x = Climate_Treatment, y = avg_mass, fill=Galling_Status)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = avg_mass - se, ymax = avg_mass + se), width = 0.2,
                position = position_dodge(0.85)) +
  labs(x = "Treatment", y = "Seeds Mass (cm)") +
  scale_fill_manual(values = c("olivedrab4", "darkseagreen2"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "right")
dev.off()

# TITLE:          REX: Greenness plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           July 2021

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(plotrix)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
green <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_greenness_L1.csv"))

# Set ggplot2 plotting
# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(14))
theme_update(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             axis.title = element_text(size=16,face="bold"))

# Boxplot
png("rex_green.png", units="in", width=8, height=6, res=300)
ggplot(green, aes(x = treatment, y = greenness, fill = gall_present)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  labs(x = "Treatment", y = "Leaf Greenness", fill = "Gall Presence") +
  scale_fill_manual(values = c("olivedrab", "olivedrab1"), labels = c("Gall", "No Gall")) +
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=16)) +
  scale_x_discrete(limits = c("irr_control", "ambient", "drought", "warmed", "warmed_drought"),
                   labels=c("ambient" = "Ambient",
                            "drought" = "Drought",
                            "irr_control" = "Irrigated \n Control",
                            "warmed" = "Warmed",
                            "warmed_drought" = "Warmed & \n Drought"))
                   #guide = guide_axis(n.dodge=2))
  #geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = gall_present)) +
  #theme_classic()
dev.off()

# Gall average plot
png("gall_green.png", units="in", width=6, height=5, res=300)
ggplot(green, aes(x = gall_present, y = greenness, fill = gall_present)) +
  geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = gall_present)) +
  scale_color_manual(values = c("gall" = "olivedrab", "no_gall" = "olivedrab")) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  labs(x = "Gall Presence", y = "Leaf Greenness") +
  scale_fill_manual(values = c("olivedrab", "olivedrab")) +
  scale_x_discrete(labels=c("gall" = "Gall",
                            "no_gall" = "No Gall")) +
  #theme_classic() +
  theme(legend.position = "none")
dev.off()

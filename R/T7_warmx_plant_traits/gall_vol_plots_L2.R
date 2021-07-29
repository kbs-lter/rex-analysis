# TITLE:          REX: Gall volume plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Figures
# PROJECT:        REX
# DATE:           July 2021

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
galls <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_gall_vol_L1.csv"))

# Take subplot average of galls
gall_avg <- galls %>%
  group_by(rep, footprint, treatment) %>%
  summarize(avg_gall = mean(sphere_vol, na.rm = TRUE))

# Set ggplot2 plotting
# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(14))
theme_update(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             axis.title = element_text(size=16,face="bold"))

# Boxplot of volume x treatment
png("vol_treatment.png", units="in", width=8, height=6, res=300)
ggplot(gall_avg, aes(x = treatment, y = avg_gall)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Treatment", y = "Gall Volume (cm^3)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("irr_control", "ambient", "drought", "warmed", "warmed_drought"),
                   labels=c("ambient" = "Ambient",
                            "drought" = "Drought",
                            "irr_control" = "Irrigated \n Control",
                            "warmed" = "Warmed",
                            "warmed_drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=2)) +
  #geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = gall_present)) +
  theme(legend.position = "none")
dev.off()

# Regression
lm_gall <- lm(sphere_vol ~ plant_height, data = galls)
png("vol_height.png", units="in", width=5, height=5, res=300)
ggplot(galls,aes(plant_height, sphere_vol)) +
  geom_point(color='olivedrab') + 
  geom_smooth(method='lm', color="olivedrab") +
  labs(x = "Plant Height (cm)", y = "Gall Volume (cm^3)")
dev.off()
summary(lm_gall)

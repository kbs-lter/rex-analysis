# TITLE:          REX: Gall volume plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Figures
# PROJECT:        REX
# DATE:           July 2021; updated Jan 2022

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
dir<-Sys.getenv("DATA_DIR")

# Read in data
galls <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_Soca_gall_vol_L1.csv"))

# Take subplot average of galls
gall_avg <- galls %>%
  group_by(treatment, drought_period) %>%
  summarize(avg_gall = mean(sphere_vol_cm3, na.rm = TRUE),
            se = std.error(sphere_vol_cm3, na.rm = TRUE))
gall_avg2 <- galls %>%
  group_by(treatment) %>%
  summarize(avg_gall = mean(sphere_vol_cm3, na.rm = TRUE),
            se = std.error(sphere_vol_cm3, na.rm = TRUE))
gall_avg3 <- galls %>%
  group_by(drought_period) %>%
  summarize(avg_gall = mean(sphere_vol_cm3, na.rm = TRUE),
            se = std.error(sphere_vol_cm3, na.rm = TRUE))

# Boxplot of volume x treatment
gall_avg <- gall_avg %>%
  mutate(across(drought_period, factor, levels=c("Pre-Drought","Drought","Post-Drought"))) # re-ordering drought period
png("gall_volume_overall.png", units="in", width=10, height=6, res=300)
ggplot(galls, aes(x = treatment, y = sphere_vol_cm3)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  facet_grid(.~drought_period) +
  labs(x = "Treatment", y = "Gall Volume (cm^3)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=2)) +
  ylim(0,100) +
  geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Boxplot of volume w/o drought period
png("gall_volume_temp.png", units="in", width=8, height=6, res=300)
ggplot(galls, aes(x = treatment, y = sphere_vol_cm3)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Treatment", y = "Gall Volume (cm^3)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=2)) +
  ylim(0,100) +
  geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Boxplot of volume only across drought periods
png("gall_volume_time.png", units="in", width=8, height=6, res=300)
ggplot(galls, aes(x = drought_period, y = sphere_vol_cm3)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Timing of Drought Period", y = "Gall Volume (cm^3)") +
  #scale_fill_manual(values = c("olivedrab")) +
  #scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
  #                 labels=c("Ambient" = "Ambient",
  #                          "Ambient Drought" = "Drought",
  #                          "Irrigated Control" = "Irrigated \n Control",
  #                          "Warm" = "Warmed",
  #                          "Warm Drought" = "Warmed & \n Drought"),
  #                 guide = guide_axis(n.dodge=2)) +
  ylim(0,100) +
  geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Barplot of volume x treatment
png("gall_volume_overall_bar.png", units="in", width=10, height=6, res=300)
ggplot(gall_avg, aes(x = treatment, y = avg_gall)) +
  facet_grid(.~drought_period) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_gall - se, ymax = avg_gall + se), width = 0.2,
                         position = "identity") +
  labs(x = "Treatment", y = "Gall Volume (cm^3)") +
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

# Barplot of only temperature treatment
png("gall_volume_temp_bar.png", units="in", width=8, height=6, res=300)
ggplot(gall_avg2, aes(x = treatment, y = avg_gall)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_gall - se, ymax = avg_gall + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Gall Volume (cm^3)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "none")
dev.off()

# Barplot of only drought period
png("gall_volume_time_bar.png", units="in", width=8, height=6, res=300)
ggplot(gall_avg3, aes(x = drought_period, y = avg_gall)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_gall - se, ymax = avg_gall + se), width = 0.2,
                position = "identity") +
  labs(x = "Timing of Drought Period", y = "Gall Volume (cm^3)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  #scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
  #                 labels=c("Ambient" = "Ambient",
  #                          "Ambient Drought" = "Drought",
  #                          "Irrigated Control" = "Irrigated \n Control",
  #                          "Warm" = "Warmed",
  #                          "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "none")
dev.off()

## Regression - data does not contain height anymore
#lm_gall <- lm(sphere_vol_cm3 ~ plant_height, data = galls)
#png("vol_height.png", units="in", width=5, height=5, res=300)
#ggplot(galls,aes(plant_height, sphere_vol)) +
#  geom_point(color='olivedrab') + 
#  geom_smooth(method='lm', color="olivedrab") +
#  labs(x = "Plant Height (cm)", y = "Gall Volume (cm^3)")
#dev.off()
#summary(lm_gall)

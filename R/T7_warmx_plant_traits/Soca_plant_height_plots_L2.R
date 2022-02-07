# TITLE:          REX: Soca plant height plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_plant_traits L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Jan 2022

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
height <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_Soca_plant_height_L1.csv"))

# Take averages
height_avg <- height %>%
  group_by(treatment, drought_period) %>%
  summarize(avg_height = mean(plant_height_cm, na.rm = TRUE),
            se = std.error(plant_height_cm, na.rm = TRUE))
height_avg2 <- height %>%
  group_by(treatment) %>%
  summarize(avg_height = mean(plant_height_cm, na.rm = TRUE),
            se = std.error(plant_height_cm, na.rm = TRUE))
height_avg3 <- height %>%
  group_by(drought_period) %>%
  summarize(avg_height = mean(plant_height_cm, na.rm = TRUE),
            se = std.error(plant_height_cm, na.rm = TRUE))
height_avg4 <- height %>%
  group_by(treatment, gall_present) %>%
  summarize(avg_height = mean(plant_height_cm, na.rm = TRUE),
            se = std.error(plant_height_cm, na.rm = TRUE))

# Boxplot of height x treatment
height_avg <- height_avg %>%
  mutate(across(drought_period, factor, levels=c("Pre-Drought","Drought","Post-Drought"))) # re-ordering drought period
height <- height %>%
  mutate(across(drought_period, factor, levels=c("Pre-Drought","Drought","Post-Drought"))) # re-ordering drought period
png("height_overall.png", units="in", width=10, height=6, res=300)
ggplot(height, aes(x = treatment, y = plant_height_cm)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  facet_grid(.~drought_period) +
  labs(x = "Treatment", y = "Plant Height (cm)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=2)) +
  #geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Boxplot of height w/o drought period
png("height_temp.png", units="in", width=8, height=6, res=300)
ggplot(height, aes(x = treatment, y = plant_height_cm)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Treatment", y = "Plant Height (cm)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  #geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Boxplot of height only across drought periods
png("height_time.png", units="in", width=8, height=6, res=300)
ggplot(height, aes(x = drought_period, y = plant_height_cm)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Timing of Drought Period", y = "Plant Height (cm)") +
  #scale_fill_manual(values = c("olivedrab")) +
  #scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
  #                 labels=c("Ambient" = "Ambient",
  #                          "Ambient Drought" = "Drought",
  #                          "Irrigated Control" = "Irrigated \n Control",
  #                          "Warm" = "Warmed",
  #                          "Warm Drought" = "Warmed & \n Drought"),
  #                 guide = guide_axis(n.dodge=2)) +
  #geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Barplot of height x treatment
png("height_overall_bar.png", units="in", width=10, height=6, res=300)
ggplot(height_avg, aes(x = treatment, y = avg_height)) +
  facet_grid(.~drought_period) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Plant Height (cm)") +
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
png("height_temp_bar.png", units="in", width=8, height=6, res=300)
ggplot(height_avg2, aes(x = treatment, y = avg_height)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Plant Height (cm)") +
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
ggplot(height_avg3, aes(x = drought_period, y = avg_height)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_height - se, ymax = avg_height + se), width = 0.2,
                position = "identity") +
  labs(x = "Timing of Drought Period", y = "Plant Height (cm)") +
  scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  #scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
  #                 labels=c("Ambient" = "Ambient",
  #                          "Ambient Drought" = "Drought",
  #                          "Irrigated Control" = "Irrigated \n Control",
  #                          "Warm" = "Warmed",
  #                          "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "none")
dev.off()

# Boxplot of height x treatment x gall
png("height_galling.png", units="in", width=9, height=6, res=300)
ggplot(height, aes(x = treatment, y = plant_height_cm, fill=gall_present)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  labs(x = "Treatment", y = "Plant Height (cm)") +
  scale_fill_manual(values = c("olivedrab4","darkseagreen2"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought"),
                   guide = guide_axis(n.dodge=1)) +
  theme(legend.position = "right")
dev.off()

# barplot of treatment X height X galling
png("height_galling_bar.png", units="in", width=9, height=6, res=300)
ggplot(height_avg4, aes(x = treatment, y = avg_height, fill=gall_present)) +
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

# TITLE:          REX: Gall chamber volume & count plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Kristin Wolford, Emily Parker, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_insect L1 folder
# DATA OUTPUT:    Plots
# PROJECT:        REX
# DATE:           Jan 2022; updated June 2024 w/ gall weight

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
count <- read.csv(file.path(dir, "T7_warmx_insect/L1/T7_warmx_Soca_gall_chmb_count_L1.csv"))
vol <- read.csv(file.path(dir, "T7_warmx_insect/L1/T7_warmx_Soca_gall_chmb_vol_L1.csv"))
weight <- read.csv(file.path(dir, "T7_warmx_insect/L1/T7_warmx_Soca_galls_weight_L1.csv"))

# Taking averages
count_avg <- count %>%
  group_by(treatment) %>%
  summarize(avg_count = mean(num_of_chambers, na.rm = TRUE),
            se = std.error(num_of_chambers, na.rm = TRUE))
vol_avg <- vol %>%
  group_by(treatment) %>%
  summarize(avg_vol = mean(chamber_volume_mm3, na.rm = TRUE),
            se = std.error(chamber_volume_mm3, na.rm = TRUE))
weight_avg <- weight %>%
  group_by(Climate_Treatment) %>%
  summarize(avg_weight = mean(Dried_Weight, na.rm = TRUE),
            se = std.error(Dried_Weight, na.rm = TRUE))

# Boxplot of chamber counts
png("gall_chamber_count_box.png", units="in", width=8, height=6, res=300)
ggplot(count, aes(x = treatment, y = num_of_chambers)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Treatment", y = "Number of Chambers per Gall") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Barplot of chamber counts
png("gall_chamber_count_bar.png", units="in", width=8, height=6, res=300)
ggplot(count_avg, aes(x = treatment, y = avg_count)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_count - se, ymax = avg_count + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Number of Chambers per Gall") +
  #scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "none")
dev.off()

# Boxplot of chamber volumne
png("gall_chamber_vol_box.png", units="in", width=8, height=6, res=300)
ggplot(vol, aes(x = treatment, y = chamber_volume_mm3)) +
  geom_boxplot(color = "black", outlier.shape = NA, fill = "olivedrab") +
  labs(x = "Treatment", y = "Gall chamber volume (mm^3)") +
  #scale_fill_manual(values = c("olivedrab")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  geom_jitter(aes(alpha=0.8), color="olivedrab", fill="olivedrab", shape=16, size=2) +
  theme(legend.position = "none")
dev.off()

# Barplot of chamber counts
png("gall_chamber_vol_bar.png", units="in", width=8, height=6, res=300)
ggplot(vol_avg, aes(x = treatment, y = avg_vol)) +
  geom_bar(position = "identity", stat = "identity", color = "black",fill = "olivedrab") +
  geom_errorbar(aes(ymin = avg_vol - se, ymax = avg_vol + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Gall chamber volume (mm^3)") +
  #scale_fill_manual(values = c("#a6bddb", "#fb6a4a")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "none")
dev.off()


# Plot of means for gall weight
ggplot(weight_avg, aes(x = Climate_Treatment, y = avg_weight)) +
  geom_pointrange(aes(ymin = avg_weight - se, ymax = avg_weight + se),pch=21,size=1, fill = "olivedrab4") +
  labs(x = "Treatment", y = "Gall Weight (g)") +
  #scale_fill_manual(values = c("olivedrab4", "darkseagreen2"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "right")

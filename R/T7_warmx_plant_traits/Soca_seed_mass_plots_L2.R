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
library(ggpubr)

# Set ggplot2 plotting
# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(14))
theme_update(axis.text.x = element_text(size = 12),
             axis.text.y = element_text(size = 16),
             axis.title = element_text(size=16,face="bold"))

# Set working directory
dir <- setwd("/Users/Emily/Documents/R/Goldenrod Project")
dir<-Sys.getenv("DATA_DIR") # Kara


# Read in data
mass <- read.csv(file.path(dir, "L1/T7_warmx_soca_seeds_mass_L1.csv"))
mass <- read.csv(file.path(dir, "T7_warmx_plant_traits/L1/T7_warmx_soca_seeds_mass_L1.csv")) # Kara


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
  scale_fill_manual(values = c("purple4", "plum1"), name="Gall Presence",labels=c("Gall","No Gall")) +
  scale_x_discrete(limits = c("Irrigated Control", "Ambient", "Ambient Drought", "Warm", "Warm Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Ambient Drought" = "Drought",
                            "Irrigated Control" = "Irrigated \n Control",
                            "Warm" = "Warmed",
                            "Warm Drought" = "Warmed & \n Drought")) +
  theme(legend.position = "right")
dev.off()



####### probability of having seed + weight of seed #######
## making binary response for if it had a seed or not
mass_binom <- mass %>%
  mutate_at(vars(contains('Seeds_Mass')), ~1 * (. != 0))
mass_binom$Seeds_Mass[mass_binom$Seeds_Mass == 1] <- "Seed"
mass_binom$Seeds_Mass[mass_binom$Seeds_Mass == 0] <- "No Seed"
mass_binom_sum <- mass_binom %>%
  group_by(Treatment,Rep,Footprint,Subplot,Climate_Treatment,Galling_Status,Seeds_Mass) %>%
  count(Climate_Treatment,Galling_Status,Seeds_Mass) %>%
  group_by(Treatment,Rep,Footprint,Subplot,Climate_Treatment,Galling_Status) %>%
  mutate(n = n/sum(n)) %>%
  group_by(Climate_Treatment,Galling_Status,Seeds_Mass) %>%
  summarize(mean_n = mean(n),
            se = std.error(n))
mass_binom_seed <- mass_binom_sum %>%
  filter(Seeds_Mass == "Seed")
# plot
level_order <- c("Irrigated Control","Ambient","Ambient Drought","Warm","Warm Drought") 
binom_plot <- ggplot(mass_binom_seed, aes(x = factor(Climate_Treatment, level = level_order), y = mean_n, fill=Galling_Status)) +
  geom_pointrange(aes(ymin = mean_n - se, ymax = mean_n + se), ,pch=21,size=1,position=position_dodge(0.2)) +
  labs(x = NULL, y = "Probability of having a seed") +
  scale_x_discrete(labels=c("Ambient" = "Ambient", "Warm" = "Warmed",
                            "Ambient Drought" = "Drought", "Irrigated Control" = "Irrigated\nControl",
                            "Warm Drought" = "Warmed &\nDrought")) +
  scale_fill_manual(name="Galling status",
                    labels=c("Galled","Non-galled"),
                    values = c("purple4", "plum2")) +
  theme_bw(14) +
  theme(axis.title = element_text(face="bold"))


## seed weight
mass_cond <- mass %>%
  group_by(Climate_Treatment, Galling_Status) %>%
  summarize(avg_weight = mean(Seeds_Mass, na.rm = TRUE),
            se = std.error(Seeds_Mass, na.rm = TRUE))
# plot
cond_plot <- ggplot(mass_cond, aes(x = factor(Climate_Treatment, level = level_order), y = avg_weight, fill=Galling_Status)) +
  geom_pointrange(aes(ymin = avg_weight - se, ymax = avg_weight + se),pch=21,size=1,position=position_dodge(0.2)) +
  labs(x = NULL, y = "Seed weight (g)", title=NULL) +
  scale_x_discrete(labels=c("Ambient" = "Ambient", "Warm" = "Warmed",
                            "Ambient Drought" = "Drought", "Irrigated Control" = "Irrigated\nControl",
                            "Warm Drought" = "Warmed &\nDrought")) +
  scale_fill_manual(name="Galling status",
                    labels=c("Galled","Non-galled"),
                    values = c("purple4", "plum2")) +
  theme_bw(14) +
  theme(axis.title = element_text(face="bold"))

# plotting binary & conditional plot on same figure
png("seed_weight.png", units="in", width=11, height=5, res=300)
ggarrange(binom_plot,cond_plot,
          ncol = 2, common.legend = T, legend="right",widths = c(1, 1))
dev.off()

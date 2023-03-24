# TITLE:          REX: Zarnetske plots ANPP analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond, Jordan Zapata
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_ANPP L1 folder
# DATA OUTPUT:    Analyses
# PROJECT:        REX
# DATE:           June 2022


### Preliminary steps & setting up data ###
# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(lmerTest)
library(car)
library(bbmle)
library(sjPlot)
library(emmeans)
library(stats)
library(plotrix)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")

# Read in data
anpp_2021 <- read.csv(file.path(dir, "T7_ANPP/L1/T7_warmx_ANPP_2021_L1.csv"))
anpp_2022 <- read.csv(file.path(dir, "T7_ANPP/L1/T7_warmx_ANPP_2022_L1.csv"))

# delete unneeded columns
anpp_2021 = subset(anpp_2021, select = -c(Replicate, Unique_ID))

# add year column
anpp_2021$Year <- "2021"
anpp_2022$Year <- "2022"

# merge both years
anpp <- rbind(anpp_2021,anpp_2022)

# Main questions for analyses:
# Does biomass differ between the climate treatment (warmed, drought, warmed+drought, ambient, irr. control)?
# Does the biomass of native + exotic species differ between these climate treatments?
# Also comparisons of growth form + species-specific comparisons

# remove climate treatment plots we aren't interested in here (i.e., insecticide plots)
unique(anpp$Subplot_Descriptions)
anpp <- anpp %>%
  filter(!(Subplot_Descriptions=='drought_insecticide' |
             Subplot_Descriptions=='insecticide' |
             Subplot_Descriptions=='warmed_drought_insecticide' |
             Subplot_Descriptions=='warmed_insecticide'))


#### Data exploration ###
# checking raw data
hist(anpp$Dried_Plant_Biomass_g)
qqnorm(anpp$Dried_Plant_Biomass_g)

# checking model with log transformation (untransformed data was right skewed)
# is it normally distributed? yes
m1 <- lm(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code + Year, data=anpp)
hist(resid(m1))
qqnorm(resid(m1))
shapiro.test(resid(m1))
# homogeneity of variance? yes if p >0.05 (no significant difference between the group variances)
leveneTest(residuals(m1) ~ anpp$Subplot_Descriptions) # yes
leveneTest(residuals(m1) ~ anpp$Species_Code) # no, but not surprising that species differ
leveneTest(residuals(m1) ~ anpp$Year) # yes

# checking model with sqrt transformation (untransformed data was right skewed)
# is it normally distributed? yes
m1.2 <- lm(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code + Year, data=anpp)
hist(resid(m1.2))
qqnorm(resid(m1.2))
shapiro.test(resid(m1.2))
# homogeneity of variance? yes if p >0.05 (no significant difference between the group variances)
leveneTest(residuals(m1.2) ~ anpp$Subplot_Descriptions) # yes
leveneTest(residuals(m1.2) ~ anpp$Species_Code) # no
leveneTest(residuals(m1.2) ~ anpp$Year) # yes



### Model exploration ###
## species-level biomass (plot-summed biomass models down below, and are likely more useful)
# year
m2 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
m3 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
anova(m2,m3)
AICctab(m2,m3) # m2
anova(m2)
summary(m2)
# origin
m4 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + origin + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
m5 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * origin + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
anova(m4,m5)
AICctab(m4,m5) # m4
anova(m4)
summary(m4)
# species
m6 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
m7 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * Species_Code + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
anova(m6,m7)
AICctab(m6,m7) # going with m7
anova(m7)
summary(m7)
# growth habit
m8 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + growth_habit + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
m9 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * growth_habit + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
anova(m8,m9)
AICctab(m8,m9) # m8
anova(m8)
summary(m8)
# rhizomatous
m10 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + MH_rhizomatous_suggestion + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
m11 <- lmer(sqrt(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * MH_rhizomatous_suggestion + Year + (1|Rep/Footprint_Location), data=anpp, REML=F)
anova(m10,m11)
AICctab(m10,m11) # m10
anova(m10)
summary(m10)


# summing data to plot-level for analyses
anpp_plot <- anpp %>%
  group_by(Rep, Footprint_Location, Subplot_Location, Subplot_Descriptions, Year) %>% 
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE))
hist(anpp_plot$sum_biomass)
# overall biomass
m_plot1 <- lmer(sum_biomass ~ Subplot_Descriptions + (1|Rep/Footprint_Location), data=anpp_plot, REML=F)
m_plot2 <- lmer(sum_biomass ~ Subplot_Descriptions + Year + (1|Rep/Footprint_Location), data=anpp_plot, REML=F)
m_plot3 <- lmer(sum_biomass ~ Subplot_Descriptions * Year + (1|Rep/Footprint_Location), data=anpp_plot, REML=F)
anova(m_plot1,m_plot2)
anova(m_plot2,m_plot3)
hist(resid(m_plot2))
shapiro.test(resid(m_plot2))
anova(m_plot2)
summary(m_plot2)
anpp_plot <- within(anpp_plot, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
# one year
anpp_plot_yr <- anpp_plot %>%
  filter(Year == "2021") # changed to the year I wanted & re-ran the model each time
m_plot_yr <- lmer(sum_biomass ~ Subplot_Descriptions + (1|Rep/Footprint_Location), data=anpp_plot_yr, REML=F)
hist(resid(m_plot_yr))
shapiro.test(resid(m_plot_yr))
anova(m_plot_yr)
summary(m_plot_yr)
anpp_plot_yr <- within(anpp_plot_yr, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated_control"))


# growth habit
anpp_plot_growth <- anpp %>%
  group_by(Rep, Footprint_Location, Subplot_Location, Subplot_Descriptions, growth_habit, Year) %>% 
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE)) %>%
  filter(!(growth_habit == "Vine" | growth_habit == "" | growth_habit == "Shurb" | growth_habit == "Subshrub/vine"))
# one year
anpp_plot_growth_yr <- anpp_plot_growth %>%
  filter(Year == "2022" & growth_habit == "Forb") # changed to the year I wanted & re-ran the model each time
m_plot_growth2 <- lmer(sum_biomass ~ Subplot_Descriptions + (1|Rep/Footprint_Location), data=anpp_plot_growth_yr, REML=F)
anova(m_plot_growth2)
summary(m_plot_growth2)
anpp_plot_growth_yr <- within(anpp_plot_growth_yr, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated_control"))


### species specific models ###
anpp_spp <- anpp %>%
  filter(Species_Code == "POAPR")
## list of species to test: ACHMI, SOOGR, ASTPI, ERIAN, SOOCA, ASTSA, HYPPE, PHLPR, DACGL, POAPR, AGRRE, POACO, TRFPR
m_spp <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_spp, REML=F)
summary(m_spp)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))

# summing to plot level
anpp_spp_plot <- anpp_spp %>%
  group_by(Treatment, Replicate, Footprint_Location, Subplot_Location, Subplot_Descriptions) %>% 
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE))
m_spp_plot <- lmer(sum_biomass~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_spp_plot, REML=F)
summary(m_spp_plot)
anpp_spp_plot <- within(anpp_spp_plot, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))


### origin specific models ###
anpp_org <- anpp %>%
  filter(origin == "Native")
m_org <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_org, REML=F)
summary(m_org)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_org <- within(anpp_org, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_org <- within(anpp_org, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_org <- within(anpp_org, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_org <- within(anpp_org, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_org <- within(anpp_org, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))


### growth form specific models ###
anpp_growth <- anpp %>%
  filter(growth_habit == "Graminoid")
m_growth <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_growth, REML=F)
summary(m_growth)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))

# summing to plot-level
anpp_growth_2 <- anpp_growth %>%
  group_by(Treatment, Replicate, Footprint_Location, Subplot_Location, Subplot_Descriptions) %>% 
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE))
m_growth2 <- lmer(log(sum_biomass) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_growth_2, REML=F)
summary(m_growth2)
anpp_growth_2 <- within(anpp_growth_2, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))


### rhizomatous specific models ###
anpp_rhiz <- anpp %>%
  filter(MH_rhizomatous_suggestion == "highly_rhizomatous")
m_rhiz <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_rhiz, REML=F)
summary(m_rhiz)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_rhiz <- within(anpp_rhiz, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_rhiz <- within(anpp_rhiz, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_rhiz <- within(anpp_rhiz, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_rhiz <- within(anpp_rhiz, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_rhiz <- within(anpp_rhiz, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))


# looking at data between treatments
ggplot(anpp, aes(x=origin, y=Dried_Plant_Biomass_g)) + 
  facet_wrap(.~Subplot_Descriptions) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  theme_classic()

anpp2 <- anpp %>%
  group_by(Field_Loc_Code, Subplot_Descriptions) %>%
  summarize(plot_sum = sum(Dried_Plant_Biomass_g)) %>%
  group_by(Subplot_Descriptions) %>%
  summarize(mean_biomass = mean(plot_sum, na.rm = TRUE),
            se = std.error(plot_sum, na.rm = TRUE))

ggplot(anpp2, aes(x=Subplot_Descriptions, y=mean_biomass)) + 
  geom_bar(position = "dodge", stat = "identity", col = "black") +
  geom_errorbar(aes(ymin = mean_biomass - se, ymax = mean_biomass + se), width = 0.2,
                position = position_dodge(0.9)) +
  theme_classic()



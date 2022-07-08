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
anpp <- read.csv(file.path(dir, "T7_ANPP/L1/T7_Zarnetske_ANPP_L1.csv"))

# Main questions for analyses:
# Does biomass differ between the climate treatment (warmed, drought, warmed+drought, ambient, irr. control)?
# Does the biomass of native + exotic species differ between these climate treatments?
# Also comparisons of growth form + species-specific comparisons

# remove climate treatment plots we aren't interested in here (i.e., insecticide plots)
anpp <- anpp %>%
  filter(!(Subplot_Descriptions=='drought_insecticide' |
             Subplot_Descriptions=='insecticide' |
             Subplot_Descriptions=='warmed_drought_insecticide' |
             Subplot_Descriptions=='warmed_insecticide'))

# remove species we aren't interested in (i.e., unknowns and unsorted plant material)
anpp <- anpp %>%
  filter(!(Species_Code=='UNASTAE' |
             Species_Code=='UNDIC' |
             Species_Code=='UNFABAE' |
             Species_Code=='UNSRT' |
             Species_Code=='CORSPP' | # removing cornus species here (for now) bc it is not designated as native/exotic
             origin == 'Both')) # removing species listed as both native and exotic


#### Data exploration ###
# checking raw data
hist(anpp$Dried_Plant_Biomass_g)
qqnorm(anpp$Dried_Plant_Biomass_g)

# checking model with log transformation (untransformed data was right skewed)
# is it normally distributed? yes
m1 <- lm(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code, data=anpp)
hist(resid(m1))
qqnorm(resid(m1))
plot(m1)
# homogeneity of variance? yes, p >0.05 (no significant difference between the group variances)
leveneTest(residuals(m1) ~ anpp$Subplot_Descriptions)
leveneTest(residuals(m1) ~ anpp$Species_Code)


### Model exploration ###
m2 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + origin + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m3 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * origin + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m4 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m5 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * Species_Code + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m6 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Species_Code) + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m7 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + growth_habit + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m8 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * growth_habit + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m9 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + MH_rhizomatous_suggestion + (1|Replicate/Footprint_Location), data=anpp, REML=F)
m10 <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * MH_rhizomatous_suggestion + (1|Replicate/Footprint_Location), data=anpp, REML=F)
AICctab(m2, m3, m4, m5, m6, m7, m8, m9, m10)
summary(m7)
emmeans(m5, list(pairwise ~ Subplot_Descriptions * Species_Code), adjust = "tukey")
emmip(m5, Species_Code ~ Subplot_Descriptions)

# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp <- within(anpp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp <- within(anpp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp <- within(anpp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp <- within(anpp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp <- within(anpp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))
anpp <- within(anpp, growth_habit <- relevel(factor(growth_habit), ref = "Graminoid"))
anpp <- within(anpp, growth_habit <- relevel(factor(growth_habit), ref = "Forb"))
anpp <- within(anpp, growth_habit <- relevel(factor(growth_habit), ref = "Vine"))
anpp <- within(anpp, MH_rhizomatous_suggestion <- relevel(factor(MH_rhizomatous_suggestion), ref = "non_rhizomatous"))
anpp <- within(anpp, MH_rhizomatous_suggestion <- relevel(factor(MH_rhizomatous_suggestion), ref = "slightly_rhizomatous"))
anpp <- within(anpp, MH_rhizomatous_suggestion <- relevel(factor(MH_rhizomatous_suggestion), ref = "highly_rhizomatous"))

# quick anova to check 
two.way <- aov(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions * origin, data = anpp)
summary(two.way)
two.way2 <- aov(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + Species_Code, data = anpp)
summary(two.way2)


### species specific models ###
anpp_spp <- anpp %>%
  filter(Species_Code == "TRFPR")
## list of species to test: ACHMI, SOOGR, ASTPI, ERIAN, SOOCA, ASTSA, HYPPE, PHLPR, DACGL, POAPR, AGRRE, POACO, TRFPR
m_spp <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_spp, REML=F)
summary(m_spp)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_spp <- within(anpp_spp, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))


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
  filter(growth_habit == "Forb")
m_growth <- lmer(log(Dried_Plant_Biomass_g) ~ Subplot_Descriptions + (1|Replicate/Footprint_Location), data=anpp_growth, REML=F)
summary(m_growth)
# re-leveling the dataframe & re-running model for post-hoc comparisions
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "drought"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "irrigated"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "warmed_drought"))
anpp_growth <- within(anpp_growth, Subplot_Descriptions <- relevel(factor(Subplot_Descriptions), ref = "ambient"))


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



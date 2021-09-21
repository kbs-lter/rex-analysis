# TITLE:          REX: Drought VOC analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_VOC L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           July 2021

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(broom)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021drought_L1.csv"))



#### VOC Composition - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:359]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
ab.div<-adonis2(ab.dist~Treatment, data=voc_transpose, permutations = 999, method="bray", strata="Rep")
ab.div

# pairwise comparisons of permanova
ab.pair<-pairwise.adonis2(ab.dist~Treatment, data=voc_transpose, method="bray", strata="Rep")
ab.pair

# testing for homogeneity of dispersion among groups
# >0.05 meets assumption of adonis permanova
# i.e. adonis may give sig. p-value even if groups overlap because within-group data is heterogenous
# this test looks to see if group data is heterogenous; p>0.05 means it is not, and therefore 
# we can assume adonis results are "real" and not a result of heterogenous dispersion
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
dispersion
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T)



#### VOC Abundance - Mixed model ####
# Data exploration
voc_transpose$rowsums <- rowSums(voc_transpose[2:359])
descdist(voc_transpose$rowsums, discrete = FALSE)
hist(voc_transpose$rowsums)
qqnorm(voc_transpose$rowsums)
shapiro.test(voc_transpose$rowsums)
# kinda right skewed, going to try a few transformations

# square root transformation
voc_transpose$sqrt_rowsums <- sqrt(voc_transpose$rowsums)
descdist(voc_transpose$sqrt_rowsums, discrete = FALSE)
hist(voc_transpose$sqrt_rowsums)
qqnorm(voc_transpose$sqrt_rowsums)
shapiro.test(voc_transpose$sqrt_rowsums)

# cubed root transformation
voc_transpose$cubed_rowsums <- (voc_transpose$rowsums)^(1/3)
descdist(voc_transpose$cubed_rowsums, discrete = FALSE)
hist(voc_transpose$cubed_rowsums)
qqnorm(voc_transpose$cubed_rowsums)
shapiro.test(voc_transpose$cubed_rowsums)
# none of these are great, come back to this 

# comparison with other models
m1 <- lmer(cubed_rowsums ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m1)
emmeans(m1, list(pairwise ~ Treatment), adjust = "tukey")
AICctab(m1, m2, weights=T)
# emmeans is a lot different than summary (?) so re-leveling to get pairwise comparisons
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
summary(m1)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
summary(m1)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated"))
summary(m1)

# specific compound test
m2 <- lmer(Caryophyllene ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated"))
summary(m2)
# specific compound test
m2 <- lmer(Ethanone..1..4.ethylphenyl.. ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated"))
summary(m2)



#### VOC Abundance - ANOVA ####
# overall treatment differences
voc_transpose$rowsums <- rowSums(voc_transpose[2:359])
anova1 <- aov(rowsums~Treatment, data = voc_transpose)
summary(anova1)
TukeyHSD(anova1)

# specific compound differences
for (i in 2:359){
  column <- names(voc_transpose[i])
  anova <- broom::tidy(aov(voc_transpose[,i] ~ Treatment, data = voc_transpose))
  
  # only want aov with P < 0.07 printed
  if(anova$p.value[1] < 0.07) {
    
    print(column)
    print(anova)
  }
}

# pairwise comparisons
for (i in 2:359){
  column <- names(voc_transpose[i])
  anova <- aov(voc_transpose[,i] ~ Treatment, data = voc_transpose)
  tukey <- TukeyHSD(anova)
  
  # only want tukey with P < 0.07 printed
  if(any(tukey$Treatment[, "p adj"] < 0.07)) {
    
    print(column)
    print(setNames(tukey, column))
    
  }
}


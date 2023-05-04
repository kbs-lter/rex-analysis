# TITLE:          REX: 2022 VOC analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_VOC L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           Aug 2022

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(broom)
library(fitdistrplus)
library(lmerTest)
library(emmeans)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_named_VOC_2022_L1.csv"))
#voc_transpose <- voc_transpose[!grepl("Irrigated_Control", voc_transpose$Treatment),] # removing irrigated control (just to test)



#### VOC Composition (treatment) - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:429]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose$Rep)
ab.div<-adonis2(ab.dist~Treatment, data=voc_transpose, permutations = perm, method="bray")
ab.div
# note: also can run models with rep as additive or interactive effect w/ treatment

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
plot(dispersion, hull=F, ellipse=T, label=F)


#### VOC Composition (rep) - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:429]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose$Treatment)
ab.div<-adonis2(ab.dist~Rep, data=voc_transpose, permutations = perm, method="bray")
ab.div

# pairwise comparisons of permanova
ab.pair<-pairwise.adonis2(ab.dist~Rep, data=voc_transpose, method="bray")
ab.pair

# testing for homogeneity of dispersion among groups
# >0.05 meets assumption of adonis permanova
# i.e. adonis may give sig. p-value even if groups overlap because within-group data is heterogenous
# this test looks to see if group data is heterogenous; p>0.05 means it is not, and therefore 
# we can assume adonis results are "real" and not a result of heterogenous dispersion
dispersion<-betadisper(ab.dist, group=voc_transpose$Rep)
dispersion
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T)



### VOC composition - ANOSIM (analysis of similarity) ###
ab = voc_transpose[,2:429]
mat_ab = as.matrix(ab)
ano = anosim(mat_ab, voc_transpose$Treatment, distance = "bray", permutations = 999)
ano


# note: update the section below with new matrix information
#### VOC Abundance - Mixed model ####
# Data exploration on raw data
voc_transpose$rowsums <- rowSums(voc_transpose[2:429])
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
hist(resid(m1))
shapiro.test(resid(m1))
# emmeans is a lot different than summary (?) so re-leveling to get pairwise comparisons
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed_Drought"))


# specific compound test
m2 <- lmer(.beta..Myrcene ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
anova(m2)
summary(m2)
# pairwise comparisons
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Ambient_Control"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))

m3 <- lmer(X4.Hexen.1.ol..acetate ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
anova(m3)
summary(m3)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))


#### VOC Abundance - ANOVA ####
# overall treatment differences
anova1 <- aov(rowsums~Treatment, data = voc_transpose)
summary(anova1)
TukeyHSD(anova1)

# specific compound differences
for (i in 2:429){
  column <- names(voc_transpose[i])
  anova <- broom::tidy(aov(voc_transpose[,i] ~ Treatment, data = voc_transpose))
  
  # only want aov with P < 0.05 printed
  if(anova$p.value[1] < 0.05) {
    
    print(column)
    print(anova)
  }
}



# pairwise comparisons
for (i in 2:429){
  column <- names(voc_transpose[i])
  anova <- aov(voc_transpose[,i] ~ Treatment, data = voc_transpose)
  tukey <- TukeyHSD(anova)
  
  # only want tukey with P < 0.05 printed
  if(any(tukey$Treatment[, "p adj"] < 0.05)) {
    
    print(column)
    print(setNames(tukey, column))
    
  }
}


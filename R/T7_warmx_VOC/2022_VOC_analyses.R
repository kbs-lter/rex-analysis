# TITLE:          REX: 2022 VOC analyses
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
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
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_total_VOC_2022_L1.csv"))
#voc_transpose <- voc_transpose[!grepl("Irrigated_Control", voc_transpose$Treatment),] # removing irrigated control (just to test)


#### VOC Composition - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:1493]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose$Rep)
ab.div<-adonis2(ab.dist~Treatment, data=voc_transpose, permutations = perm, method="bray")
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


# note: update the section below with new matrix information
#### VOC Abundance - Mixed model ####
# Data exploration
voc_transpose$rowsums <- rowSums(voc_transpose[2:1455])
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
m1 <- lmer(rowsums ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m1)
# emmeans is a lot different than summary (?) so re-leveling to get pairwise comparisons
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
summary(m1)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
summary(m1)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated"))
summary(m1)

# specific compound test
m2 <- lmer(X.Z.Z...alpha..Farnesene ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m2)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated_Control"))

m3 <- lmer(endo.Borneol ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
summary(m3)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Irrigated_Control"))


#### VOC Abundance - ANOVA ####
# overall treatment differences
voc_transpose$rowsums <- rowSums(voc_transpose[2:1455])
anova1 <- aov(rowsums~Treatment, data = voc_transpose)
summary(anova1)
TukeyHSD(anova1)

# specific compound differences
for (i in 2:1455){
  column <- names(voc_transpose[i])
  anova <- broom::tidy(aov(voc_transpose[,i] ~ Treatment, data = voc_transpose))
  
  # only want aov with P < 0.05 printed
  if(anova$p.value[1] < 0.05) {
    
    print(column)
    print(anova)
  }
}

# Benzaldehyde..3.ethyl. 0.0000568
# endo.Borneol 0.000223
# Germacrene.D 0.0288
# Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester 0.00929
# .alpha..Bourbonene 0.00305
# X.Z.Z...alpha..Farnesene 0.0123
# o.Xylene 0.0302
# .beta..Myrcene 0.0246
# Camphor 0.0420
# dl.Menthol 0.0107
# Ethanone..1..4.ethylphenyl.. 0.00613
# Pyridine 0.0485
# Salicylic.acid..tert..butyl.ester 0.00468

# pairwise comparisons
for (i in 2:1455){
  column <- names(voc_transpose[i])
  anova <- aov(voc_transpose[,i] ~ Treatment, data = voc_transpose)
  tukey <- TukeyHSD(anova)
  
  # only want tukey with P < 0.05 printed
  if(any(tukey$Treatment[, "p adj"] < 0.05)) {
    
    print(column)
    print(setNames(tukey, column))
    
  }
}


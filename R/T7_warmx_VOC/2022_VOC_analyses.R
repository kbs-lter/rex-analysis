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
library(knitr)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_named_VOC_2022_L1.csv"))
# add check.names=F into the read.csv function to have full compound names (i.e., w/o periods for spaces)
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

# how do specific compounds respond to these treatments?
# using indicator species approach
# https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html#:~:text=Indicator%20species%20are%20often%20determined,types%2C%20disturbance%20states%2C%20etc.
indval = multipatt(ab, voc_transpose$Treatment, 
                   control = how(nperm=999)) 
summary(indval)

# same test, but this time not allowing for the grouping of treatments
indval2 = multipatt(ab, voc_transpose$Treatment, duleg=T,
                   control = how(nperm=999)) 
summary(indval2)

# same test, but this time allowing for max of 2 treatments to be grouped
indval3 = multipatt(ab, voc_transpose$Treatment, max.order=2,
                   control = how(nperm=999)) 
summary(indval3)

# is the association random? last column is p-value
prefsign = signassoc(ab, cluster=voc_transpose$Treatment, alternative = "two.sided", 
                     control = how(nperm=199)) 

# compound combinations as indicators?
# note: this output it really large
ab.comb = combinespecies(ab, max.order = 2)$XC
indvalspcomb = multipatt(ab.comb, voc_transpose$Treatment, duleg = TRUE, 
                         control = how(nperm=999))
summary(indvalspcomb, alpha=0.05)

# why were these compounds chosen? the code below gives A and B values
# A=1 means that this compound was only found in that treatment, while B=1 means that this compound
# was found across all reps of that treatment
# so A=1, B=0.3 for warmed means it was only found in warmed plots, but not all warmed plots
summary(indval, indvalcomp=TRUE)
# can also change the alpha level
summary(indval, alpha=0.08, indvalcomp=TRUE)

# presence absence testing
ab.pa = ifelse(ab>0,1,0)
phi = multipatt(ab.pa, voc_transpose$Treatment, func = "r", 
                control = how(nperm=999)) 
summary(phi)
round(head(phi$str),3)

### testing for homogeneity of dispersion among groups
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
shapiro.test(voc_transpose$sqrt_rowsums) # pretty close

# cubed root transformation
voc_transpose$cubed_rowsums <- (voc_transpose$rowsums)^(1/3)
descdist(voc_transpose$cubed_rowsums, discrete = FALSE)
hist(voc_transpose$cubed_rowsums)
qqnorm(voc_transpose$cubed_rowsums)
shapiro.test(voc_transpose$cubed_rowsums) # best 

# comparison with other models
m1 <- lmer(cubed_rowsums ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
anova(m1)
summary(m1)
hist(resid(m1))
shapiro.test(resid(m1))
# comparisons
emmeans(m1, list(pairwise ~ Treatment), adjust = "tukey")



#### VOC Abundance - ANOVA ####
# overall treatment differences
anova1 <- aov(rowsums~Treatment, data = voc_transpose)
summary(anova1)
TukeyHSD(anova1)



### VOC abundance - specific compound differences ###
# pulls out compounds that had their abundances significantly different between treatments (p<0.05)
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



### VOC abundance - mixed models for compounds identified above ###
# Beta Myrcene
m2 <- lmer(.beta..Myrcene ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
anova(m2)
summary(m2)
# pairwise comparisons
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Ambient_Control"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))

# 4-Hexen-1-ol acetate
m3 <- lmer(X4.Hexen.1.ol..acetate ~ Treatment + (1|Rep), data = voc_transpose, REML=FALSE)
anova(m3)
summary(m3)
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
voc_transpose <- within(voc_transpose, Treatment <- relevel(factor(Treatment), ref = "Drought"))




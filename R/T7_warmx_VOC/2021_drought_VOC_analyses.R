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



#### PERMANOVA ####
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



#### ANOVA ####
# overall treatment differences
voc_transpose$rowsums <- rowSums(voc_transpose[2:359])
anova1 <- aov(rowsums~Treatment, data = voc_transpose)
summary(anova1)
tuk <- TukeyHSD(anova1)

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

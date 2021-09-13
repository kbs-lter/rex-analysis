# TITLE:          REX: Pre-drought VOC analyses
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

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021predrought_L1.csv"))



#### PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:300]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(1)
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
voc_transpose$rowsums <- rowSums(voc_transpose[2:300])
anova_one_way <- aov(rowsums~Treatment+Rep, data = voc_transpose)
summary(anova_one_way)
TukeyHSD(anova_one_way)

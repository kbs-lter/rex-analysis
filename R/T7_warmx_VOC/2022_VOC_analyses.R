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
# removing rep 1 - data is weird
voc_transpose_rm <- voc_transpose %>%
  filter(!(Rep == 1))
# removing samples w/ abnormally high abundances - see notes below in "abundance" category
voc_transpose_rm <- voc_transpose_rm %>%
  filter(!(Unique_ID == 79)) %>%
  filter(!(Unique_ID == 39)) %>%
  filter(!(Unique_ID == 62))
# how many individuals were measure per treatment + rep?
# I used this info in the leaf biomass L2 script to determine total biomass for only measured individuals
voc_transpose_rm %>% 
  count(Treatment,Rep)
voc_biomass <- read.csv(file.path(dir, "T7_warmx_VOC/L1/VOC_biomass_2022_L1.csv"))
# making biomass treatments match voc data
voc_biomass$Treatment[voc_biomass$Treatment == "Ambient"] <- "Ambient_Control"
voc_biomass$Treatment[voc_biomass$Treatment == "Irrigated"] <- "Irrigated_Control"

# merge voc w/ biomass data
voc_transpose_rm <- left_join(voc_transpose_rm,voc_biomass,by=c("Treatment","Rep"))
# divide voc abundances by indiv plant biomass per treatment/rep
# first remove compound name column & meta info columns
voc_sample_names <- voc_transpose_rm[,1, drop=FALSE]
voc_meta_info <- voc_transpose_rm[,430:436, drop=FALSE]
voc_transpose_rm2 = subset(voc_transpose_rm, select = -c(Sample_ID,Unique_ID,Rep,Footprint,Subplot,Treatment,Notes,Weight_g))
# divide
voc_weighted_abun <- voc_transpose_rm2/voc_transpose_rm2[,429]
# remerging with meta info
voc_transpose_rm3 <- cbind(voc_sample_names,voc_weighted_abun,voc_meta_info)
# removing indiv weight column 
voc_transpose_rm3 = subset(voc_transpose_rm3, select = -c(Weight_indiv_g))
# so now when I want to do analyses for abundances/individual, I use the voc_transpose_rm3 dataframe



#### VOC Composition (treatment) - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose_rm3[,2:429]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose_rm3$Rep)
ab.div<-adonis2(ab.dist~Treatment, data=voc_transpose_rm3, permutations = perm, method="bray")
ab.div
# note: also can run models with rep as additive or interactive effect w/ treatment

# pairwise comparisons of permanova
ab.pair<-pairwise.adonis2(ab.dist~Treatment, data=voc_transpose_rm3, method="bray", strata="Rep")
ab.pair

# testing for homogeneity of dispersion among groups
# >0.05 meets assumption of adonis permanova
# i.e. adonis may give sig. p-value even if groups overlap because within-group data is heterogenous
# this test looks to see if group data is heterogenous; p>0.05 means it is not, and therefore 
# we can assume adonis results are "real" and not a result of heterogenous dispersion
dispersion<-betadisper(ab.dist, group=voc_transpose_rm3$Treatment)
dispersion
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T)


#### VOC Composition (rep) - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose_rm3[,2:429]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose_rm3$Treatment)
ab.div<-adonis2(ab.dist~Rep, data=voc_transpose_rm3, permutations = perm, method="bray")
ab.div

# pairwise comparisons of permanova
ab.pair<-pairwise.adonis2(ab.dist~Rep, data=voc_transpose_rm3, method="bray")
ab.pair

# testing for homogeneity of dispersion among groups
# >0.05 meets assumption of adonis permanova
# i.e. adonis may give sig. p-value even if groups overlap because within-group data is heterogenous
# this test looks to see if group data is heterogenous; p>0.05 means it is not, and therefore 
# we can assume adonis results are "real" and not a result of heterogenous dispersion
dispersion<-betadisper(ab.dist, group=voc_transpose_rm3$Rep)
dispersion
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T)



### VOC composition - ANOSIM (analysis of similarity) ###
ab = voc_transpose_rm3[,2:429]
mat_ab = as.matrix(ab)
ano = anosim(mat_ab, voc_transpose_rm3$Treatment, distance = "bray", permutations = 999)
ano


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


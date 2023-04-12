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
voc_meta_info2 <- voc_transpose_rm[,438, drop=FALSE]
voc_transpose_rm2 = subset(voc_transpose_rm, select = -c(Sample_ID,Unique_ID,Rep,Footprint,Subplot,Treatment,Notes,Weight_g,time_sampled))
# divide
voc_weighted_abun <- voc_transpose_rm2/voc_transpose_rm2[,429]
# remerging with meta info
voc_transpose_rm3 <- cbind(voc_sample_names,voc_weighted_abun,voc_meta_info,voc_meta_info2)
# removing indiv weight column 
voc_transpose_rm3 = subset(voc_transpose_rm3, select = -c(Weight_indiv_g))

# divide by hours sampled
# first remove compound name column & meta info columns
voc_sample_names_time <- voc_transpose_rm3[,1, drop=FALSE]
voc_meta_info_time <- voc_transpose_rm3[,430:436, drop=FALSE]
voc_transpose_rm_time = subset(voc_transpose_rm3, select = -c(Sample_ID,Unique_ID,Rep,Footprint,Subplot,Treatment,Notes,Weight_g))
# divide
voc_weighted_abun_time <- voc_transpose_rm_time/voc_transpose_rm_time[,429]
# remerging with meta info
voc_transpose_rm4 <- cbind(voc_sample_names_time,voc_weighted_abun_time,voc_meta_info_time)
# removing indiv weight column 
voc_transpose_rm4 = subset(voc_transpose_rm4, select = -c(time_sampled))
# so now when I want to do analyses for abundances/individual, I use the voc_transpose_rm4 dataframe



#### VOC Composition (treatment) - PERMANOVA ####
# make community matrix - extract columns with abundance information
ab = voc_transpose_rm4[,2:429]

# dissimilarity matrix
ab.dist<-vegdist(ab, method='bray')

# run permanova
set.seed(123)
perm <- how(nperm = 999, blocks=voc_transpose_rm4$Rep)
ab.div<-adonis2(ab.dist~Treatment, data=voc_transpose_rm4, permutations = perm, method="bray")
ab.div
# note: also can run models with rep as additive or interactive effect w/ treatment

# pairwise comparisons of permanova
ab.pair<-pairwise.adonis2(ab.dist~Treatment, data=voc_transpose_rm4, method="bray", strata="Rep")
ab.pair

# testing for homogeneity of dispersion among groups
# >0.05 meets assumption of adonis permanova
# i.e. adonis may give sig. p-value even if groups overlap because within-group data is heterogenous
# this test looks to see if group data is heterogenous; p>0.05 means it is not, and therefore 
# we can assume adonis results are "real" and not a result of heterogenous dispersion
dispersion<-betadisper(ab.dist, group=voc_transpose_rm4$Treatment)
dispersion
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T, label=F)


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
# Data exploration on raw data
voc_transpose_rm4$rowsums <- rowSums(voc_transpose_rm4[2:429])
descdist(voc_transpose_rm4$rowsums, discrete = FALSE)
hist(voc_transpose_rm4$rowsums)
qqnorm(voc_transpose_rm4$rowsums)
shapiro.test(voc_transpose_rm4$rowsums)
# kinda right skewed, going to try a few transformations

# square root transformation
voc_transpose_rm4$sqrt_rowsums <- sqrt(voc_transpose_rm4$rowsums)
descdist(voc_transpose_rm4$sqrt_rowsums, discrete = FALSE)
hist(voc_transpose_rm4$sqrt_rowsums)
qqnorm(voc_transpose_rm4$sqrt_rowsums)
shapiro.test(voc_transpose_rm4$sqrt_rowsums)

# cubed root transformation
voc_transpose_rm4$cubed_rowsums <- (voc_transpose_rm4$rowsums)^(1/3)
descdist(voc_transpose_rm4$cubed_rowsums, discrete = FALSE)
hist(voc_transpose_rm4$cubed_rowsums)
qqnorm(voc_transpose_rm4$cubed_rowsums)
shapiro.test(voc_transpose_rm4$cubed_rowsums)
# none of these are great, come back to this 

# comparison with other models
m1 <- lmer(cubed_rowsums ~ Treatment + (1|Rep), data = voc_transpose_rm4, REML=FALSE)
summary(m1)
hist(resid(m1))
shapiro.test(resid(m1))
# emmeans is a lot different than summary (?) so re-leveling to get pairwise comparisons
voc_transpose_rm4 <- within(voc_transpose_rm4, Treatment <- relevel(factor(Treatment), ref = "Warmed_Drought"))


# specific compound test
m2 <- lmer(.beta..Myrcene ~ Treatment + (1|Rep), data = voc_transpose_rm5, REML=FALSE)
anova(m2)
summary(m2)
# pairwise comparisons
voc_transpose_rm5 <- within(voc_transpose_rm5, Treatment <- relevel(factor(Treatment), ref = "Ambient_Control"))
voc_transpose_rm5 <- within(voc_transpose_rm5, Treatment <- relevel(factor(Treatment), ref = "Drought"))

m3 <- lmer(X4.Hexen.1.ol..acetate ~ Treatment + (1|Rep), data = voc_transpose_rm5, REML=FALSE)
anova(m3)
summary(m3)
voc_transpose_rm5 <- within(voc_transpose_rm5, Treatment <- relevel(factor(Treatment), ref = "Warmed"))
voc_transpose_rm5 <- within(voc_transpose_rm5, Treatment <- relevel(factor(Treatment), ref = "Drought"))


#### VOC Abundance - ANOVA ####
# overall treatment differences
anova1 <- aov(rowsums~Treatment, data = voc_transpose_rm4)
summary(anova1)
TukeyHSD(anova1)

# specific compound differences
for (i in 2:429){
  column <- names(voc_transpose_rm4[i])
  anova <- broom::tidy(aov(voc_transpose_rm4[,i] ~ Treatment, data = voc_transpose_rm4))
  
  # only want aov with P < 0.05 printed
  if(anova$p.value[1] < 0.05) {
    
    print(column)
    print(anova)
  }
}



# pairwise comparisons
for (i in 2:429){
  column <- names(voc_transpose_rm4[i])
  anova <- aov(voc_transpose_rm4[,i] ~ Treatment, data = voc_transpose_rm4)
  tukey <- TukeyHSD(anova)
  
  # only want tukey with P < 0.05 printed
  if(any(tukey$Treatment[, "p adj"] < 0.05)) {
    
    print(column)
    print(setNames(tukey, column))
    
  }
}


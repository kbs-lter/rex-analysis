# TITLE:          REX: 2022 Solidago leaves analysis
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_VOC L0 folder
# DATA OUTPUT:    Calculating total biomass of leaves measured for VOCs in each REX treatment
# PROJECT:        REX
# DATE:           July 2022



# Jan 2023 notes:
# also calculate biomass per rep, per treatment
# use this on figures to add in error bars for each treatment


### Data importing & checking ###
# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Read in data
length_weight <- read.csv(file.path(dir, "T7_warmx_VOC/L0/REX_T7_VOC_2022_Solidago_length_weight_L0.csv"))
voc_leaves <- read.csv(file.path(dir, "T7_warmx_VOC/L0/REX_T7_VOC_2022_Solidago_leaves_L0.csv"))

# checking spreadsheets
unique(length_weight$Treatment)
unique(voc_leaves$Treatment)

# removing plants from voc_leaves dataframe
# removing these because I don't use their VOC data, so I need to make sure I am not calculating more biomass than what was measured for VOCs
voc_leaves <- voc_leaves %>%
  filter(!(Rep == 5 & Treatment == "Warmed" & Plant_Number == 1)) %>%
  filter(!(Rep == 4 & Treatment == "Ambient" & Plant_Number == 1)) %>%
  filter(!(Rep == 4 & Treatment == "Warmed" & Plant_Number == 1)) %>%
  filter(!(Rep == 4 & Treatment == "Warmed_Drought" & Plant_Number == 1)) %>%
  filter(!(Rep == 3 & Treatment == "Warmed_Drought" & Plant_Number == 1)) %>%
  filter(!(Rep == 3 & Treatment == "Irrigated" & Plant_Number == 1)) %>%
  filter(!(Rep == 2 & Treatment == "Ambient" & Plant_Number == 1)) %>%
  filter(!(Rep == 1))


### Getting regression btwn leaf length and weight for each treatment ###
# making dataframes for each treatment in regression dataframe
warmed <- length_weight %>%
  filter(Treatment == "Warmed")
ambient <- length_weight %>%
  filter(Treatment == "Ambient")
irr <- length_weight %>%
  filter(Treatment == "Irrigated")
drought <- length_weight %>%
  filter(Treatment == "Drought")
warm_drought <- length_weight %>%
  filter(Treatment == "Warmed_Drought")

# regression btwn length and weight for leaves in each treatment
w_mod <- lm(Weight_g ~ Length_cm, data = warmed)
summary(w_mod)$coef
plot(Weight_g ~ Length_cm, data=warmed)
abline(w_mod)
a_mod <- lm(Weight_g ~ Length_cm, data = ambient)
summary(a_mod)$coef
plot(Weight_g ~ Length_cm, data=ambient)
abline(a_mod)
i_mod <- lm(Weight_g ~ Length_cm, data = irr)
summary(i_mod)$coef
plot(Weight_g ~ Length_cm, data=irr)
abline(i_mod)
d_mod <- lm(Weight_g ~ Length_cm, data = drought)
summary(d_mod)$coef
plot(Weight_g ~ Length_cm, data=drought)
abline(d_mod)
wd_mod <- lm(Weight_g ~ Length_cm, data = warm_drought)
summary(wd_mod)$coef
plot(Weight_g ~ Length_cm, data=warm_drought)
abline(wd_mod)

# making dataframe for each treatment & rep in voc_leaves dataframe
trt_rep <- function(df, trt, rep){
  df2 <- df %>%
    filter(Treatment == trt & Rep == rep) %>%
    select(-Treatment, -Rep, -Plant_Number)
  return(df2)
}
w_2 <- trt_rep(voc_leaves,"Warmed",2)
w_3 <- trt_rep(voc_leaves,"Warmed",3)
w_4 <- trt_rep(voc_leaves,"Warmed",4)
w_5 <- trt_rep(voc_leaves,"Warmed",5)

a_2 <- trt_rep(voc_leaves,"Ambient",2)
a_3 <- trt_rep(voc_leaves,"Ambient",3)
a_4 <- trt_rep(voc_leaves,"Ambient",4)
a_5 <- trt_rep(voc_leaves,"Ambient",5)

ir_2 <- trt_rep(voc_leaves,"Irrigated",2)
ir_3 <- trt_rep(voc_leaves,"Irrigated",3)
ir_4 <- trt_rep(voc_leaves,"Irrigated",4)
ir_5 <- trt_rep(voc_leaves,"Irrigated",5)

d_2 <- trt_rep(voc_leaves,"Drought",2)
d_3 <- trt_rep(voc_leaves,"Drought",3)
d_4 <- trt_rep(voc_leaves,"Drought",4)
d_5 <- trt_rep(voc_leaves,"Drought",5)

wd_2 <- trt_rep(voc_leaves,"Warmed_Drought",2)
wd_3 <- trt_rep(voc_leaves,"Warmed_Drought",3)
wd_4 <- trt_rep(voc_leaves,"Warmed_Drought",4)
wd_5 <- trt_rep(voc_leaves,"Warmed_Drought",5)

# making a dataframe to store biomass measurements in
voc_biomass <- data.frame(Treatment = c("Warmed","Warmed","Warmed","Warmed",
                                        "Ambient","Ambient","Ambient","Ambient",
                                        "Irrigated","Irrigated","Irrigated","Irrigated",
                                        "Drought","Drought","Drought","Drought",
                                        "Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought"),
                          Rep = c(2,3,4,5,
                                  2,3,4,5,
                                  2,3,4,5,
                                  2,3,4,5,
                                  2,3,4,5),
                          Weight_g = c(NA))

# predicting biomass per treatment + rep from leaf lengths
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 2, sum(predict(w_mod, w_2)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 3, sum(predict(w_mod, w_3)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 4, sum(predict(w_mod, w_4)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 5, sum(predict(w_mod, w_5)), voc_biomass$Weight_g)

voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 2, sum(predict(a_mod, a_2)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 3, sum(predict(a_mod, a_3)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 4, sum(predict(a_mod, a_4)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 5, sum(predict(a_mod, a_5)), voc_biomass$Weight_g)

voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 2, sum(predict(i_mod, ir_2)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 3, sum(predict(i_mod, ir_3)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 4, sum(predict(i_mod, ir_4)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 5, sum(predict(i_mod, ir_5)), voc_biomass$Weight_g)

voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 2, sum(predict(d_mod, d_2)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 3, sum(predict(d_mod, d_3)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 4, sum(predict(d_mod, d_4)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 5, sum(predict(d_mod, d_5)), voc_biomass$Weight_g)

voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 2, sum(predict(wd_mod, wd_2)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 3, sum(predict(wd_mod, wd_3)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 4, sum(predict(wd_mod, wd_4)), voc_biomass$Weight_g)
voc_biomass$Weight_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 5, sum(predict(wd_mod, wd_5)), voc_biomass$Weight_g)

# making a column for per-plant biomass
# because I can't tie back individual plany's biomass to their VOC emissions, I'm going to
# divide total biomass per rep + treatment by the number of plants measured in that plot
# to get a proxy for per-individual biomass
voc_biomass$Weight_indiv_g <- NA

voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 2, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 3, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 4, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed" & voc_biomass$Rep == 5, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)

voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 2, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 3, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 4, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Ambient" & voc_biomass$Rep == 5, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)

voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 2, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 3, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 4, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Irrigated" & voc_biomass$Rep == 5, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)

voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 2, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 3, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 4, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Drought" & voc_biomass$Rep == 5, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)

voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 2, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 3, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 4, voc_biomass$Weight_g/4, voc_biomass$Weight_indiv_g)
voc_biomass$Weight_indiv_g <- ifelse(voc_biomass$Treatment == "Warmed_Drought" & voc_biomass$Rep == 5, voc_biomass$Weight_g/5, voc_biomass$Weight_indiv_g)

# save output
write.csv(voc_biomass, file.path(dir,"T7_warmx_VOC/L1/VOC_biomass_2022_L1.csv"), row.names=F)


# old code for total treatment biomass
warmed_voc <- voc_leaves %>%
  filter(Treatment == "Warmed") %>%
  select(-Treatment, -Rep, -Plant_Number)
ambient_voc <- voc_leaves %>%
  filter(Treatment == "Ambient") %>%
  select(-Treatment, -Rep, -Plant_Number)
irr_voc <- voc_leaves %>%
  filter(Treatment == "Irrigated") %>%
  select(-Treatment, -Rep, -Plant_Number)
drought_voc <- voc_leaves %>%
  filter(Treatment == "Drought") %>%
  select(-Treatment, -Rep, -Plant_Number)
warm_drought_voc <- voc_leaves %>%
  filter(Treatment == "Warmed_Drought") %>%
  select(-Treatment, -Rep, -Plant_Number)

# predicting biomass from leaf lengths
sum(predict(w_mod, warmed_voc)) # 65.61699 g sampled
sum(predict(a_mod, ambient_voc)) # 60.93225 g sampled
sum(predict(i_mod, irr_voc)) # 65.71546 g sampled
sum(predict(d_mod, drought_voc)) # 57.64733 g sampled
sum(predict(wd_mod, warm_drought_voc)) # 50.26868 g sampled



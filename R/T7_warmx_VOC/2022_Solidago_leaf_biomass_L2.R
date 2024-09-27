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
summary(w_mod)$coef # y = 0.0235x - 0.0654
plot(Weight_g ~ Length_cm, data=warmed)
abline(w_mod)
a_mod <- lm(Weight_g ~ Length_cm, data = ambient)
summary(a_mod)$coef # y = 0.0217x - 0.0589
plot(Weight_g ~ Length_cm, data=ambient)
abline(a_mod)
i_mod <- lm(Weight_g ~ Length_cm, data = irr)
summary(i_mod)$coef # y = 0.0222x - 0.0558
plot(Weight_g ~ Length_cm, data=irr)
abline(i_mod)
d_mod <- lm(Weight_g ~ Length_cm, data = drought)
summary(d_mod)$coef # y = 0.0192x - 0.0458
plot(Weight_g ~ Length_cm, data=drought)
abline(d_mod)
wd_mod <- lm(Weight_g ~ Length_cm, data = warm_drought)
summary(wd_mod)$coef # y = 0.0221x - 0.0653
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

# adding in a column for total # hours samples for each rep, so that in my plots and analyses
# I can divide the abundance data by weight and hours to get emissions/g/hr
voc_biomass$time_sampled <- NA
voc_biomass$time_sampled[voc_biomass$Rep == 1] <- 5
voc_biomass$time_sampled[voc_biomass$Rep == 2] <- 5
voc_biomass$time_sampled[voc_biomass$Rep == 3] <- 7
voc_biomass$time_sampled[voc_biomass$Rep == 4] <- 7
voc_biomass$time_sampled[voc_biomass$Rep == 5] <- 7

# save output
write.csv(voc_biomass, file.path(dir,"T7_warmx_VOC/L1/VOC_biomass_2022_L1.csv"), row.names=F)



# calculating plant weight variability per plot
# doing this to see how much weights vary between individuals in a plot
# if variance is low, this validates my methods of dividing VOC abundances by avg. weight per plot
# making dataframe for each treatment, rep, and plant ID in voc_leaves dataframe
trt_rep_id <- function(df, trt, rep, id){
  df2 <- df %>%
    filter(Treatment == trt & Rep == rep & Plant_Number == id) %>%
    select(-Treatment, -Rep, -Plant_Number)
  return(df2)
}
a_2_2 <- trt_rep_id(voc_leaves,"Ambient",2,2)
a_2_3 <- trt_rep_id(voc_leaves,"Ambient",2,3)
a_2_4 <- trt_rep_id(voc_leaves,"Ambient",2,4)
a_2_5 <- trt_rep_id(voc_leaves,"Ambient",2,5)

a_3_1 <- trt_rep_id(voc_leaves,"Ambient",3,1)
a_3_2 <- trt_rep_id(voc_leaves,"Ambient",3,2)
a_3_3 <- trt_rep_id(voc_leaves,"Ambient",3,3)
a_3_4 <- trt_rep_id(voc_leaves,"Ambient",3,4)
a_3_5 <- trt_rep_id(voc_leaves,"Ambient",3,5)

a_4_2 <- trt_rep_id(voc_leaves,"Ambient",4,2)
a_4_3 <- trt_rep_id(voc_leaves,"Ambient",4,3)
a_4_4 <- trt_rep_id(voc_leaves,"Ambient",4,4)
a_4_5 <- trt_rep_id(voc_leaves,"Ambient",4,5)

a_5_1 <- trt_rep_id(voc_leaves,"Ambient",5,1)
a_5_2 <- trt_rep_id(voc_leaves,"Ambient",5,2)
a_5_3 <- trt_rep_id(voc_leaves,"Ambient",5,3)
a_5_4 <- trt_rep_id(voc_leaves,"Ambient",5,4)
a_5_5 <- trt_rep_id(voc_leaves,"Ambient",5,5)

w_2_1 <- trt_rep_id(voc_leaves,"Warmed",2,1)
w_2_2 <- trt_rep_id(voc_leaves,"Warmed",2,2)
w_2_3 <- trt_rep_id(voc_leaves,"Warmed",2,3)
w_2_4 <- trt_rep_id(voc_leaves,"Warmed",2,4)
w_2_5 <- trt_rep_id(voc_leaves,"Warmed",2,5)

w_3_1 <- trt_rep_id(voc_leaves,"Warmed",3,1)
w_3_2 <- trt_rep_id(voc_leaves,"Warmed",3,2)
w_3_3 <- trt_rep_id(voc_leaves,"Warmed",3,3)
w_3_4 <- trt_rep_id(voc_leaves,"Warmed",3,4)
w_3_5 <- trt_rep_id(voc_leaves,"Warmed",3,5)

w_4_2 <- trt_rep_id(voc_leaves,"Warmed",4,2)
w_4_3 <- trt_rep_id(voc_leaves,"Warmed",4,3)
w_4_4 <- trt_rep_id(voc_leaves,"Warmed",4,4)
w_4_5 <- trt_rep_id(voc_leaves,"Warmed",4,5)

w_5_2 <- trt_rep_id(voc_leaves,"Warmed",5,2)
w_5_3 <- trt_rep_id(voc_leaves,"Warmed",5,3)
w_5_4 <- trt_rep_id(voc_leaves,"Warmed",5,4)
w_5_5 <- trt_rep_id(voc_leaves,"Warmed",5,5)

d_2_1 <- trt_rep_id(voc_leaves,"Drought",2,1)
d_2_2 <- trt_rep_id(voc_leaves,"Drought",2,2)
d_2_3 <- trt_rep_id(voc_leaves,"Drought",2,3)
d_2_4 <- trt_rep_id(voc_leaves,"Drought",2,4)
d_2_5 <- trt_rep_id(voc_leaves,"Drought",2,5)

d_3_1 <- trt_rep_id(voc_leaves,"Drought",3,1)
d_3_2 <- trt_rep_id(voc_leaves,"Drought",3,2)
d_3_3 <- trt_rep_id(voc_leaves,"Drought",3,3)
d_3_4 <- trt_rep_id(voc_leaves,"Drought",3,4)
d_3_5 <- trt_rep_id(voc_leaves,"Drought",3,5)

d_4_1 <- trt_rep_id(voc_leaves,"Drought",4,1)
d_4_2 <- trt_rep_id(voc_leaves,"Drought",4,2)
d_4_3 <- trt_rep_id(voc_leaves,"Drought",4,3)
d_4_4 <- trt_rep_id(voc_leaves,"Drought",4,4)
d_4_5 <- trt_rep_id(voc_leaves,"Drought",4,5)

d_5_1 <- trt_rep_id(voc_leaves,"Drought",5,1)
d_5_2 <- trt_rep_id(voc_leaves,"Drought",5,2)
d_5_3 <- trt_rep_id(voc_leaves,"Drought",5,3)
d_5_4 <- trt_rep_id(voc_leaves,"Drought",5,4)
d_5_5 <- trt_rep_id(voc_leaves,"Drought",5,5)

wd_2_1 <- trt_rep_id(voc_leaves,"Warmed_Drought",2,1)
wd_2_2 <- trt_rep_id(voc_leaves,"Warmed_Drought",2,2)
wd_2_3 <- trt_rep_id(voc_leaves,"Warmed_Drought",2,3)
wd_2_4 <- trt_rep_id(voc_leaves,"Warmed_Drought",2,4)
wd_2_5 <- trt_rep_id(voc_leaves,"Warmed_Drought",2,5)

wd_3_2 <- trt_rep_id(voc_leaves,"Warmed_Drought",3,2)
wd_3_3 <- trt_rep_id(voc_leaves,"Warmed_Drought",3,3)
wd_3_4 <- trt_rep_id(voc_leaves,"Warmed_Drought",3,4)
wd_3_5 <- trt_rep_id(voc_leaves,"Warmed_Drought",3,5)

wd_4_2 <- trt_rep_id(voc_leaves,"Warmed_Drought",4,2)
wd_4_3 <- trt_rep_id(voc_leaves,"Warmed_Drought",4,3)
wd_4_4 <- trt_rep_id(voc_leaves,"Warmed_Drought",4,4)
wd_4_5 <- trt_rep_id(voc_leaves,"Warmed_Drought",4,5)

wd_5_1 <- trt_rep_id(voc_leaves,"Warmed_Drought",5,1)
wd_5_2 <- trt_rep_id(voc_leaves,"Warmed_Drought",5,2)
wd_5_3 <- trt_rep_id(voc_leaves,"Warmed_Drought",5,3)
wd_5_4 <- trt_rep_id(voc_leaves,"Warmed_Drought",5,4)
wd_5_5 <- trt_rep_id(voc_leaves,"Warmed_Drought",5,5)

i_2_1 <- trt_rep_id(voc_leaves,"Irrigated",2,1)
i_2_2 <- trt_rep_id(voc_leaves,"Irrigated",2,2)
i_2_3 <- trt_rep_id(voc_leaves,"Irrigated",2,3)
i_2_4 <- trt_rep_id(voc_leaves,"Irrigated",2,4)
i_2_5 <- trt_rep_id(voc_leaves,"Irrigated",2,5)

i_3_2 <- trt_rep_id(voc_leaves,"Irrigated",3,2)
i_3_3 <- trt_rep_id(voc_leaves,"Irrigated",3,3)
i_3_4 <- trt_rep_id(voc_leaves,"Irrigated",3,4)
i_3_5 <- trt_rep_id(voc_leaves,"Irrigated",3,5)

i_4_1 <- trt_rep_id(voc_leaves,"Irrigated",4,1)
i_4_2 <- trt_rep_id(voc_leaves,"Irrigated",4,2)
i_4_3 <- trt_rep_id(voc_leaves,"Irrigated",4,3)
i_4_4 <- trt_rep_id(voc_leaves,"Irrigated",4,4)
i_4_5 <- trt_rep_id(voc_leaves,"Irrigated",4,5)

i_5_1 <- trt_rep_id(voc_leaves,"Irrigated",5,1)
i_5_2 <- trt_rep_id(voc_leaves,"Irrigated",5,2)
i_5_3 <- trt_rep_id(voc_leaves,"Irrigated",5,3)
i_5_4 <- trt_rep_id(voc_leaves,"Irrigated",5,4)
i_5_5 <- trt_rep_id(voc_leaves,"Irrigated",5,5)




# making a dataframe to store biomass measurements in
voc_biomass2 <- data.frame(Treatment = c("Ambient","Ambient","Ambient","Ambient",
                                        "Ambient","Ambient","Ambient","Ambient","Ambient",
                                        "Ambient","Ambient","Ambient","Ambient",
                                        "Ambient","Ambient","Ambient","Ambient","Ambient",
                                        "Warmed","Warmed","Warmed","Warmed","Warmed",
                                        "Warmed","Warmed","Warmed","Warmed","Warmed",
                                        "Warmed","Warmed","Warmed","Warmed",
                                        "Warmed","Warmed","Warmed","Warmed",
                                        "Drought","Drought","Drought","Drought","Drought",
                                        "Drought","Drought","Drought","Drought","Drought",
                                        "Drought","Drought","Drought","Drought","Drought",
                                        "Drought","Drought","Drought","Drought","Drought",
                                        "Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought",
                                        "Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought",
                                        "Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought",
                                        "Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought","Warmed_Drought",
                                        "Irrigated","Irrigated","Irrigated","Irrigated","Irrigated",
                                        "Irrigated","Irrigated","Irrigated","Irrigated",
                                        "Irrigated","Irrigated","Irrigated","Irrigated","Irrigated",
                                        "Irrigated","Irrigated","Irrigated","Irrigated","Irrigated"),
                          Rep = c(2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,5,
                                  2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,
                                  2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,
                                  2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,
                                  2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5),
                          Plant_ID = c(2,3,4,5,1,2,3,4,5,2,3,4,5,1,2,3,4,5,
                                       1,2,3,4,5,1,2,3,4,5,2,3,4,5,2,3,4,5,
                                       1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,
                                       1,2,3,4,5,2,3,4,5,2,3,4,5,1,2,3,4,5,
                                       1,2,3,4,5,2,3,4,5,1,2,3,4,5,1,2,3,4,5),
                          Weight_g = c(NA))

# predicting biomass per treatment + rep from leaf lengths
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 2, sum(predict(a_mod, a_2_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 3, sum(predict(a_mod, a_2_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 4, sum(predict(a_mod, a_2_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 5, sum(predict(a_mod, a_2_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 1, sum(predict(a_mod, a_3_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 2, sum(predict(a_mod, a_3_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 3, sum(predict(a_mod, a_3_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 4, sum(predict(a_mod, a_3_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 5, sum(predict(a_mod, a_3_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 2, sum(predict(a_mod, a_4_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 3, sum(predict(a_mod, a_4_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 4, sum(predict(a_mod, a_4_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 5, sum(predict(a_mod, a_4_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 1, sum(predict(a_mod, a_5_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 2, sum(predict(a_mod, a_5_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 3, sum(predict(a_mod, a_5_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 4, sum(predict(a_mod, a_5_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Ambient" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 5, sum(predict(a_mod, a_5_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 1, sum(predict(w_mod, w_2_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 2, sum(predict(w_mod, w_2_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 3, sum(predict(w_mod, w_2_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 4, sum(predict(w_mod, w_2_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 5, sum(predict(w_mod, w_2_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 1, sum(predict(w_mod, w_3_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 2, sum(predict(w_mod, w_3_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 3, sum(predict(w_mod, w_3_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 4, sum(predict(w_mod, w_3_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 5, sum(predict(w_mod, w_3_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 2, sum(predict(w_mod, w_4_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 3, sum(predict(w_mod, w_4_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 4, sum(predict(w_mod, w_4_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 5, sum(predict(w_mod, w_4_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 2, sum(predict(w_mod, w_5_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 3, sum(predict(w_mod, w_5_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 4, sum(predict(w_mod, w_5_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 5, sum(predict(w_mod, w_5_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 1, sum(predict(d_mod, d_2_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 2, sum(predict(d_mod, d_2_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 3, sum(predict(d_mod, d_2_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 4, sum(predict(d_mod, d_2_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 5, sum(predict(d_mod, d_2_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 1, sum(predict(d_mod, d_3_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 2, sum(predict(d_mod, d_3_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 3, sum(predict(d_mod, d_3_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 4, sum(predict(d_mod, d_3_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 5, sum(predict(d_mod, d_3_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 1, sum(predict(d_mod, d_4_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 2, sum(predict(d_mod, d_4_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 3, sum(predict(d_mod, d_4_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 4, sum(predict(d_mod, d_4_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 5, sum(predict(d_mod, d_4_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 1, sum(predict(d_mod, d_5_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 2, sum(predict(d_mod, d_5_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 3, sum(predict(d_mod, d_5_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 4, sum(predict(d_mod, d_5_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 5, sum(predict(d_mod, d_5_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 1, sum(predict(wd_mod, wd_2_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 2, sum(predict(wd_mod, wd_2_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 3, sum(predict(wd_mod, wd_2_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 4, sum(predict(wd_mod, wd_2_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 5, sum(predict(wd_mod, wd_2_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 2, sum(predict(wd_mod, wd_3_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 3, sum(predict(wd_mod, wd_3_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 4, sum(predict(wd_mod, wd_3_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 5, sum(predict(wd_mod, wd_3_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 2, sum(predict(wd_mod, wd_4_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 3, sum(predict(wd_mod, wd_4_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 4, sum(predict(wd_mod, wd_4_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 5, sum(predict(wd_mod, wd_4_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 1, sum(predict(wd_mod, wd_5_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 2, sum(predict(wd_mod, wd_5_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 3, sum(predict(wd_mod, wd_5_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 4, sum(predict(wd_mod, wd_5_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Warmed_Drought" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 5, sum(predict(wd_mod, wd_5_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 1, sum(predict(i_mod, i_2_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 2, sum(predict(i_mod, i_2_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 3, sum(predict(i_mod, i_2_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 4, sum(predict(i_mod, i_2_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 2 & voc_biomass2$Plant_ID == 5, sum(predict(i_mod, i_2_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 2, sum(predict(i_mod, i_3_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 3, sum(predict(i_mod, i_3_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 4, sum(predict(i_mod, i_3_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 3 & voc_biomass2$Plant_ID == 5, sum(predict(i_mod, i_3_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 1, sum(predict(i_mod, i_4_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 2, sum(predict(i_mod, i_4_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 3, sum(predict(i_mod, i_4_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 4, sum(predict(i_mod, i_4_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 4 & voc_biomass2$Plant_ID == 5, sum(predict(i_mod, i_4_5)), voc_biomass2$Weight_g)

voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 1, sum(predict(i_mod, i_5_1)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 2, sum(predict(i_mod, i_5_2)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 3, sum(predict(i_mod, i_5_3)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 4, sum(predict(i_mod, i_5_4)), voc_biomass2$Weight_g)
voc_biomass2$Weight_g <- ifelse(voc_biomass2$Treatment == "Irrigated" & voc_biomass2$Rep == 5 & voc_biomass2$Plant_ID == 5, sum(predict(i_mod, i_5_5)), voc_biomass2$Weight_g)

# calculating mean +/- SD for each plot
voc_biomass_var <- voc_biomass2 %>%
  group_by(Rep,Treatment) %>%
  summarize(avg_weight = mean(Weight_g),
            sd_weight = sd(Weight_g))





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



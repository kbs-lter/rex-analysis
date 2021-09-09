# TITLE:          REX: Soca VOC plots
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

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021predrought_L1.csv"))

# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:269]

# turn abundance data frame into a matrix
mat_ab = as.matrix(ab)

# generate nmds plot
set.seed(123)
nmds = metaMDS(mat_ab, distance = "bray")
nmds
plot(nmds)

# extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

# add columns to data frame 
data.scores$Compound = voc_transpose$Compound
data.scores$Treatment = voc_transpose$Treatment
head(data.scores)

# plot
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = Treatment))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")
  #scale_colour_manual(values = c("#009E73", "#E69F00")) 

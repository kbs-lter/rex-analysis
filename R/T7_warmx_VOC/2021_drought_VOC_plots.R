# TITLE:          REX: Drought VOC plots
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
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021drought_L1.csv"))

# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:292]

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
data.scores$Group_treat = voc_transpose$Group_treat
data.scores$Rep = as.character(voc_transpose$Rep)
head(data.scores)
str(data.scores)

# plot
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = Treatment, shape = Rep))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")

# plot with hulls - groups for each treatment
W <- data.scores[data.scores$Treatment == "Warmed", ][chull(data.scores[data.scores$Treatment == 
                                                                          "Warmed", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
A <- data.scores[data.scores$Treatment == "Ambient", ][chull(data.scores[data.scores$Treatment == 
                                                                           "Ambient", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
D <- data.scores[data.scores$Treatment == "Drought", ][chull(data.scores[data.scores$Treatment == 
                                                                           "Drought", c("NMDS1", "NMDS2")]), ]  # hull values for grp C
WD <- data.scores[data.scores$Treatment == "Warmed_Drought", ][chull(data.scores[data.scores$Treatment == 
                                                                                   "Warmed_Drought", c("NMDS1", "NMDS2")]), ]  # hull values for grp D
I <- data.scores[data.scores$Treatment == "Irrigated", ][chull(data.scores[data.scores$Treatment == 
                                                                             "Irrigated", c("NMDS1", "NMDS2")]), ]  # hull values for grp E
hull.data <- rbind(W, A, D, WD, I)  #combine groups
hull.data

ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Treatment,group=Treatment),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Treatment),size=4) + # add the point markers
  coord_equal() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")

# plot with hulls - groups for each rep
Rep1 <- data.scores[data.scores$Rep == "1", ][chull(data.scores[data.scores$Rep == 
                                                                "1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
Rep3 <- data.scores[data.scores$Rep == "3", ][chull(data.scores[data.scores$Rep == 
                                                                "3", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
Rep5 <- data.scores[data.scores$Rep == "5", ][chull(data.scores[data.scores$Rep == 
                                                                "5", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

hull.data2 <- rbind(Rep1, Rep3, Rep5)  #combine groups
hull.data2

ggplot() + 
  geom_polygon(data=hull.data2,aes(x=NMDS1,y=NMDS2,fill=Rep,group=Rep),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Rep),size=4) + # add the point markers
  coord_equal() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")


# plot with grouped treatments for warming and control
# plot with hulls - groups for each treatment
W2 <- data.scores[data.scores$Group_treat == "Warmed", ][chull(data.scores[data.scores$Group_treat == 
                                                                 "Warmed", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
A2 <- data.scores[data.scores$Group_treat == "Control", ][chull(data.scores[data.scores$Group_treat == 
                                                                  "Control", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
D2 <- data.scores[data.scores$Group_treat == "Drought", ][chull(data.scores[data.scores$Group_treat == 
                                                                           "Drought", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

hull.data3 <- rbind(W2, A2, D2)  #combine groups
hull.data3

ggplot() + 
  geom_polygon(data=hull.data3,aes(x=NMDS1,y=NMDS2,fill=Group_treat,group=Group_treat),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Group_treat),size=4) + # add the point markers
  coord_equal() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")

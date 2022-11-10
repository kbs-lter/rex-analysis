# TITLE:          REX: 2022 VOC plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_VOC L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           aug 2022

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(vegan)
library(plotrix)
library(plotly)
library(broom)
library(Rtsne)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_named_VOC_2022_L1.csv"))


#### NMDS ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:783]

# turn abundance data frame into a matrix
mat_ab = as.matrix(ab)

# generate nmds plot
set.seed(1)
nmds = metaMDS(mat_ab, distance = "bray",k = 2)
nmds
plot(nmds)

# extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

# add columns to data frame 
data.scores$Sample_ID = voc_transpose$Sample_ID
data.scores$Treatment = voc_transpose$Treatment
data.scores$Rep = as.character(voc_transpose$Rep)
head(data.scores)
str(data.scores)

# plot
png("rep_nmds.png", units="in", width=6, height=5, res=300)
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = Rep))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  stat_ellipse(aes(fill=Rep), alpha=.2,type='t',size =1, geom="polygon")+
  labs(x = "NMDS1", colour = "Rep", y = "NMDS2")
dev.off()

#### PCoA - Treatment ####
ab = voc_transpose[,2:783]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
png("climate_pcoa.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, color=grps),size=4,shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                     values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                    values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()


#### PCoA - Grouped Treatments ####
voc_transpose2 <- voc_transpose
#voc_transpose2 <- voc_transpose2[!(voc_transpose2$Treatment == "Irrigated_Control"),] 
#voc_transpose2 <- voc_transpose2[!(voc_transpose2$Treatment == "Ambient_Control"),] 
#voc_transpose2 <- voc_transpose2[!(voc_transpose2$Treatment == "Warmed_Drought"),] 
#voc_transpose2 <- voc_transpose2[!(voc_transpose2$Treatment == "Warmed"),] 
#voc_transpose2 <- voc_transpose2[!(voc_transpose2$Treatment == "Drought"),] 
voc_transpose2$Treatment <- gsub('Warmed_Drought', 'Drought', voc_transpose2$Treatment)
voc_transpose2$Treatment <- gsub('Ambient_Control', 'Not_Drought', voc_transpose2$Treatment)
voc_transpose2$Treatment <- gsub('Irrigated_Control', 'Not_Drought', voc_transpose2$Treatment)
voc_transpose2$Treatment <- gsub('Warmed', 'Not_Drought', voc_transpose2$Treatment)
ab = voc_transpose2[,2:783]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose2$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
png("climate_pcoa_grouped.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, color=grps),size=4,shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  #scale_color_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
  #                   values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  #scale_fill_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
  #                  values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()


#### PCoA - Rep ####
ab = voc_transpose[,2:783]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Rep)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
png("rep_pcoa.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, color=grps),size=4,shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("1", "2", "3", "4", "5"),
                     values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1","#d7191c"))+
  scale_fill_manual(labels=c("1", "2", "3", "4", "5"),
                    values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Field Rep", fill="Field Rep") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))
dev.off()


#### ABUNDANCE ####
# calculating total abundance for each sample
voc_transpose$rowsums <- rowSums(voc_transpose[2:1455])

# calculating total abundance for each treatment
voc_transpose_sum <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(abun = sum(rowsums))
voc_transpose_sum2 <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(abun = mean(rowsums),
            se = std.error(rowsums))

# manually making a dataframe that takes abundance sums divided by plant biomass for each treatment
voc_transpose_w <- data.frame(Treatment = c("Ambient","Irrigated","Drought","Warmed","Warmed_Drought"),
                  weighted_abun = c(4414484,2202304,1920915,3319851,6272371))
# total abundance from voc_transpose_sum, plant biomass from the 2022_Solidago_leaf_biomass_L2.R script
# ambient 248766556/56.35235
# irrigated 144725416/65.71546
# drought 105300748/54.81802
# warmed 217838644/65.61699
# warmed drought 315303824/50.26868

# total abundances
level_order <- c('Ambient_Control', 'Irrigated_Control', 'Drought', "Warmed", "Warmed_Drought")
png("climate_ab.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum2, aes(x = factor(Treatment, level = level_order), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Average VOC Abundance")
dev.off()

# abundances for specific compounds, found in the analyses script
# alpha farnesene
voc_transpose_cmpd <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, X.Z.Z...alpha..Farnesene) %>%
  group_by(Treatment) %>%
  summarize(abun = sum(X.Z.Z...alpha..Farnesene))
voc_transpose_cmpd_a <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, X.Z.Z...alpha..Farnesene) %>%
  group_by(Treatment) %>%
  summarize(abun = mean(X.Z.Z...alpha..Farnesene),
            se = std.error(X.Z.Z...alpha..Farnesene))
# manually making a dataframe that takes abundance sums divided by plant biomass for each treatment
voc_transpose_cmpd_w <- data.frame(Treatment = c("Ambient","Irrigated","Drought","Warmed","Warmed_Drought"),
                              weighted_abun = c(30471.56,3146.155,5789.228,19820.35,1954.955))
# total abundance from voc_transpose_sum, plant biomass from the 2022_Solidago_leaf_biomass_L2.R script
# ambient 1717144/56.35235
# irrigated 206751/65.71546
# drought 317354/54.81802
# warmed 1300552/65.61699
# warmed drought 98273/50.26868
png("climate_farn.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose_cmpd_a, aes(x = factor(Treatment, level = level_order), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Average VOC Abundance")
dev.off()

# endo Borneol
voc_transpose_cmpd2 <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, endo.Borneol) %>%
  group_by(Treatment) %>%
  summarize(abun = sum(endo.Borneol))
voc_transpose_cmpd2_a <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, endo.Borneol) %>%
  group_by(Treatment) %>%
  summarize(abun = mean(endo.Borneol),
            se = std.error(endo.Borneol))
voc_transpose_cmpd2_w <- data.frame(Treatment = c("Ambient","Irrigated","Drought","Warmed","Warmed_Drought"),
                                   weighted_abun = c(6189.467,209.2354,745.1382,9477.088,688.2019))
# total abundance from voc_transpose_sum, plant biomass from the 2022_Solidago_leaf_biomass_L2.R script
# ambient 348791/56.35235
# irrigated 13750/65.71546
# drought 40847/54.81802
# warmed 621858/65.61699
# warmed drought 34595/50.26868
png("climate_endo.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose_cmpd2_a, aes(x = factor(Treatment, level = level_order), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Average VOC Abundance")
dev.off()



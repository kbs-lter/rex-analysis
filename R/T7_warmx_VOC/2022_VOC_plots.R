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
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2022_L1.csv"))
voc_transpose <- voc_transpose %>%
  filter(!(Treatment == "Bag"))


#### NMDS ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:3011]

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
  stat_ellipse(aes(fill=Treatment), alpha=.2,type='t',size =1, geom="polygon")+
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")


#### PCoA ####
ab = voc_transpose[,2:3011]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
#png("drought_comp.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, color=grps),size=4,shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                     values=c('#2c7bb6','#abd9e9',"#ffffbf","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                    values=c('#2c7bb6','#abd9e9',"#ffffbf","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
#dev.off()


#### ABUNDANCE ####
voc_transpose$rowsums <- rowSums(voc_transpose[2:3011])

voc_transpose_sum <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(abun = sum(rowsums),
            se = std.error(rowsums))

voc_transpose_sum_rep <- voc_transpose %>%
  group_by(Rep) %>%
  summarize(abun = sum(rowsums),
            se = std.error(rowsums))

level_order <- c('Ambient_Control', 'Irrigated_Control', 'Warmed', "Warmed_Drought", "Drought")
png("drought_ab.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum, aes(x = factor(Treatment, level = level_order), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "goldenrod") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Treatment")
dev.off()
png("drought_tot_ab.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose, aes(x = factor(Treatment, level = level_order), y = rowsums)) + 
  geom_boxplot(color = 'black', fill = "cornflowerblue", outlier.shape=NA) +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Abundance", fill = "Treatment") +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15,face="bold")) +
  scale_x_discrete(limits = c("Ambient_Control", "Irrigated_Control", "Warmed", "Warmed_Drought", "Drought"),
                   labels=c("Ambient_Control" = "Ambient",
                            "Drought" = "Drought",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed" = "Warmed",
                            "Warmed_Drought" = "Warmed + Drought"),
                   guide = guide_axis(n.dodge=2)) +
  ylim(0e+00,3e+07)
dev.off()


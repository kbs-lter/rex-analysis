# TITLE:          REX: Pre-drought VOC plots
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
library(plotrix)
library(ggfortify)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021predrought_L1.csv"))



#### NMDS ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:300]

# turn abundance data frame into a matrix
mat_ab = as.matrix(ab)

# generate nmds plot
set.seed(1)
nmds = metaMDS(mat_ab, distance = "bray")
nmds
plot(nmds)

# the plot below is nightmarish, would be nice with fewer compounds/samples
ordiplot(nmds,type="n")
orditorp(nmds,display="species",col="red",air=0.01)
orditorp(nmds,display="sites",cex=1.25,air=0.01)

# Shepard plot: Large scatter around the line suggests that original dissimilarities
# are not well preserved in the reduced number of dimensions
stressplot(nmds)

# ggplot version
# extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
# add columns to data frame 
data.scores$Compound = voc_transpose$Compound
data.scores$Treatment = voc_transpose$Treatment
data.scores$Rep = as.character(voc_transpose$Rep)
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
  stat_ellipse(aes(fill=Treatment), alpha=.2,type='t',size =0.5, level=0.7, geom="polygon")+
  scale_x_discrete(labels=c("Ambient" = "Ambient",
                            "Drought" = "Drought",
                            "Irrigated" = "Irrigated",
                            "Warmed" = "Warmed",
                            "WarmedDrought" = "Warmed + Drought")) +
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")
  #scale_colour_manual(values = c("#009E73", "#E69F00")) 

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



#### PCoA ####
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
dispersion
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
png("predrought_comp.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, color=grps),size=4,shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                     values=c('#2c7bb6','#abd9e9',"#ffffbf","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
                    values=c('#2c7bb6','#abd9e9',"#ffffbf","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()




#### ABUNDANCE ####
voc_transpose$rowsums <- rowSums(voc_transpose[2:300])

voc_transpose_sum <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(abun = sum(rowsums),
            se = std.error(rowsums))

level_order <- c('Ambient', 'Irrigated', 'Warmed', "WarmedDrought", "Drought")
png("predrought_comp.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum, aes(x = factor(Treatment, level = level_order), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "goldenrod") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Treatment")
dev.off()
png("predrought_tot_ab.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose, aes(x = factor(Treatment, level = level_order), y = rowsums)) + 
  geom_boxplot(color = 'black', fill = "cornflowerblue") +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Relative Abundance", fill = "Treatment") +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15,face="bold")) +
  scale_x_discrete(limits = c("Ambient", "Irrigated", "Warmed", "WarmedDrought", "Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Drought" = "Drought",
                            "Irrigated" = "Irrigated",
                            "Warmed" = "Warmed",
                            "WarmedDrought" = "Warmed + Drought"),
                   guide = guide_axis(n.dodge=2))
dev.off()


#### PCA ####
voc_transpose2 <- voc_transpose
row.names(voc_transpose2) <- paste(voc_transpose2$Treatment, row.names(voc_transpose2), sep="_") 
voc_transpose2$Treatment <- NULL

pca_res <- prcomp(ab, scale. = F)

autoplot(pca_res, data=voc_transpose, color="Treatment")
plot(pca_res$x[,1], pca_res$x[,2]) # same plot as above

df_out <- as.data.frame(pca_res$x)
df_out$group <- sapply( strsplit(as.character(row.names(voc_transpose2)), "_"), "[[", 1 )
head(df_out)

ggplot(df_out, aes(x=PC1, y=PC2, color=group)) +
  geom_point() +
  theme_classic()

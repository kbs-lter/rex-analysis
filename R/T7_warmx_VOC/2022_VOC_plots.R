# TITLE:          REX: 2022 VOC plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_warmx_VOC L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           Aug 2022; edited Jan 2023

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(vegan)
library(plotrix)
library(plotly)
library(broom)
library(Rtsne)
library(gridExtra)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_named_VOC_2022_L1.csv"))



#### NMDS ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:429]

# turn abundance data frame into a matrix
mat_ab = as.matrix(ab)

# generate nmds plot
set.seed(1)
nmds = metaMDS(mat_ab, distance = "bray",k = 2)
nmds
plot(nmds)

# extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds)) # broken?

# add columns to data frame 
data.scores$Sample_ID = voc_transpose_rm$Sample_ID
data.scores$Treatment = voc_transpose_rm$Treatment
data.scores$Rep = as.character(voc_transpose_rm$Rep)
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
# selecting compound columns
ab = voc_transpose[,2:429]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(group=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
#seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3]) #wrong
seg.data<-merge(vectors[,1:3],centroids[,1:3], by="group")
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# specify order of treatments
seg.data$group <- factor(seg.data$group, levels = c("Ambient_Control","Irrigated_Control","Warmed","Drought","Warmed_Drought"))
levels(seg.data$group)
# make figure
png("climate_pcoa.png", units="in", width=6.5, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.8, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=3,shape=16) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) + 
  scale_color_manual(labels=c("Ambient (A)", "Irrigated (I)", "Warmed (W)", "Drought (D)", "Warmed + Drought \n (WD)"),
                     values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient (A)", "Irrigated (I)", "Warmed (W)", "Drought (D)", "Warmed + Drought \n (WD)"),
                    values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()


# w/ no irrigated
voc_transpose_rm5 <- voc_transpose %>%
  filter(!(Treatment == "Irrigated_Control"))
# selecting compound columns
ab = voc_transpose_rm5[,2:429]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose_rm5$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(group=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))
seg.data<-merge(vectors[,1:3],centroids[,1:3], by="group")
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# specify order of treatments
seg.data$group <- factor(seg.data$group, levels = c("Ambient_Control","Warmed","Drought","Warmed_Drought"))
levels(seg.data$group)
# make figure
png("climate_pcoa_no_ir.png", units="in", width=6.5, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.8, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=3,shape=16) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) + 
  scale_color_manual(labels=c("Ambient (A)", "Warmed (W)", "Drought (D)", "Warmed + Drought \n (WD)"),
                     values=c('#2c7bb6',"khaki1","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient (A)", "Warmed (W)", "Drought (D)", "Warmed + Drought \n (WD)"),
                    values=c('#2c7bb6',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()

# rep w/ no irrigated
ab = voc_transpose_rm5[,2:429]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose_rm5$Rep)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(group=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))
seg.data<-merge(vectors[,1:3],centroids[,1:3], by="group")
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# specify order of treatments
levels(seg.data$group)
# make figure
png("rep_pcoa_no_ir.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) +   
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("2", "3", "4", "5"),
                     values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1"))+
  scale_fill_manual(labels=c("2", "3", "4", "5"),
                    values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1"))+
  labs(x="PCoA 1",y="PCoA 2", color="Field Rep", fill="Field Rep") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))
dev.off()

# why does the permanova show sig. results when there's so much overlap?
# https://chrischizinski.github.io/rstats/adonis/
dispersion
# create the convex hulls of the outermost points
grp1.hull<-seg.data[seg.data$group=="Ambient_Control",1:3][chull(seg.data[seg.data$group=="Ambient_Control",2:3]),]
grp2.hull<-seg.data[seg.data$group=="Drought",1:3][chull(seg.data[seg.data$group=="Drought",2:3]),]
grp3.hull<-seg.data[seg.data$group=="Irrigated_Control",1:3][chull(seg.data[seg.data$group=="Irrigated_Control",2:3]),]
grp4.hull<-seg.data[seg.data$group=="Warmed",1:3][chull(seg.data[seg.data$group=="Warmed",2:3]),]
grp5.hull<-seg.data[seg.data$group=="Warmed_Drought",1:3][chull(seg.data[seg.data$group=="Warmed_Drought",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull,grp4.hull,grp5.hull)

# plot data for each treatment
panel.a<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="Ambient_Control",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group=="Ambient_Control",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[1,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) + 
  geom_point(data=seg.data[seg.data$group=="Ambient_Control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16) +
  labs(title="Ambient_Control",x="",y="") +
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme(legend.position="none")

panel.b<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="Drought",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group=="Drought",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[seg.data$group=="Drought",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17) +
  labs(title="Drought",x="",y="") +
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme(legend.position="none")

panel.c<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="Irrigated_Control",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group=="Irrigated_Control",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group=="Irrigated_Control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15) +
  labs(title="Irrigated Control",x="",y="") +
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme(legend.position="none")

panel.d<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="Warmed",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group=="Warmed",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[4,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=14) + 
  geom_point(data=seg.data[seg.data$group=="Warmed",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=14) +
  labs(title="Warmed",x="",y="") +
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme(legend.position="none")

panel.e<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="Warmed_Drought",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group=="Warmed_Drought",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[5,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=13) + 
  geom_point(data=seg.data[seg.data$group=="Warmed_Drought",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=13) +
  labs(title="Warmed Drought",x="",y="") +
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme(legend.position="none")

# all treatments (most informative)
panel.f<-ggplot() + 
  geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=group),size=4,colour="red") + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=group),size=2) +
  labs(title="All",x="",y="")
  #coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  #theme(legend.position="none")

grid.arrange(panel.a,panel.b,panel.c,panel.d,panel.e,panel.f,nrow=1)

# so what these figures show, especially panel.f, is the variance of the data for each treatment
# as well as the centroids of the data
# so we can see that drought and warmed drought's centroids are further from the other treatments,
# hence their difference seen in permanova. The variance of each treatment seems similar



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
ab = voc_transpose2[,2:429]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose2$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(group=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
#seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
seg.data<-merge(vectors[,1:3],centroids[,1:3], by="group")
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# specify order of treatments
levels(seg.data$group)
# make figure
png("climate_pcoa_grouped.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) +  
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  #scale_color_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
  #                   values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  #scale_fill_manual(labels=c("Ambient (A)", "Drought (D)", "Irrigated (I)", "Warmed (W)", "Warmed + Drought \n (WD)"),
  #                  values=c('#2c7bb6','#abd9e9',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  theme_classic()
dev.off()



#### PCoA - Rep ####
ab = voc_transpose[,2:429]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Rep)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(group=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
#seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
seg.data<-merge(vectors[,1:3],centroids[,1:3], by="group")
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# specify order of treatments
levels(seg.data$group)
# make figure
png("rep_pcoa.png", units="in", width=6, height=5, res=300)
ggplot() + 
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.7, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) +   
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=2.5,shape=16) +
  scale_color_manual(labels=c("2", "3", "4", "5"),
                     values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1"))+
  scale_fill_manual(labels=c("2", "3", "4", "5"),
                    values=c('#2c7bb6','#abd9e9',"#fdae61","khaki1"))+
  labs(x="PCoA 1",y="PCoA 2", color="Field Rep", fill="Field Rep") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))
dev.off()



#### ABUNDANCE ####
# calculating total abundance for each sample
voc_transpose$rowsums <- rowSums(voc_transpose[2:429])
voc_transpose_rm5$rowsums <- rowSums(voc_transpose_rm5[2:429])

# calculating total abundance for each treatment $ rep
#voc_transpose_sum <- voc_transpose %>%
#  group_by(Treatment, Rep) %>%
#  summarize(abun = sum(rowsums))
# avg abundance per treaatment biomass
voc_transpose_sum2 <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(avg_abun = mean(rowsums),
            se = std.error(rowsums))
voc_transpose_sum3 <- voc_transpose_rm5 %>%
  group_by(Treatment) %>%
  summarize(avg_abun = mean(rowsums),
            se = std.error(rowsums))
# plot - abundance w/o considering biomass
level_order2 <- c('Ambient_Control', 'Irrigated_Control', 'Drought', "Warmed", "Warmed_Drought")
png("climate_ab_no_ir.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum3, aes(x = factor(Treatment, level = level_order2), y = avg_abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  geom_errorbar(aes(ymin = avg_abun - se, ymax = avg_abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Average VOC Abundance \n (peak area/g/hr)")
dev.off()



### abundances for specific compounds, found in the analyses script ###
level_order2 <- c('Ambient_Control', 'Irrigated_Control', 'Drought', "Warmed", "Warmed_Drought")
level_order3 <- c('Ambient_Control', 'Warmed','Drought', "Warmed_Drought")

# o xylene
voc_transpose_xylene <- voc_transpose_rm5 %>%
  dplyr::select(Sample_ID, Treatment, Rep, o.Xylene)
# taking average per treatment
voc_transpose_xylene_a <- voc_transpose_xylene %>%
  group_by(Treatment) %>%
  summarize(abun = mean(o.Xylene),
            se = std.error(o.Xylene))

# Hepten 2 one 6 methyl
voc_transpose_hepten <- voc_transpose_rm5 %>%
  dplyr::select(Sample_ID, Treatment, Rep, X5.Hepten.2.one..6.methyl.)
# taking average per treatment
voc_transpose_hepten_a <- voc_transpose_hepten %>%
  group_by(Treatment) %>%
  summarize(abun = mean(X5.Hepten.2.one..6.methyl.),
            se = std.error(X5.Hepten.2.one..6.methyl.))

# beta Myrcene
voc_transpose_myrcene <- voc_transpose_rm5 %>%
  dplyr::select(Sample_ID, Treatment, Rep, .beta..Myrcene)
# taking average per treatment
voc_transpose_myrcene_a <- voc_transpose_myrcene %>%
  group_by(Treatment) %>%
  summarize(abun = mean(.beta..Myrcene),
            se = std.error(.beta..Myrcene))

# Hexen 1 ol acetate
voc_transpose_hexen <- voc_transpose_rm5 %>%
  dplyr::select(Sample_ID, Treatment, Rep, X4.Hexen.1.ol..acetate)
# taking average per treatment
voc_transpose_hexen_a <- voc_transpose_hexen %>%
  group_by(Treatment) %>%
  summarize(abun = mean(X4.Hexen.1.ol..acetate),
            se = std.error(X4.Hexen.1.ol..acetate))

# plot
png("climate_hexen.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose_hexen_a, aes(x = factor(Treatment, level = level_order3), y = abun)) + 
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
  labs(x = "Treatment", y = "Average VOC Abundance \n (peak area/g/h)")
dev.off()








# old code using biomass totals per treatment
# manually making a dataframe that takes abundance sums divided by plant biomass for each treatment
voc_transpose_w <- data.frame(Treatment = c("Ambient","Irrigated","Drought","Warmed","Warmed_Drought"),
                              weighted_abun = c(6.471801,6.714231,4.172067,6.209308,4.201861))
# total abundance from voc_transpose_sum, plant biomass from the 2022_Solidago_leaf_biomass_L2.R script
# ambient 394.3414/60.93225
# irrigated 441.2288/65.71546
# drought 240.5085/57.64733
# warmed 407.4361/65.61699
# warmed drought 211.2220/50.26868

# total abundance - including biomass
level_order <- c('Ambient', 'Irrigated', 'Drought', "Warmed", "Warmed_Drought")
ggplot(voc_transpose_w, aes(x = factor(Treatment, level = level_order), y = weighted_abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Total VOC abundance")


# manually making a dataframe that takes abundance sums divided by plant biomass for each treatment
voc_transpose_cmpd_w <- data.frame(Treatment = c("Ambient","Irrigated","Drought","Warmed","Warmed_Drought"),
                                   weighted_abun = c(0.06388608,0.008622839,0.01526052,0.09104409,0))
# total abundance from voc_transpose_sum, plant biomass from the 2022_Solidago_leaf_biomass_L2.R script
# ambient 3.8927223/60.93225
# irrigated 0.5666538/65.71546
# drought 0.8797281/57.64733
# warmed 5.9740389/65.61699
# warmed drought 0/50.26868

# total abundance - including biomass
ggplot(voc_transpose_cmpd_w, aes(x = factor(Treatment, level = level_order), y = weighted_abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed + \n Drought")) +
  labs(x = "Treatment", y = "Total VOC abundance")


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



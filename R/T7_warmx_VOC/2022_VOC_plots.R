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
library(paletteer)
library(ggthemes)

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
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.95, geom="polygon")+
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
  stat_ellipse(data=seg.data,aes(x=v.PCoA1,y=v.PCoA2,fill=group), alpha=.4,type='t',size =0.5, level=0.95, geom="polygon")+
  #geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4.7,color="black",shape=16) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color=group),alpha=0.7,size=3,shape=16) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2, fill=group),color="black",size=5,shape=21) + 
  scale_color_manual(labels=c("Ambient", "Warmed", "Drought", "Warmed + Drought"),
                     values=c('#2c7bb6',"khaki1","#fdae61","#d7191c"))+
  scale_fill_manual(labels=c("Ambient", "Warmed", "Drought", "Warmed + Drought"),
                    values=c('#2c7bb6',"khaki1","#fdae61","#d7191c"))+
  labs(x="PCoA 1",y="PCoA 2", color="Treatment", fill="Treatment") +
  annotate("text", x = -0.125, y=0.04, label = "A", size=5) +
  annotate("text", x = -0.07, y=0.06, label = "A", size=5) +
  annotate("text", x = 0.04, y=0.07, label = "B", size=5) +
  annotate("text", x = 0.1, y=0.09, label = "B", size=5) +
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




#### stacked bar plot - composition #####
# transforming from wide to long
voc_long <- voc_transpose_rm5 %>%
  pivot_longer(names_to = "compound", values_to = "avg_abun", cols = -c(Sample_ID,Unique_ID,Rep,Footprint,Subplot,Treatment,Notes))
# adding in full compounds names
voc_long$full_name <- NA
voc_long$full_name[voc_long$compound == "Ethanone..1..4.ethylphenyl.."] <- "Ethanone, 1-(4-ethylphenyl)-"
voc_long$full_name[voc_long$compound == "Salicylic.acid..tert..butyl.ester"] <- "Salicylic acid, tert.-butyl ester"
voc_long$full_name[voc_long$compound == "Butanenitrile..2.hydroxy.3.methyl."] <- "Butanenitrile, 2-hydroxy-3-methyl-"
voc_long$full_name[voc_long$compound == "X1.3.Bis.cyclopentyl..1.cyclopentanone"] <- "1,3-Bis(cyclopentyl)-1-cyclopentanone"
voc_long$full_name[voc_long$compound == "Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester"] <- "Propanoic acid, 2-methyl-, 3-hydroxy-2,2,4-trimethylpentyl ester"
voc_long$full_name[voc_long$compound == "X1.7.Nonadiene..4.8.dimethyl."] <- "1,7-Nonadiene, 4,8-dimethyl-"
voc_long$full_name[voc_long$compound == "X5.Hepten.2.one..6.methyl."] <- "5-Hepten-2-one, 6-methyl-"
voc_long$full_name[voc_long$compound == "dl.Menthol"] <- "dl-Menthol"
voc_long$full_name[voc_long$compound == "Pentane..2.bromo."] <- "Pentane, 2-bromo-"
voc_long$full_name[voc_long$compound == "X3.Heptanone..2.methyl."] <- "3-Heptanone, 2-methyl-"
voc_long$full_name[voc_long$compound == "Benzoic.acid..2.ethylhexyl.ester"] <- "Benzoic acid, 2-ethylhexyl ester"
voc_long$full_name[voc_long$compound == "X3.Butenoic.acid..ethyl.ester"] <- "3-Butenoic acid, ethyl ester"
voc_long$full_name[voc_long$compound == "Acetic.acid..1.1.dimethylethyl.ester"] <- "Acetic acid, 1,1-dimethylethyl ester"
voc_long$full_name[voc_long$compound == "Decane..2.4.dimethyl."] <- "Decane, 2,4-dimethyl-"
voc_long$full_name[voc_long$compound == "endo.Borneol"] <- "endo-Borneol"
voc_long$full_name[voc_long$compound == "X.Z.Z...alpha..Farnesene"] <- "(Z,Z)-alpha-Farnesene"
voc_long$full_name[voc_long$compound == "p.Cymene"] <- "p-Cymene"
voc_long$full_name[voc_long$compound == "X.....beta..Bourbonene"] <- "(-)-beta-Bourbonene"
voc_long$full_name[voc_long$compound == "X4.tert.Butylcyclohexyl.acetate"] <- "4-tert-Butylcyclohexyl acetate"
voc_long$full_name[voc_long$compound == "X6.10.Dimethyl.3..1.methylethylidene..1.cyclodecene"] <- "6,10-Dimethyl-3-(1-methylethylidene)-1-cyclodecene"
voc_long$full_name[voc_long$compound == "X2.Ethylhexyl.salicylate"] <- "2-Ethylhexyl salicylate"
voc_long$full_name[voc_long$compound == "Diisopropyl.adipate"] <- "Diisopropyl adipate"
voc_long$full_name[voc_long$compound == "X2.Cyclohexen.1.one"] <- "2-Cyclohexen-1-one"
voc_long$full_name[voc_long$compound == "o.Xylene"] <- "o-Xylene"
voc_long$full_name[voc_long$compound == "Styrene"] <- "Styrene"
voc_long$full_name[voc_long$compound == ".alpha..Bourbonene"] <- "alpha-Bourbonene"
voc_long$full_name[voc_long$compound == "X2.Hexene..2.5.dimethyl."] <- "2-Hexene, 2,5-dimethyl-"
voc_long$full_name[voc_long$compound == "X3.Hexen.1.ol"] <- "3-Hexen-1-ol"
voc_long$full_name[voc_long$compound == "Butane..1.ethoxy."] <- "Butane, 1-ethoxy-"
# adding in climate association - from the multiplatt comparisons in the VOC analyses script
voc_long$climate_group <- NA
voc_long$climate_group[voc_long$compound == "Ethanone..1..4.ethylphenyl.."] <- "Ambient"
voc_long$climate_group[voc_long$compound == "Salicylic.acid..tert..butyl.ester"] <- "Ambient"
voc_long$climate_group[voc_long$compound == "Butanenitrile..2.hydroxy.3.methyl."] <- "Ambient"
voc_long$climate_group[voc_long$compound == "X1.3.Bis.cyclopentyl..1.cyclopentanone"] <- "Warmed"
voc_long$climate_group[voc_long$compound == "Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester"] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "X1.7.Nonadiene..4.8.dimethyl."] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "X5.Hepten.2.one..6.methyl."] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "dl.Menthol"] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "Pentane..2.bromo."] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "X3.Heptanone..2.methyl."] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "Benzoic.acid..2.ethylhexyl.ester"] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "X3.Butenoic.acid..ethyl.ester"] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "Acetic.acid..1.1.dimethylethyl.ester"] <- "Warmed Drought"
voc_long$climate_group[voc_long$compound == "Decane..2.4.dimethyl."] <- "Ambient & Drought"
voc_long$climate_group[voc_long$compound == "endo.Borneol"] <- "Ambient & Warmed"
voc_long$climate_group[voc_long$compound == "X.Z.Z...alpha..Farnesene"] <- "Ambient & Warmed"
voc_long$climate_group[voc_long$compound == "p.Cymene"] <- "Ambient & Warmed"
voc_long$climate_group[voc_long$compound == "X.....beta..Bourbonene"] <- "Ambient & Warmed"
voc_long$climate_group[voc_long$compound == "X4.tert.Butylcyclohexyl.acetate"] <- "Drought & Warmed Drought"
voc_long$climate_group[voc_long$compound == "X6.10.Dimethyl.3..1.methylethylidene..1.cyclodecene"] <- "Drought & Warmed Drought"
voc_long$climate_group[voc_long$compound == "X2.Ethylhexyl.salicylate"] <- "Drought & Warmed Drought"
voc_long$climate_group[voc_long$compound == "Diisopropyl.adipate"] <- "Drought & Warmed Drought"
voc_long$climate_group[voc_long$compound == "X2.Cyclohexen.1.one"] <- "Drought & Warmed Drought"
voc_long$climate_group[voc_long$compound == "o.Xylene"] <- "Warmed & Warmed Drought"
voc_long$climate_group[voc_long$compound == "Styrene"] <- "Warmed & Warmed Drought"
voc_long$climate_group[voc_long$compound == ".alpha..Bourbonene"] <- "Ambient & Drought & Warmed"
voc_long$climate_group[voc_long$compound == "X2.Hexene..2.5.dimethyl."] <- "Drought & Warmed & Warmed Drought"
voc_long$climate_group[voc_long$compound == "X3.Hexen.1.ol"] <- "Drought & Warmed & Warmed Drought"
voc_long$climate_group[voc_long$compound == "Butane..1.ethoxy."] <- "Drought & Warmed & Warmed Drought"

# selecting compounds that popped up in indicator species analysis
# focusing on these compounds because including all of them would be too confusing on a fig
# first, all compounds associated w/ 1 treatment (but still allowing grouping, just not included here)
voc_long_cmpd <- voc_long %>%
  filter(compound == "Ethanone..1..4.ethylphenyl.." |
           compound == "Salicylic.acid..tert..butyl.ester" |
           compound == "Butanenitrile..2.hydroxy.3.methyl." |
           compound == "Butane..2.cyclopropyl." |
           compound == "X1.3.Bis.cyclopentyl..1.cyclopentanone" |
           compound == "X2.3.4.Trimethyl.1.pentanol" |
           compound == "Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester" |
           compound == "X1.7.Nonadiene..4.8.dimethyl." |
           compound == "X5.Hepten.2.one..6.methyl." |
           compound == "dl.Menthol" |
           compound == "Pentane..2.bromo." |
           compound == "X3.Heptanone..2.methyl." |
           compound == "Benzoic.acid..2.ethylhexyl.ester" |
           compound == "X3.Butenoic.acid..ethyl.ester " |
           compound == "Acetic.acid..1.1.dimethylethyl.ester") %>%
  group_by(Treatment, compound) %>%
  summarize(avg_abundance = mean(avg_abun*10000))
# second, all compounds associated w/ 1 treatment (not allowing grouping)
voc_long_cmpd2 <- voc_long %>%
  filter(compound == "Ethanone..1..4.ethylphenyl.." |
           compound == "Salicylic.acid..tert..butyl.ester" |
           compound == "Butanenitrile..2.hydroxy.3.methyl." |
           compound == "Bicyclo.3.2.0.hepta.2.6.diene" |
           compound == "X1.3.Bis.cyclopentyl..1.cyclopentanone" |
           compound == "X2.3.4.Trimethyl.1.pentanol" |
           compound == "Acetyl.valeryl" |
           compound == "Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester" |
           compound == "X4.tert.Butylcyclohexyl.acetate" |
           compound == "X1.7.Nonadiene..4.8.dimethyl." |
           compound == "Diisopropyl.adipate" |
           compound == "X5.Hepten.2.one..6.methyl." |
           compound == "X1.Hexanol..2.ethyl." |
           compound == "o.Xylene" |
           compound == "dl.Menthol" |
           compound == "Pentane..2.bromo." |
           compound == "Linalyl.acetate" |
           compound == "X3.Heptanone..2.methyl." |
           compound == "Benzoic.acid..2.ethylhexyl.ester" |
           compound == "Octane..2.methyl." |
           compound == "Acetic.acid..1.1.dimethylethyl.ester") %>%
  group_by(Treatment, compound) %>%
  summarize(avg_abundance = mean(avg_abun*10000))
# third, all compounds associated with 1, 2, or 3 treatments (allowing grouping) * using this one *
voc_long_cmpd3 <- voc_long %>%
  drop_na(full_name) %>%
  group_by(Treatment, full_name, climate_group) %>%
  summarize(avg_abundance = mean(avg_abun*10000))
# custom axis order
voc_long_cmpd2$Treatment <- factor(voc_long_cmpd2$Treatment, levels=c("Ambient_Control","Warmed",
                                                                      "Drought","Warmed_Drought"))
voc_long_cmpd3$Treatment <- factor(voc_long_cmpd3$Treatment, levels=c("Ambient_Control","Warmed",
                                                                      "Drought","Warmed_Drought"))

# stacked bar plot
# https://jkzorz.github.io/2019/06/05/stacked-bar-plots.html
ggplot(voc_long_cmpd2, aes(x = Treatment, fill = compound, y = avg_abundance)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        legend.position="none",
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Average VOC Abundance \n (peak area*10,000/g/hr)", fill = "Compound")

# bubble plot
# https://jkzorz.github.io/2019/06/05/Bubble-plots.html
colours = c( "#fe0000",  "#800001", "#fe6a00", "#803400","#ffd800","#806b00", "#00fe21", "#007f0e", "#0094fe",
             "#00497e","#0026ff","#001280","#b100fe","#590080")
colours2 = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC",
             "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
             "#771122", "#AA4455", "#DD7788")
colours3 = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 30)
# ordering compounds from highest abundance to lowest
voc_long_cmpd2$compound = reorder(voc_long_cmpd2$compound, voc_long_cmpd2$avg_abundance, FUN = mean)
voc_long_cmpd3$full_name = reorder(voc_long_cmpd3$full_name, voc_long_cmpd3$avg_abundance, FUN = mean)
# plot
png("bubble_vocs.png", units="in", width=13, height=7, res=300)
ggplot(voc_long_cmpd3, aes(x = Treatment, y = full_name)) + 
  geom_point(aes(size = avg_abundance, fill = full_name), alpha = 0.75, shape = 21) + 
  #scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Average VOC Abundance \n (peak area*10,000/g/hr)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey92"), 
        panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey92"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +
  scale_fill_manual(values = colours3, guide="none") +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Warmed" = "Warmed",
                            "Drought" = "Drought",
                            "Warmed_Drought" = "Warmed\nDrought"))
dev.off()
  #scale_y_discrete(labels=c("Acetic.acid..1.1.dimethylethyl.ester" = "Acetic acid, 1,1-dimethylethyl ester",
  #                          "Acetyl.valeryl" = "Acetyl valeryl",
  #                          "Benzoic.acid..2.ethylhexyl.ester" = "Benzoic acid, 2-ethylhexyl ester",
  #                          "Bicyclo.3.2.0.hepta.2.6.diene" = "Bicyclo[3.2.0]hepta-2,6-diene",
  #                          "Butanenitrile..2.hydroxy.3.methyl." = "Butanenitrile, 2-hydroxy-3-methyl-",
  #                          "Diisopropyl.adipate" = "Diisopropyl adipate",
  #                          "dl.Menthol" = "dl-Menthol",
  #                          "Ethanone..1..4.ethylphenyl.." = "Ethanone, 1-(4-ethylphenyl)-",
  #                          "Linalyl.acetate" = "Linalyl acetate",
  #                          "o.Xylene" = "o-Xylene",
  #                          "Octane..2.methyl." = "Octane, 2-methyl-",
  #                          "Pentane..2.bromo." = "Pentane, 2-bromo-",
  #                          "Propanoic.acid..2.methyl...3.hydroxy.2.2.4.trimethylpentyl.ester" = "Propanoic acid, 2-methyl-, 3-hydroxy-2,2,4-trimethylpentyl ester",
  #                          "Salicylic.acid..tert..butyl.ester" = "Salicylic acid, tert.-butyl ester",
  #                          "X1.3.Bis.cyclopentyl..1.cyclopentanone" = "1,3-Bis(cyclopentyl)-1-cyclopentanone",
  #                          "X1.7.Nonadiene..4.8.dimethyl." = "1,7-Nonadiene, 4,8-dimethyl-",
  #                          "X1.Hexanol..2.ethyl." = "1-Hexanol, 2-ethyl-",
  #                          "X2.3.4.Trimethyl.1.pentanol" = "2,3,4-Trimethyl-1-pentanol",
  #                          "X3.Heptanone..2.methyl." = "3-Heptanone, 2-methyl-",
  #                          "X4.tert.Butylcyclohexyl.acetate" = "4-tert-Butylcyclohexyl acetate",
  #                          "X5.Hepten.2.one..6.methyl." = "5-Hepten-2-one, 6-methyl-"))




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
            se = std.error(rowsums),
            count = n())
voc_transpose_sum_CI <- voc_transpose_rm5 %>%
  group_by(Treatment) %>%
  summarize(avg_abun = mean(rowsums),
            CI = mean(rowsums)-(qnorm(0.975)*sd(rowsums)/sqrt(length(rowsums))),
            CI_total = avg_abun-CI,
            count = n())
voc_transpose_sum3 <- voc_transpose_rm5 %>%
  group_by(Treatment) %>%
  summarize(avg_abun = mean(rowsums),
            se = std.error(rowsums))
# plot - abundance w/o considering biomass
level_order2 <- c('Ambient_Control', 'Irrigated_Control', 'Warmed', "Drought", "Warmed_Drought")
png("climate_ab_no_ir.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum2, aes(x = factor(Treatment, level = level_order2), y = avg_abun)) + 
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
                            "Warmed_Drought" = "Warmed\nDrought")) +
  labs(x = "Treatment", y = "Average VOC Abundance \n (peak area/g/hr)")
dev.off()

# dot plot
png("climate_ab_95.png", units="in", width=6, height=5, res=300)
ggplot(voc_transpose_sum_CI,aes(x = factor(Treatment, level = level_order2), y = avg_abun)) +
  geom_pointrange(aes(ymin=avg_abun-CI_total, ymax=avg_abun+CI_total), ,pch=21,size=1,position=position_dodge(0.3),fill="lightsteelblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  ylim(0.25,1.25) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
                            "Irrigated_Control" = "Irrigated",
                            "Warmed_Drought" = "Warmed\nDrought")) +
  labs(x = "Treatment", y = "Average VOC Abundance \n (peak area/g/hr)")
dev.off()



### abundances for specific compounds, found in the analyses script ###
level_order2 <- c('Ambient_Control', 'Irrigated_Control', 'Drought', "Warmed", "Warmed_Drought")
level_order3 <- c('Ambient_Control', 'Warmed','Drought', "Warmed_Drought")

# o xylene
voc_transpose_xylene <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, Rep, o.Xylene)
# taking average per treatment
voc_transpose_xylene_a <- voc_transpose_xylene %>%
  group_by(Treatment) %>%
  summarize(abun = mean(o.Xylene),
            se = std.error(o.Xylene))

# X.Z.Z...alpha..Farnesene
voc_transpose_farn <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, Rep, X.Z.Z...alpha..Farnesene)
# taking average per treatment
voc_transpose_farn_a <- voc_transpose_farn %>%
  group_by(Treatment) %>%
  summarize(abun = mean(X.Z.Z...alpha..Farnesene),
            se = std.error(X.Z.Z...alpha..Farnesene))

# beta Myrcene
voc_transpose_myrcene <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, Rep, .beta..Myrcene)
# taking average per treatment
voc_transpose_myrcene_a <- voc_transpose_myrcene %>%
  group_by(Treatment) %>%
  summarize(abun = mean(.beta..Myrcene),
            se = std.error(.beta..Myrcene))

# Butanenitrile..2.hydroxy.3.methyl.
voc_transpose_butane <- voc_transpose %>%
  dplyr::select(Sample_ID, Treatment, Rep, Butanenitrile..2.hydroxy.3.methyl.)
# taking average per treatment
voc_transpose_butane_a <- voc_transpose_butane %>%
  group_by(Treatment) %>%
  summarize(abun = mean(Butanenitrile..2.hydroxy.3.methyl.),
            se = std.error(Butanenitrile..2.hydroxy.3.methyl.))

# plot
png("climate_hexen.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose_xylene_a, aes(x = factor(Treatment, level = level_order2), y = abun)) + 
  geom_bar(position = "identity", stat = "identity", color = 'black', fill = "lightsteelblue3") +
  geom_errorbar(aes(ymin = abun - se, ymax = abun + se), width = 0.2,
                position = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15)) +
  scale_x_discrete(labels=c("Ambient_Control" = "Ambient",
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



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
library(plotrix)
library(plotly)
library(broom)
library(Rtsne)

# Set working directory
dir<-Sys.getenv("DATA_DIR")

# Read in data
voc_transpose <- read.csv(file.path(dir, "T7_warmx_VOC/L1/T7_VOC_2021drought_L1.csv"))



#### NMDS ####
# make community matrix - extract columns with abundance information
ab = voc_transpose[,2:359]

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
data.scores$Compound = voc_transpose$Compound
data.scores$Treatment = voc_transpose$Treatment
data.scores$Group_treat = voc_transpose$Group_treat
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

# 3D plot of data - to do this one, change K in nmds to 3
plot_ly(data.scores, x=~NMDS1, y=~NMDS2, z=~NMDS3, type="scatter3d", mode="markers", color=~Treatment)

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
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Treatment),size=3) + # add the point markers
  coord_equal() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())

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



#### PCoA ####
ab = voc_transpose[,2:359]
ab.dist<-vegdist(ab, method='bray')
dispersion<-betadisper(ab.dist, group=voc_transpose$Treatment)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion$centroids),data.frame(dispersion$centroids))
vectors<-data.frame(group=dispersion$group,data.frame(dispersion$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
png("drought_comp.png", units="in", width=6, height=5, res=300)
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



##### tSNE #####
# https://datavizpyr.com/how-to-make-tsne-plot-in-r/
set.seed(143)
# making rep a character column so it isn't included in numerical columns for tSNE
voc_transpose$Rep <- as.character(voc_transpose$Rep)
# making ID column w/ row identifiers
voc_transpose2 <- voc_transpose %>% 
  mutate(ID=row_number()) 
# selecting all non-numeric columns as the meta data
voc_meta <- voc_transpose2 %>%
  select(ID,Treatment,Rep,Group_treat)
# run tSNE
# calculating perplexity: 3 x perplexity < nrows - 1
tSNE_fit <- voc_transpose2 %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  Rtsne(perplexity=15)
# making dataframe
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())
# joining w meta data
tSNE_df <- tSNE_df %>%
  inner_join(voc_meta, by="ID")
# tSNE plots
ggplot(tSNE_df, aes(x = tSNE1, y = tSNE2, color = Treatment)) +
  geom_point() +
  theme(legend.position="right") +
  theme_classic()
ggplot(tSNE_df, aes(x = tSNE1, y = tSNE2, color = Group_treat)) +
  geom_point() +
  theme(legend.position="right") +
  theme_classic()



#### ABUNDANCE ####
voc_transpose$rowsums <- rowSums(voc_transpose[2:359])

voc_transpose_sum <- voc_transpose %>%
  group_by(Treatment) %>%
  summarize(abun = sum(rowsums),
            se = std.error(rowsums))

level_order <- c('Ambient', 'Irrigated', 'Warmed', "WarmedDrought", "Drought")
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
  geom_boxplot(color = 'black', fill = "cornflowerblue") +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Abundance", fill = "Treatment") +
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

# plots for each compound - compounds picked based on analyses in analyses script
level_order <- c('Ambient', 'Irrigated', 'Warmed', "WarmedDrought", "Drought")
png("drought_caro.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose, aes(x = factor(Treatment, level = level_order), y = Caryophyllene)) + 
  geom_boxplot(color = 'black', fill = "cornflowerblue") +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Abundance: Caryophyllene", fill = "Treatment") +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13,face="bold")) +
  scale_x_discrete(limits = c("Ambient", "Irrigated", "Warmed", "WarmedDrought", "Drought"),
                   labels=c("Ambient" = "Ambient",
                            "Drought" = "Drought",
                            "Irrigated" = "Irrigated",
                            "Warmed" = "Warmed",
                            "WarmedDrought" = "Warmed + Drought"),
                   guide = guide_axis(n.dodge=2))
dev.off()

ggplot(voc_transpose, aes(x = factor(Treatment, level = level_order), y = Cyclopentanone..2.cyclopentylidene.)) + 
  geom_boxplot(color = 'black', fill = "lightskyblue1") +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Relative Abundance - Cyclopentanone 2-cyclopentylidene-", fill = "Treatment")

png("drought_eth.png", units="in", width=6, height=4, res=300)
ggplot(voc_transpose, aes(x = factor(Treatment, level = level_order), y = Ethanone..1..4.ethylphenyl..)) + 
  geom_boxplot(color = 'black', fill = "cornflowerblue") +
  #geom_jitter(alpha = 0.5, color = "goldenrod") +
  theme_classic() +
  labs(x = "Treatment", y = "Relative Abundance: Ethanone, \n 1-(4-ethylphenyl)-", fill = "Treatment") +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13,face="bold")) +
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
  geom_point(size=3) +
  theme_classic()

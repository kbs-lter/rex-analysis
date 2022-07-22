# TITLE:          REX: T2 Wheat Measurements - height, greenness, biomass
# AUTHORS:        Moriah Young, Lisa Leonard
# COLLABORATORS:  Gran Falvo (biomass data)
# DATA INPUT:     Data imported as csv files from shared REX Google drive 
# DATA OUTPUT:    
# PROJECT:        REX
# DATE:           July 2022

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(plotrix)

# Set working directory from .Renviron
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Read in data
#wheat_data <-read.csv("C:\\Users\\lisal\\Downloads\\T2_height_greenness_2022_L0 - Sheet1 (2).csv")
wheat_data <- read.csv(file.path(dir, "T2_height_greenness_2022_L0.csv")) # Moriah
#biomass_data<-read.csv("C:\\Users\\lisal\\Downloads\\T2_biomass_2022_L0 - Sheet1.csv")
biomass_data <- read.csv(file.path(dir, "T2_biomass_2022_L0.csv"))

# remove Year 3 drought from the data from the wheat and biomass dataframes
wheat_data <- filter(wheat_data, Subplot_Descriptions != "drought_corn_control")
biomass_data <- filter(biomass_data, Subplot_Descriptions != "drought_corn_control" & 
                               Subplot_Descriptions != "drought_corn_fungicide")

# make height and greenness into two different data frames
height <- select(wheat_data, -9)
greenness <- select(wheat_data, -8)

# remove 6/13/2022 data - just want to look at the 6/27/2022 height data
height1 <- filter(height, Date == "6/27/2022")

# remove 6/27/2022 data - just want to look at the 6/13/2022 height data
greenness <- filter(greenness, Date == "6/13/2022") 

# Adjusting biomass data to reflect 1m2 - the frame used was 115 cm by 75 cm, which is 0.8625m2 in area
biomass_data$biomass_meter2 <- biomass_data$anpp/0.8625

# remove outlier from height data (row 73 - very low height value)
# height1 <- height1[-73,]

# look at data
View(height1)
str(height1)
summary(height1$Height_cm)
summary(greenness$Greenness)

#taking height avg with %>% and summarizing all the heights used in the datasheet, 
height_avg <-height1 %>% 
        group_by(Subplot_Descriptions) %>%
        summarize(height_avg=mean(Height_cm,na.rm=TRUE),
                  se_height=std.error(Height_cm,na.rm=TRUE))

#now do this for avg_greenness by plugging in greenness. 
avg_greenness <-greenness %>%
        group_by(Subplot_Descriptions)%>%
        summarize(avg_greenness=mean(Greenness,na.rm=TRUE),
                  se_greenness=std.error(Greenness,na.rm=TRUE))

#lastly,take biomass avg and do standard error
biomass_avg<-biomass_data %>%
        group_by(Subplot_Descriptions) %>%
        summarize(biomass_avg=mean(biomass_meter2,na.rm=TRUE),
                  se_biomass=std.error(biomass_meter2,na.rm=TRUE))

#merge all three datasets together
wheat_green_height<-merge.data.frame(height_avg,avg_greenness)

#now merge biomass as well
wheat_green_height_biomass<-merge(wheat_green_height,biomass_avg)

#Steal Moriah's code from Github,thx
#Plot for height
ggplot(height_avg, aes(x = Subplot_Descriptions, y = height_avg, fill = Subplot_Descriptions)) +
  geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = Subplot_Descriptions)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = height_avg - se_height, ymax = height_avg + se_height), width = 0.2,
                position = "identity") +
  scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
       "drought_control" = "deepskyblue", "drought_corn_control" = "chartreuse1", "drought_fungicide" = "lightsteelblue2",
       "drought_legacy_control"= "seagreen2", "drought_legacy_fungicide"= "turquoise")) +
  labs(x = "Treatment", y = "Average height (cm)") +
  scale_x_discrete(labels=c("control" = "Control",
                            "control_fungicide" = "Fungicide",
                          "drought_control" = "Drought",
                           "drought_corn_control" = "Year 3 \n Drought",
                           "drought_fungicide" = "Drought \n Fungicide",
                            "drought_legacy_control"= "Drought \n Legacy",   
                             "drought_legacy_fungicide"="Drought \n Legacy Fungicide")) +
  theme_classic() +
  theme(legend.position = "none")

# box plot for height
png("T2_height_2022.png", units="in", width=6, height=5, res=300)
ggplot(height1, aes(x = Subplot_Descriptions, y = Height_cm, fill = Subplot_Descriptions)) +
        geom_boxplot(color = "black", outlier.shape = NA) +
        labs(x = "Treatment", y = "Height (cm)", fill = "Subplot_Descriptions") +
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
                                     "drought_control" = "deepskyblue", "drought_fungicide" = "deeppink",
                                     "drought_legacy_control"= "plum4", "drought_legacy_fungicide"= "blue2")) +
        theme(legend.text = element_text(size=16),
              legend.title = element_text(size=16)) +
        scale_x_discrete(labels=c("control" = "Control",
                                  "control_fungicide" = "Fungicide",
                                  "drought_control" = "Drought",
                                  "drought_fungicide" = "Drought \n Fungicide",
                                  "drought_legacy_control"= "Drought \n Legacy",   
                                  "drought_legacy_fungicide"="Drought \n Legacy \n Fungicide")) +
        theme_classic() +
        theme(legend.position="none")
dev.off()
  
# now do this for greenness 
ggplot(wheat_green_height_biomass, aes(x = Subplot_Descriptions, y = avg_greenness, fill = Subplot_Descriptions)) +
  geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = Subplot_Descriptions)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_greenness- se_greenness, ymax = avg_greenness + se_greenness), width = 0.2,
                position = "identity") +
  scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
                               "drought_control" = "deepskyblue", "drought_corn_control" = "lightcoral", "drought_fungicide" = "deeppink",
                               "drought_legacy_control"= "plum4", "drought_legacy_fungicide"= "blue2")) +
  labs(x = "Treatment", y = "Average Greenness") +
  scale_x_discrete(labels=c("control" = "Control",
                            "control_fungicide" = "Fungicide",
                            "drought_control" = "Drought",
                            "drought_corn_control" = "Year 3 \n Drought",
                            "drought_fungicide" = "Drought \n Fungicide",
                            "drought_legacy_control"= "Drought \n Legacy",   
                            "drought_legacy_fungicide"="Drought \n Legacy Fungicide")) +
  theme_classic()+
  theme(legend.position="none")
  #Fungicide=red colors, drought=blue, both purple idk wjnrw;vofr.erw

# box plot for greenness
ggplot(greenness, aes(x = Subplot_Descriptions, y = Greenness, fill = Subplot_Descriptions)) +
        geom_boxplot(color = "black", outlier.shape = NA) +
        labs(x = "Treatment", y = "Greenness", fill = "Subplot_Descriptions") +
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
                                     "drought_control" = "deepskyblue", "drought_fungicide" = "deeppink",
                                     "drought_legacy_control"= "plum4", "drought_legacy_fungicide"= "blue2")) +
        theme(legend.text = element_text(size=16),
              legend.title = element_text(size=16)) +
        scale_x_discrete(labels=c("control" = "Control",
                                  "control_fungicide" = "Fungicide",
                                  "drought_control" = "Drought",
                                  "drought_fungicide" = "Drought \n Fungicide",
                                  "drought_legacy_control"= "Drought \n Legacy",   
                                  "drought_legacy_fungicide"="Drought \n Legacy Fungicide")) +
        theme_classic() +
        theme(legend.position="none")

#lastly, do this for biomass
ggplot(wheat_green_height_biomass, aes(x = Subplot_Descriptions, y = biomass_avg, fill = Subplot_Descriptions)) +
        geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.6, aes(colour = "Subplot_Descriptions")) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = biomass_avg - se_biomass, ymax = biomass_avg + se_biomass), width = 0.2,
                      position = "identity") +
        scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
                                     "drought_control" = "deepskyblue", "drought_corn_control" = "lightcoral", "drought_fungicide" = "deeppink",
                                     "drought_legacy_control"= "plum4", "drought_legacy_fungicide"= "blue2")) +
        labs(x = "Treatment", y = "Average Biomass (g)") +
        scale_x_discrete(labels=c("control" = "Control",
                                  "control_fungicide" = "Fungicide",
                                  "drought_control" = "Drought",
                                  "drought_corn_control" = "Year 3 \n Drought",
                                  "drought_fungicide" = "Drought \n Fungicide",
                                  "drought_legacy_control"= "Drought \n Legacy",   
                                  "drought_legacy_fungicide"="Drought \n Legacy Fungicide")) +
        theme_classic()+
        theme(legend.position="none")
#Try to do some boxplots for the height, greenness, biomass
#geom_box?

# box plot for biomass
png("T2_biomass_2022.png", units="in", width=6, height=5, res=300)
ggplot(biomass_data, aes(x = Subplot_Descriptions, y = biomass_meter2, fill = Subplot_Descriptions)) +
        geom_boxplot(color = "black", outlier.shape = NA) +
        labs(x = "Treatment", y = "Average Biomass (g)", fill = "Subplot_Descriptions") +
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        scale_fill_manual(values = c("control" = "orchid1", "control_fungicide" = "cyan2", 
                                     "drought_control" = "deepskyblue", "drought_fungicide" = "deeppink",
                                     "drought_legacy_control"= "plum4", "drought_legacy_fungicide"= "blue2")) +
        theme(legend.text = element_text(size=16),
             legend.title = element_text(size=16)) +
        scale_x_discrete(labels=c("control" = "Control",
                                 "control_fungicide" = "Fungicide",
                                 "drought_control" = "Drought",
                                 "drought_fungicide" = "Drought \n Fungicide",
                                "drought_legacy_control"= "Drought \n Legacy",   
                                 "drought_legacy_fungicide"="Drought \n Legacy \n Fungicide")) +
        theme_classic() +
        theme(legend.position="none")
dev.off()


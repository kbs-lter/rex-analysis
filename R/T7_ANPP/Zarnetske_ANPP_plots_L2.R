# TITLE:          REX: Biomass plots
# AUTHORS:        Jordan Zapata, Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Moriah Young, Mark Hammond
# DATA INPUT:     Data imported as csv files from shared REX Google drive T7_ANPP L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           July 2022


### Preliminary steps & setting up data ###
# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse) # ggplot2 and dplyr are contained within tidyverse
library(plotrix) # I use this package for the standard error function

# read in data
#data = read.csv("T7.csv") # Jordan's filepath
dir <- Sys.getenv("DATA_DIR")
anpp <- read.csv(file.path(dir, "T7_ANPP/L1/T7_Zarnetske_ANPP_L1.csv")) # Kara's filepath
str(data)

# removing treatments we aren't interested in
anpp <- anpp %>%
  filter(!(Subplot_Descriptions=='drought_insecticide' |
             Subplot_Descriptions=='insecticide' |
             Subplot_Descriptions=='warmed_drought_insecticide' |
             Subplot_Descriptions=='warmed_insecticide'))

# remove species we aren't interested in (i.e., unknowns and unsorted plant material)
anpp <- anpp %>%
  filter(!(Species_Code=='UNASTAE' |
             Species_Code=='UNDIC' |
             Species_Code=='UNFABAE' |
             Species_Code=='UNSRT' |
             Species_Code=='CORSPP' | # removing cornus species here (for now) bc it is not designated as native/exotic
             origin == 'Both' | # removing species listed as both native and exotic
             growth_habit=='Vine')) # removing vine growth form

# Specifying the order that we want the treatments to be listed in the figures on the x-axis
level_order <- c('irrigated', 'ambient', 'drought','warmed','warmed_drought')


### Figures ###
# Kara's code - used Jordan's code as a template
# dot plot
ggplot(anpp, aes(x = Subplot_Descriptions, y = Dried_Plant_Biomass_g)) +
  geom_point() +
  theme_classic()

# bar plot
ggplot(anpp, aes(x = Subplot_Descriptions, y = Dried_Plant_Biomass_g)) +
  geom_bar(colour = "sienna1", fill = "sienna1", stat='identity') + 
  ggtitle("Total Biomass Per Treatment") +
  theme_classic()

# bar plot 2
# to get average biomass per treatment, first sum total biomass per plot
# then, take the average of that
anpp_bar <- anpp %>%
  group_by(Field_Loc_Code, Subplot_Descriptions) %>%
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE)) %>%
  group_by(Subplot_Descriptions) %>%
  summarize(average_biomass = mean(sum_biomass, na.rm = TRUE),
            se = std.error(sum_biomass, na.rm = TRUE))
png("ANPP_L2_barplot_overall.png", units="in", width=8, height=6, res=300)
ggplot(anpp_bar, aes(x = Subplot_Descriptions, y = average_biomass)) +
  geom_bar(position = "identity", stat = "identity", col = "black", fill="sienna1") +
  geom_errorbar(aes(ymin = average_biomass - se, ymax = average_biomass + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Average Biomass (g)") +
  scale_x_discrete(limits = level_order,
                   labels=c("ambient" = "Ambient",
                            "warmed" = "Warmed",
                            "irrigated" = "Irrigated",
                            "drought" = "Drought",
                            "warmed_drought" = "Warmed Drought")) +
  theme_classic()
dev.off()

# box plot
anpp_box <- anpp %>%
  group_by(Field_Loc_Code, Subplot_Descriptions) %>%
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE))
ggplot(anpp_box, aes(x = Subplot_Descriptions, y = sum_biomass)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7, fill="sienna1") +
  geom_jitter(aes(alpha=0.6), shape=16, size=2, color="sienna") +
  labs(x = "Treatment", y = "Average Biomass (g)") +
  scale_x_discrete(limits = level_order,
                   labels=c("ambient" = "Ambient",
                            "warmed" = "Warmed",
                            "irrigated" = "Irrigated",
                            "drought" = "Drought",
                            "warmed_drought" = "Warmed Drought")) +
  theme_classic()

# biomass between forb/graminoid
anpp_growth_bar <- anpp %>%
  group_by(Field_Loc_Code, Subplot_Descriptions, growth_habit) %>%
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE)) %>%
  group_by(Subplot_Descriptions, growth_habit) %>%
  summarize(average_biomass = mean(sum_biomass, na.rm = TRUE),
            se = std.error(sum_biomass, na.rm = TRUE))
png("ANPP_L2_barplot_growthform.png", units="in", width=8, height=6, res=300)
ggplot(anpp_growth_bar, aes(x = Subplot_Descriptions, y = average_biomass, fill = growth_habit)) +
  geom_bar(position = "dodge", stat = "identity", col = "black") +
  geom_errorbar(aes(ymin = average_biomass - se, ymax = average_biomass + se), width = 0.2,
                position = position_dodge(0.9)) +
  labs(x = "Treatment", y = "Average Biomass (g)", fill="Growth Habit") +
  scale_x_discrete(limits = level_order,
                   labels=c("ambient" = "Ambient",
                            "warmed" = "Warmed",
                            "irrigated" = "Irrigated",
                            "drought" = "Drought",
                            "warmed_drought" = "Warmed Drought")) +
  scale_fill_manual(values=c("sienna1", "sienna")) +
  theme_classic()
dev.off()

# specific species
anpp_sooca <- anpp %>%
  filter(Species_Code == "SOOCA") %>%
  group_by(Field_Loc_Code, Subplot_Descriptions) %>%
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE)) %>%
  group_by(Subplot_Descriptions) %>%
  summarize(average_biomass = mean(sum_biomass, na.rm = TRUE),
            se = std.error(sum_biomass, na.rm = TRUE))
png("ANPP_L2_barplot_sooca.png", units="in", width=8, height=6, res=300)
ggplot(anpp_sooca, aes(x = Subplot_Descriptions, y = average_biomass)) +
  geom_bar(position = "identity", stat = "identity", col = "black", fill="sienna1") +
  geom_errorbar(aes(ymin = average_biomass - se, ymax = average_biomass + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Average Biomass of S. canadensis (g)") +
  scale_x_discrete(limits = level_order,
                   labels=c("ambient" = "Ambient",
                            "warmed" = "Warmed",
                            "irrigated" = "Irrigated",
                            "drought" = "Drought",
                            "warmed_drought" = "Warmed Drought")) +
  theme_classic()
dev.off()

anpp_phlpr <- anpp %>%
  filter(Species_Code == "PHLPR") %>%
  group_by(Field_Loc_Code, Subplot_Descriptions) %>%
  summarize(sum_biomass = sum(Dried_Plant_Biomass_g, na.rm = TRUE)) %>%
  group_by(Subplot_Descriptions) %>%
  summarize(average_biomass = mean(sum_biomass, na.rm = TRUE),
            se = std.error(sum_biomass, na.rm = TRUE))
png("ANPP_L2_barplot_phlpr.png", units="in", width=8, height=6, res=300)
ggplot(anpp_phlpr, aes(x = Subplot_Descriptions, y = average_biomass)) +
  geom_bar(position = "identity", stat = "identity", col = "black", fill="sienna1") +
  geom_errorbar(aes(ymin = average_biomass - se, ymax = average_biomass + se), width = 0.2,
                position = "identity") +
  labs(x = "Treatment", y = "Average Biomass of P. pratense (g)") +
  scale_x_discrete(limits = level_order,
                   labels=c("ambient" = "Ambient",
                            "warmed" = "Warmed",
                            "irrigated" = "Irrigated",
                            "drought" = "Drought",
                            "warmed_drought" = "Warmed Drought")) +
  theme_classic()
dev.off()


# Jordan's code
ggplot(data,aes(x=subplot,y=biomass))+
  geom_point() +
  ggtitle("Total Biomass per Treatment") +    
  easy_center_title()

ggplot(data,  
       aes(x=subplot,y=biomass))+                                            
  geom_bar(colour = "sienna1", fill = "sienna1", stat='identity') + 
  ggtitle("Total Biomass Per Treatment") +   
  easy_center_title()

data$growth = as.character(data$growth_habit)

datcl=data%>%            
  filter(growth == "Forb" |        
           growth == "Graminoid")

ggplot(datcl,  # tells ggplot what is the dataset (dat_sum, filtered)
       aes(x=growth,y=biomass))+                                            
  geom_bar(colour = "sienna1", fill = "sienna1", stat='identity') + 
  ggtitle("Total Biomass Between Forbs and Graminoids") +    
  easy_center_title()

datg=data%>%            
  filter(Species_Code == "SOOCA" |        
           Species_Code == "ASTSA" |        
           Species_Code == "PHLPR"|        
           Species_Code == "POAPR")

ggplot(datg,  
       aes(x=scientific_name,y=biomass))+                                            
  geom_bar(colour = "sienna1", fill = "sienna1", stat='identity') + 
  ggtitle("Total Biomass Between Select Species") +    
  easy_center_title()

# TITLE:          REX: Abiotic plots for goldenrod paper
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Mark Hammond, Moriah Young, Emily Parker
# DATA INPUT:     Data imported as csv files from shared REX Google drive sensor data L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           May 2024



# Clear all existing data
rm(list=ls())

# Set working directory
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Load packages
library(tidyverse)
library(plotrix)
library(ggpubr)

# Read in data
hobo_data <- read.csv(file.path(dir, "sensors/OTC Footprints/L1/T7_warmx_HOBO_L1.csv"))
soil_data <- read.csv(file.path(dir, "sensors/OTC Footprints/L1/T7_warmx_soil_sensors_L1.csv"))
str(hobo_data)
str(soil_data)
# re-making the data column a date
hobo_data$Date_Time <- as.POSIXct(hobo_data$Date_Time, format = "%Y-%m-%d %H:%M:%S")
soil_data$Date_Time <- as.POSIXct(soil_data$sample_datetime, format = "%Y-%m-%d %H:%M:%S")



### subsetting data ###
# create new dataframe with only temp data from May - October
hobo_season <- hobo_data
hobo_season$month <- format(hobo_season$Date_Time,format="%m")
hobo_season$year <- format(hobo_season$Date_Time,format="%Y")
hobo_season$hour <- format(hobo_season$Date_Time, format="%H")
hobo_season$day <- format(hobo_season$Date_Time, format="%d")
hobo_season$year_month_day <- format(hobo_season$Date_Time, format="%Y-%m-%d")
hobo_season$date <- paste0(hobo_season$month,"",hobo_season$day)
hobo_season$date <- as.numeric(hobo_season$date)
hobo_season_sum <- hobo_season %>%
  filter(!(year == "2023")) %>%
  filter(month >= "07") %>%
  filter(month <= "10")
# create new dataframe with only soil data from May - October
soil_season <- soil_data
soil_season$month <- format(soil_season$Date_Time,format="%m")
soil_season$year <- format(soil_season$Date_Time,format="%Y")
soil_season$hour <- format(soil_season$Date_Time, format="%H")
soil_season$day <- format(soil_season$Date_Time, format="%d")
soil_season$monthday <- format(soil_season$Date_Time, format="%m%d")
soil_season$month_day <- format(soil_season$Date_Time, format="%m-%d")
soil_season_sum <- soil_season %>%
  filter(!(year == "2023" | year == "2024")) %>%
  filter(month >= "07") %>%
  filter(month <= "10")
soil_drought_check <- soil_season %>%
  filter(!(year == "2023" | year == "2024")) %>%
  filter(monthday >= "0620") %>%
  filter(monthday <= "0831") %>%
  filter(Subplot_Descriptions == "drought")



### temperature average from July - Oct ###
# take average temp per rep, per treatment
hobo_temp_rep_avg <- hobo_season_sum %>%
  group_by(Rep, Treatment) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE))
# averaging over each rep so n=6 for each treatment
hobo_temp_avg <- hobo_temp_rep_avg %>%
  group_by(Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))

# take average temp per rep, per treatment, per year
hobo_temp_rep_avg_year <- hobo_season_sum %>%
  group_by(Rep, Treatment, year) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE))
# averaging over each rep so n=6 for each treatment
hobo_temp_avg_year <- hobo_temp_rep_avg_year %>%
  group_by(Treatment, year) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))



### soil temp and moisture average from July - Oct ###
# take average per day, per rep, per treatment
soil_sampling_avg_rep <- soil_season_sum %>%
  group_by(Rep, Subplot_Descriptions) %>%
  summarize(average_temp = mean(temperature, na.rm = TRUE),
            average_moist = mean(vwc, na.rm=T))
# averaging over each rep, n varies between treatments
soil_sampling_avg <- soil_sampling_avg_rep %>%
  group_by(Subplot_Descriptions) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se_temp = std.error(average_temp, na.rm = TRUE),
            avg_moist = mean(average_moist, na.rm=T),
            se_moist = std.error(average_moist, na.rm=T))

# take average per rep, per treatment, per year
soil_sampling_avg_rep_year <- soil_season_sum %>%
  group_by(Rep, Subplot_Descriptions, year) %>%
  summarize(average_temp = mean(temperature, na.rm = TRUE),
            average_moist = mean(vwc, na.rm=T))
# averaging over each rep, n varies between treatments
soil_sampling_avg_year <- soil_sampling_avg_rep_year %>%
  group_by(Subplot_Descriptions, year) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se_temp = std.error(average_temp, na.rm = TRUE),
            avg_moist = mean(average_moist, na.rm=T),
            se_moist = std.error(average_moist, na.rm=T))



### daily soil moisture from july-aug in drought treatment ###
dr_check <- soil_drought_check %>%
  group_by(month_day,year,Rep) %>%
  summarize(average_moist = mean(vwc, na.rm = TRUE),
            se_moist = std.error(vwc, na.rm=T))



### figures ###
# plot - air temp, soil temp, and soil moisture averaged over both years
level_order1 <- c("Ambient", 'Warmed', 'Drought',"Warmed_Drought") 
level_order2 <- c("irrigated_control","ambient", 'warmed', 'drought',"warmed_drought") 
air_temp <- ggplot(hobo_temp_avg, aes(x = factor(Treatment, level = level_order1), y = avg_temp)) +
  geom_pointrange(aes(ymin=avg_temp-se, ymax=avg_temp+se),pch=21,size=1,fill="#AE1F00") +
  labs(y="Air temperature (°C)", x=NULL) +
  scale_x_discrete(labels=c("Ambient" = "A", "Drought" = "D",
                            "Warmed" = "W",
                            "Warmed_Drought" = "WD")) +
  theme_bw() +
  #annotate("text", x = 0.6, y=23, label = "A", size=6) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
soil_temp <- ggplot(soil_sampling_avg, aes(x = factor(Subplot_Descriptions, level=level_order2), y = avg_temp)) +
  geom_pointrange(aes(ymin=avg_temp-se_temp, ymax=avg_temp+se_temp), pch=21,size=1,fill="#AE1F00") +
  #geom_errorbar(aes(ymin=avg_temp-se_temp, ymax=avg_temp+se_temp),width=0.1,color="black",linetype="solid") +
  #geom_point(size = 2) +
  labs(y="Soil temperature (°C)", x=NULL) +
  scale_x_discrete(labels=c("ambient" = "A", "drought" = "D",
                            "irrigated_control" = "I", "warmed" = "W",
                            "warmed_drought" = "WD")) +
  theme_bw() +
 # annotate("text", x = 0.7, y=21.3, label = "B", size=6) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
soil_moist <- ggplot(soil_sampling_avg, aes(x = factor(Subplot_Descriptions, level=level_order2), y = avg_moist)) +
  geom_pointrange(aes(ymin=avg_moist-se_moist, ymax=avg_moist+se_moist),pch=21,size=1,fill="#AE1F00") +
  #geom_errorbar(aes(ymin=avg_moist-se_moist, ymax=avg_moist+se_moist),width=0.1,color="black",linetype="solid") +
  #geom_point(size = 2) +
  labs(y=bquote("Soil moisture " (m^3/m^3)), x=NULL) +
  scale_x_discrete(labels=c("ambient" = "A", "drought" = "D",
                            "irrigated_control" = "I", "warmed" = "W",
                            "warmed_drought" = "WD")) +
  theme_bw() +
  #annotate("text", x = 0.7, y=0.28, label = "C", size=6) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
# combine plots
abiotic_comb <- ggpubr::ggarrange(air_temp,soil_temp,soil_moist,
                                  ncol = 3,widths = c(0.9,0.9,1))
png("rex_abiotic.png", units="in", width=12, height=5, res=300)
annotate_figure(abiotic_comb,
                bottom = text_grob("Treatment", color = "black", size=17))
dev.off()



# plot - air temp, soil temp, and soil moisture separate for both years
level_order1 <- c("Ambient", 'Warmed', 'Drought',"Warmed_Drought") 
level_order2 <- c("irrigated_control","ambient", 'warmed', 'drought',"warmed_drought") 
air_temp_yearly <- ggplot(hobo_temp_avg_year, aes(x = factor(Treatment, level = level_order1), y = avg_temp, fill=year)) +
  geom_pointrange(aes(ymin=avg_temp-se, ymax=avg_temp+se),pch=21,size=1,position=position_jitter(w=0.1)) +
  labs(y="Air temperature (°C)", x=NULL) +
  scale_x_discrete(labels=c("Ambient" = "A", "Drought" = "D",
                            "Warmed" = "W",
                            "Warmed_Drought" = "WD")) +
  theme_bw() +
  #annotate("text", x = 0.6, y=23, label = "A", size=6) +
  scale_fill_manual(name="Year",
                    values = c("#FFB451", "#0b0055")) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
soil_temp_yearly <- ggplot(soil_sampling_avg_year, aes(x = factor(Subplot_Descriptions, level=level_order2), y = avg_temp, fill=year)) +
  geom_pointrange(aes(ymin=avg_temp-se_temp, ymax=avg_temp+se_temp), pch=21,size=1,position=position_jitter(w=0.1)) +
  #geom_errorbar(aes(ymin=avg_temp-se_temp, ymax=avg_temp+se_temp),width=0.1,color="black",linetype="solid") +
  #geom_point(size = 2) +
  labs(y="Soil temperature (°C)", x=NULL) +
  scale_x_discrete(labels=c("ambient" = "A", "drought" = "D",
                            "irrigated_control" = "I", "warmed" = "W",
                            "warmed_drought" = "WD")) +
  theme_bw() +
  # annotate("text", x = 0.7, y=21.3, label = "B", size=6) +
  scale_fill_manual(name="Year",
                    values = c("#FFB451", "#0b0055")) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
soil_moist_yearly <- ggplot(soil_sampling_avg_year, aes(x = factor(Subplot_Descriptions, level=level_order2), y = avg_moist, fill=year)) +
  geom_pointrange(aes(ymin=avg_moist-se_moist, ymax=avg_moist+se_moist),pch=21,size=1,position=position_jitter(w=0.1)) +
  #geom_errorbar(aes(ymin=avg_moist-se_moist, ymax=avg_moist+se_moist),width=0.1,color="black",linetype="solid") +
  #geom_point(size = 2) +
  labs(y=bquote("Soil moisture " (m^3/m^3)), x=NULL) +
  scale_x_discrete(labels=c("ambient" = "A", "drought" = "D",
                            "irrigated_control" = "I", "warmed" = "W",
                            "warmed_drought" = "WD")) +
  theme_bw() +
  #annotate("text", x = 0.7, y=0.28, label = "C", size=6) +
  scale_fill_manual(name="Year",
                    values = c("#FFB451", "#0b0055")) +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
# combine plots
abiotic_comb_yearly <- ggpubr::ggarrange(air_temp_yearly,soil_temp_yearly,soil_moist_yearly,
                                  ncol = 3,widths = c(1,1,1), common.legend=T, legend="right")
png("rex_abiotic_yearly.png", units="in", width=12, height=5, res=300)
annotate_figure(abiotic_comb_yearly,
                bottom = text_grob("Treatment", color = "black", size=17))
dev.off()



# plot - daily avg. soil moisture from july-august 2021
dr_check$Rep <- as.character(dr_check$Rep)
dr_check_2021 <- dr_check %>%
  filter(year == "2021")
everyother <- function(x) x[seq_along(x) %% 4 == 0]
png("rex_drought_2021.png", units="in", width=12, height=5, res=300)
ggplot(dr_check_2021, aes(x = month_day, y = average_moist, group=Rep, color=Rep)) +
  geom_errorbar(aes(ymin=average_moist-se_moist, ymax=average_moist+se_moist),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(y="VWC", x="Date", title="2021") +
  scale_x_discrete(breaks=everyother) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        title=element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))
dev.off()



# plot - daily avg. soil moisture from july-august 2022
dr_check$Rep <- as.character(dr_check$Rep)
dr_check_2022 <- dr_check %>%
  filter(year == "2022")
everyother <- function(x) x[seq_along(x) %% 6 == 0]
png("rex_drought_2022.png", units="in", width=12, height=5, res=300)
ggplot(dr_check_2022, aes(x = month_day, y = average_moist, group=Rep, color=Rep)) +
  geom_errorbar(aes(ymin=average_moist-se_moist, ymax=average_moist+se_moist),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(y="VWC", x="Date", title="2022") +
  scale_x_discrete(breaks=everyother) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))
dev.off()






# TITLE:          REX: HOBO plots
# AUTHORS:        Kara Dobson
# COLLABORATORS:  Phoebe Zarnetske, Mark Hammond, Moriah Young, Emily Parker
# DATA INPUT:     Data imported as csv files from shared REX Google drive sensor data L1 folder
# DATA OUTPUT:    Plots of data
# PROJECT:        REX
# DATE:           Jan 2023


# Clear all existing data
rm(list=ls())

# Set working directory
dir <- Sys.getenv("DATA_DIR")
list.files(dir)

# Load packages
library(tidyverse)

# Read in data
hobo_data <- read.csv(file.path(dir, "sensors/Phoebe Footprints/L1/T7_warmx_HOBO_L1.csv"))
str(hobo_data)
# re-making the data column a date
hobo_data$Date_Time <- as.POSIXct(hobo_data$Date_Time, format = "%Y-%m-%d %H:%M:%S")

# create new dataframe with only data from april - august from 7 AM - 7 PM (growing season during the day)
# this is mainly what we'll work with for figures (could always adjust the parameters as needed)
hobo_season <- hobo_data
hobo_season$month <- format(hobo_season$Date_Time,format="%m")
hobo_season$year <- format(hobo_season$Date_Time,format="%Y")
hobo_season$hour <- format(hobo_season$Date_Time, format="%H")
hobo_season$day <- format(hobo_season$Date_Time, format="%d")

hobo_season <- hobo_season %>%
  filter(month > "03") %>%
  filter(month < "09") %>%
  filter(hour > "06") %>%
  filter(hour < "20")

# adding another column to make drought + warmed drought just ambient and warmed
# to see warming effect w/o drought
hobo_season$Treatment2 <- NA
hobo_season$Treatment2[hobo_season$Treatment == 'Drought'] = "Ambient"
hobo_season$Treatment2[hobo_season$Treatment == 'Warmed_Drought'] = "Warmed"
hobo_season$Treatment2[hobo_season$Treatment == 'Ambient'] = "Ambient"
hobo_season$Treatment2[hobo_season$Treatment == 'Warmed'] = "Warmed"

# taking monthly average
# create new dataframes for temperatures averaged for plotting
hobo_monthly <- hobo_season %>%
  group_by(month, Rep, Treatment) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE),
            se = std.error(Temperature_C, na.rm = TRUE))
# averaging over each rep, so n=6 for each data point (n=5 for D and WD)
hobo_monthly_avg <- hobo_monthly %>%  # by month and year, 1m temp only
  group_by(month, Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))

# plot
png("rex_hobo.png", units="in", width=6, height=4, res=300)
ggplot(hobo_monthly_avg, aes(x = month, y = avg_temp, group=Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin=avg_temp-se, ymax=avg_temp+se),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(name="Treatment",
                     values = c("#a6bddb", "#687689", "#fb6a4a", "#9D422E"),
                     labels=c("Ambient","Drought","Warmed", "Warmed Drought")) +
  labs(y="Temperature (°C)", x="Month") +
  theme_classic()
dev.off()


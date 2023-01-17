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


# plot
ggplot(hobo_data, aes(x = Date_Time, y = Temperature_C)) +
  geom_point(aes(col = Treatment), size = 1) +
  geom_smooth(method = "loess", aes(linetype = Treatment, col = Treatment), size = 1.2,
              fill = "grey90") +
  labs(x = "Month",y = "Temperature °C") +
  theme_minimal() + 
  theme(legend.position = "bottom")


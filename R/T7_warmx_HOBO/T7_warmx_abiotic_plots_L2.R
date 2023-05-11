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
library(plotrix)

# Read in data
hobo_data <- read.csv(file.path(dir, "sensors/OTC Footprints/L1/T7_warmx_HOBO_L1.csv"))
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
  filter(hour > "06") %>%
  filter(hour < "20")
hobo_season_sum <- hobo_season %>%
  filter(hour > "06") %>%
  filter(hour < "20") %>%
  filter(month > "05") %>%
  filter(month < "09")

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
hobo_monthly_avg <- hobo_monthly %>%
  group_by(month, Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))

# taking 2022 average from june 1 - july 15
hobo_season$date <- paste0(hobo_season$month,"",hobo_season$day)
hobo_season$date <- as.numeric(hobo_season$date)
# selecting 2022 and june 1 - july 15
hobo_2022 <- hobo_season %>%
  filter(year == 2022) %>%
  filter(date > "531") %>%
  filter(date < "716")
# limit to the reps I used
hobo_2022 <- hobo_2022 %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)
# take average per day, per rep, per treatment
hobo_2022_avg <- hobo_2022 %>%
  group_by(Rep, Treatment) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE))
# averaging over each rep so n=4 for each treatment
hobo_2022_avg <- hobo_2022_avg %>%
  group_by(Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))


# taking 2022 average for july 11 - july 15 (sampling period)
hobo_season$date <- paste0(hobo_season$month,"",hobo_season$day)
hobo_season$date <- as.numeric(hobo_season$date)
# selecting 2022 and june 1 - july 15
hobo_sampling <- hobo_season %>%
  filter(year == 2022) %>%
  filter(date > "710") %>%
  filter(date < "716")
# limit to the reps I used
hobo_sampling <- hobo_sampling %>%
  filter(Rep == 2 | Rep == 3 | Rep == 4 | Rep == 5)
# take average per day, per rep, per treatment
hobo_sampling_avg <- hobo_sampling %>%
  group_by(Rep, Treatment) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE))
# averaging over each rep so n=4 for each treatment
hobo_sampling_avg <- hobo_sampling_avg %>%
  group_by(Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))


# taking monthly average
# create new dataframes for lux averaged for plotting
# note: lux probably isn't super accurate here bc we use solar shields over the pendants
hobo_monthly_lux <- hobo_season %>%
  group_by(month, Rep, Treatment) %>%
  summarize(average_lux = mean(Light_lux, na.rm = TRUE),
            se = std.error(Light_lux, na.rm = TRUE))
# averaging over each rep, so n=6 for each data point (n=5 for D and WD)
hobo_monthly_avg_lux <- hobo_monthly_lux %>%
  group_by(month, Treatment) %>%
  summarize(avg_lux = mean(average_lux, na.rm = TRUE),
            se = std.error(average_lux, na.rm = TRUE))

# taking daily average for the drought period
# create new dataframes for temperatures averaged for plotting
hobo_daily_sum <- hobo_season_sum %>%
  group_by(month,day, Rep, Treatment) %>%
  summarize(average_temp = mean(Temperature_C, na.rm = TRUE),
            se = std.error(Temperature_C, na.rm = TRUE))
# averaging over each rep, so n=6 for each data point (n=5 for D and WD)
hobo_daily_avg <- hobo_daily_sum %>%
  group_by(month, day,Treatment) %>%
  summarize(avg_temp = mean(average_temp, na.rm = TRUE),
            se = std.error(average_temp, na.rm = TRUE))
# merge month+day into one column
hobo_daily_avg$date <- paste0(hobo_daily_avg$month,"_",hobo_daily_avg$day)
# selecting only June 20th - August 10th (drought window)
hobo_daily_avg <- hobo_daily_avg[97:284, ]


# subtracting ambient from warmed and drought from warmed drought to get difference - monthly
hobo_monthly_diff <- hobo_monthly_avg
hobo_monthly_diff <- hobo_monthly_diff %>%
  group_by(month) %>%
  mutate(diff = avg_temp[Treatment == "Warmed"] - avg_temp[Treatment == "Ambient"]) %>%
  mutate(diff2 = avg_temp[Treatment == "Warmed_Drought"] - avg_temp[Treatment == "Drought"])
# removing values that don't match each treatment
hobo_monthly_diff$diff2[hobo_monthly_diff$Treatment == 'Ambient'] = NA
hobo_monthly_diff$diff2[hobo_monthly_diff$Treatment == 'Warmed'] = NA
hobo_monthly_diff$diff[hobo_monthly_diff$Treatment == 'Warmed_Drought'] = NA
hobo_monthly_diff$diff[hobo_monthly_diff$Treatment == 'Drought'] = NA
# merging difference columns
hobo_monthly_diff <- hobo_monthly_diff %>%
  mutate(diff = coalesce(diff,diff2)) %>%
  select(month,Treatment,diff)
# adding new name column
hobo_monthly_diff$treat_diff <- NA
hobo_monthly_diff$treat_diff[hobo_monthly_diff$Treatment == 'Ambient'] = "Warmed - Ambient"
hobo_monthly_diff$treat_diff[hobo_monthly_diff$Treatment == 'Drought'] = "Warmed Drought - Drought"
hobo_monthly_diff <- hobo_monthly_diff %>%
  filter(!is.na(treat_diff)) %>%
  select(month, diff, treat_diff)


# subtracting ambient from warmed and drought from warmed drought to get difference - daily drought
hobo_daily_diff <- hobo_daily_avg
hobo_daily_diff <- hobo_daily_diff %>%
  group_by(date) %>%
  mutate(diff = avg_temp[Treatment == "Warmed"] - avg_temp[Treatment == "Ambient"]) %>%
  mutate(diff2 = avg_temp[Treatment == "Warmed_Drought"] - avg_temp[Treatment == "Drought"])
# removing values that don't match each treatment
hobo_daily_diff$diff2[hobo_daily_diff$Treatment == 'Ambient'] = NA
hobo_daily_diff$diff2[hobo_daily_diff$Treatment == 'Warmed'] = NA
hobo_daily_diff$diff[hobo_daily_diff$Treatment == 'Warmed_Drought'] = NA
hobo_daily_diff$diff[hobo_daily_diff$Treatment == 'Drought'] = NA
# merging difference columns
hobo_daily_diff <- hobo_daily_diff %>%
  mutate(diff = coalesce(diff,diff2)) %>%
  select(date,Treatment,diff)
# adding new name column
hobo_daily_diff$treat_diff <- NA
hobo_daily_diff$treat_diff[hobo_daily_diff$Treatment == 'Ambient'] = "Warmed - Ambient"
hobo_daily_diff$treat_diff[hobo_daily_diff$Treatment == 'Drought'] = "Warmed Drought - Drought"
hobo_daily_diff <- hobo_daily_diff %>%
  filter(!is.na(treat_diff)) %>%
  select(date, diff, treat_diff)


# plot - monthly
png("rex_hobo.png", units="in", width=7, height=4, res=300)
ggplot(hobo_monthly_avg_2022, aes(x = month, y = avg_temp, group=Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin=avg_temp-se, ymax=avg_temp+se),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(name="Treatment",
                     values = c("#a6bddb", "#687689", "#fb6a4a", "#9D422E"),
                     labels=c("Ambient","Drought","Warmed", "Warmed Drought")) +
  labs(y="Temperature (°C)", x="Month") +
  theme_bw() +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=15),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))
dev.off()

# plot - monthly lux
png("rex_hobo.png", units="in", width=6, height=4, res=300)
ggplot(hobo_monthly_avg_lux, aes(x = month, y = avg_lux, group=Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin=avg_lux-se, ymax=avg_lux+se),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(name="Treatment",
                     values = c("#a6bddb", "#687689", "#fb6a4a", "#9D422E"),
                     labels=c("Ambient","Drought","Warmed", "Warmed Drought")) +
  labs(y="Light availability (Lux)", x="Month") +
  theme_classic()
dev.off()

# plot - daily in the summer
png("rex_hobo.png", units="in", width=6, height=4, res=300)
ggplot(hobo_daily_avg, aes(x = date, y = avg_temp, group=Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin=avg_temp-se, ymax=avg_temp+se),width=0.1,color="black",linetype="solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(name="Treatment",
                     values = c("#a6bddb", "#687689", "#fb6a4a", "#9D422E"),
                     labels=c("Ambient","Drought","Warmed", "Warmed Drought")) +
  labs(y="Temperature (°C)", x="Month") +
  theme_classic()
dev.off()

# difference plot
png("rex_hobo_diff.png", units="in", width=6, height=4, res=300)
ggplot(hobo_monthly_diff, aes(x = month, y = diff, group=treat_diff, fill = treat_diff)) +
  geom_bar(position = "dodge", stat = "identity", color = 'black') +
  scale_fill_manual(name="Treatment",
                     values = c("#fb6a4a", "#9D422E")) +
  labs(y="Temperature difference (°C)", x="Month") +
  theme_classic()
dev.off()

# difference plot - daily summer
ggplot(hobo_daily_diff, aes(x = date, y = diff, group=treat_diff, fill = treat_diff)) +
  geom_bar(position = "dodge", stat = "identity", color = 'black') +
  scale_fill_manual(name="Treatment",
                    values = c("#fb6a4a", "#9D422E")) +
  labs(y="Temperature (°C)", x="Month") +
  theme_classic()


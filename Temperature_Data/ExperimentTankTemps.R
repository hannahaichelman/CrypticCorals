# This script plots the tank temperatures from the TVE experiment, both looking at temperature recorded from the Apex
# and from the Hobo loggers in each tank

# read in data and set working directory
library(dplyr)
library(gridExtra)
library(ggplot2)
library(scales)
library(psych)
library(tidyverse)
library(ggpubr)

setwd("/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/")

#### Apex Temperature Data ####
# using the master doc that doesn't have the first few days of data, the probes were mis-labeled and didn't match later measurements
apex <- read.csv("/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Apex/Apex_Master_Doc.csv")

head(apex)
tail(apex)
str(apex)

# we have data from 10/06/2016 to 12/11/2016

# rename columns to match treatments
# 1 = control 1, 2 = low var, 3 = moderate var, 4 = high var, Tmpx8 = control 2

apex <- apex %>%
  dplyr::rename(HighVarTemp = Temp4, HighVarpH = pH.4, ControlpH = pH.1, ControlTemp = Temp1, LowVarpH = pH.2,
                LowVarTemp = Temp2, ModVarpH = pH.3, ModVarTemp = Temp3, Control2pH = pHx8, Control2Temp = Tmpx8,
                datetime = Date)

head(apex)

apex$datetime<-strptime(apex$datetime, format="%m/%d/%y %H:%M")
apex$datetime_ct <- as.POSIXct(apex$datetime, format="%Y-%m-%dT%H:%M:%S")
# there are two internal implementations of date/time: POSIXct, which stores seconds since UNIX epoch
# (+some other data), and POSIXlt, which stores a list of day, month, year, hour, minute, second, etc.

# transform the data to long format so treatment is its own column
apex_long = apex %>%
  select(datetime_ct, HighVarTemp, ControlTemp, LowVarTemp, ModVarTemp, Control2Temp) %>%
  gather(treatment, Temp, HighVarTemp, ControlTemp, LowVarTemp, ModVarTemp, Control2Temp) %>%
  dplyr::filter(Temp > 10) # remove strange value near 0 in High Var treatment that seems to be when probe came out of water

apex_long$treatment <- factor(apex_long$treatment, levels = c("ControlTemp", "Control2Temp", "LowVarTemp", "ModVarTemp", "HighVarTemp"))
apex_long$Day<-format(apex_long$datetime_ct,"%D")
apex_long$Day<-as.POSIXct(apex_long$Day, format="%m/%d/%y")

head(apex_long)
str(apex_long)

#### Plot Apex Temps ####
cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

apex_plot <- apex_long %>%
  dplyr::filter(treatment != "Control2Temp") %>% # ignore this, and add it back in to scale_color_manual with color blue to include this treatment
  ggplot(aes(x = datetime_ct, y = Temp, color = treatment)) +
  geom_line() +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat_reds) +
  xlab('DateTime') +
  ylab(expression(paste("Temp (",degree,"C)")))+
  scale_x_datetime(breaks = as.POSIXct(c("2016-11-09 00:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-13 00:00:00")),
                   date_labels = "%m/%d/%y")+ # to plot x axis labels at important dates, same as hobo loggers
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  #ylim(26,34) +
  ggtitle("Tank Temperatures - Apex") +
  theme_classic()
apex_plot

ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Apex/TankTempApex_nocontrol2.pdf", apex_plot, width=10, height=5, units=c("in"), useDingbats=FALSE)

# plot a subset of dates that look nice as an example of what we were going for
apex_clipped <- apex_long %>%
  dplyr::filter(datetime_ct > "2016-10-30 00:00:00", datetime_ct < "2016-11-02 00:00:00") %>%
  dplyr::filter(treatment != "Control2Temp")

apex_clipped_plot <- ggplot(apex_clipped, aes(x = datetime_ct, y = Temp, color = treatment)) +
  geom_line() +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat_reds) +
  xlab('DateTime') +
  ylab(expression(paste("Temp (",degree,"C)")))+
  scale_y_continuous(breaks = seq(28,33,by = 0.5))+
  scale_x_datetime(breaks = as.POSIXct(c("2016-10-30 00:00:00","2016-10-31 00:00:00","2016-11-01 00:00:00","2016-11-02 00:00:00")),
                   labels = c("10/30/16 00:00","10/31/16 00:00","11/1/16 00:00","11/2/16 00:00")) +
  #geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  #ylim(26,34) +
  ggtitle("Tank Temperatures - Apex") +
  theme_classic()
apex_clipped_plot

ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Apex/TankTempApex_clipped_nocontrol2_oldcols.pdf", apex_clipped_plot, width=6, height=5, units=c("in"), useDingbats=FALSE)

# When originally plotting after combining individual files into the master doc, something was off -
# even after double checking and removing the first few days worth of data, the High Var and Control 2
# treatments switch after a few days in.
# Looked at the Apex_Master_Doc.csv, High Var Temp and Control 2 Temp switch on 10.9.16 between 12:30 and 12:40
# Manually fixed this issue and read back in the csv file

#### Apex Summary Stats ####
# Some summary stats for Apex temp data
describe.by(apex_long, group = "treatment")

#### Hobo Logger Temperature Data ####
Tank1.1 <- read.csv("Tank_Hobo_Loggers/TVE_1_1_clean.csv") # log every 5 mins, 9/5/16 17:00 to 12/16/16 14:15
str(Tank1.1)
head(Tank1.1)
# TVE_1_2 lots of missing data
# TVE_1_3 lots of missing data

Tank2.1 <- read.csv("Tank_Hobo_Loggers/TVE_2_1_clean.csv") # log every 5 mins, 9/5/16 17:00 to 12/16/16 14:30
str(Tank2.1)
head(Tank2.1)
# TVE_2_2 lots of missing data
# TVE_2_3 lots of missing data at the end of the experiment

Tank3.1 <- read.csv("Tank_Hobo_Loggers/TVE_3_1_clean.csv") # log every 5 mins, 9/5/16 17:00 to 11/15/16, but
# starts logging at odd intervals and eventually 1 minute starting 11/14/16 12:30, so trim here.
str(Tank3.1)
head(Tank3.1)
# TVE_3_1 lots of missing data starting 11/14
# TVE_3_3 lots of missing data, only through 9/13/2016
# No TVE_3_2

Tank4.2 <- read.csv("Tank_Hobo_Loggers/TVE_4_2_clean.csv") # log every 5 mins, 9/5/16 17:00 to 12/16/16 14:50
str(Tank4.2)
head(Tank4.2)
# TVE_4_1 lots of missing data, unusable
# no TVE_4_3 exists

Tank5.3 <- read.csv("Tank_Hobo_Loggers/TVE_5_3_clean.csv") # log every 5 mins, 9/5/16 17:00 to 12/16/16 14:50
str(Tank5.3)
head(Tank5.3)
# TVE_5_1 lots of missing data, unusable
# TVE_5_2 lots of missing data at the end of the experiment

# convert datetime format for all loggers and add Day column
Tank1.1$DateTime<-strptime(Tank1.1$DateTime, format="%m/%d/%y %H:%M")
Tank1.1$DateTime_ct <- as.POSIXct(Tank1.1$DateTime, format="%Y-%m-%dT%H:%M:%S")
Tank1.1$Day<-format(Tank1.1$DateTime,"%D")
Tank1.1$Day<-as.POSIXct(Tank1.1$Day, format="%m/%d/%y")
head(Tank1.1)

Tank2.1$DateTime<-strptime(Tank2.1$DateTime, format="%m/%d/%y %H:%M")
Tank2.1$DateTime_ct <- as.POSIXct(Tank2.1$DateTime, format="%Y-%m-%dT%H:%M:%S")
Tank2.1$Day<-format(Tank2.1$DateTime,"%D")
Tank2.1$Day<-as.POSIXct(Tank2.1$Day, format="%m/%d/%y")
head(Tank2.1)

Tank3.1$DateTime<-strptime(Tank3.1$DateTime, format="%m/%d/%y %H:%M")
Tank3.1$DateTime_ct <- as.POSIXct(Tank3.1$DateTime, format="%Y-%m-%dT%H:%M:%S")
Tank3.1$Day<-format(Tank3.1$DateTime,"%D")
Tank3.1$Day<-as.POSIXct(Tank3.1$Day, format="%m/%d/%y")
head(Tank3.1)

Tank4.2$DateTime<-strptime(Tank4.2$DateTime, format="%m/%d/%y %H:%M")
Tank4.2$DateTime_ct <- as.POSIXct(Tank4.2$DateTime, format="%Y-%m-%dT%H:%M:%S")
Tank4.2$Day<-format(Tank4.2$DateTime,"%D")
Tank4.2$Day<-as.POSIXct(Tank4.2$Day, format="%m/%d/%y")
head(Tank4.2)

Tank5.3$DateTime<-strptime(Tank5.3$DateTime, format="%m/%d/%y %H:%M")
Tank5.3$DateTime_ct <- as.POSIXct(Tank5.3$DateTime, format="%Y-%m-%dT%H:%M:%S")
Tank5.3$Day<-format(Tank5.3$DateTime,"%D")
Tank5.3$Day<-as.POSIXct(Tank5.3$Day, format="%m/%d/%y")
head(Tank5.3)

#### Calibrating Hobos to Apex ####
# Need to try and find a way to calibrate the Hobo data to the Apex, since the Apex probes were calibrated with
# NIST-certified glass thermometers but the Hobo loggers were not calibrated

# Apex: we have data from 10/06/2016 to 12/11/2016
# Hobo: data from 9/5/16 to 12/16/16, but for moderate variability treatment, data stops on 11/15/16

# Cut all of the hobo logger data to be only the same time frame as the Apex data to find an average difference in temp
# then calculate means by day for Hobo loggers and apex and calculate average difference

# Hobo data:
Tank1.1_dayavg <- Tank1.1 %>%
  dplyr::filter(DateTime >= "2016-10-06 00:00:00", DateTime <= "2016-11-15 23:50:00") %>%
  group_by(Day) %>%
  summarise_at(vars(Temp), list(mean_temp = mean), na.rm=TRUE)

Tank2.1_dayavg <- Tank2.1 %>%
  dplyr::filter(DateTime >= "2016-10-06 00:00:00", DateTime <= "2016-11-15 23:50:00") %>%
  group_by(Day) %>%
  summarise_at(vars(Temp), list(mean_temp = mean), na.rm=TRUE)

Tank3.1_dayavg <- Tank3.1 %>%
  dplyr::filter(DateTime >= "2016-10-06 00:00:00", DateTime <= "2016-11-15 23:50:00") %>%
  group_by(Day) %>%
  summarise_at(vars(Temp), list(mean_temp = mean), na.rm=TRUE)

Tank4.2_dayavg <- Tank4.2 %>%
  dplyr::filter(DateTime >= "2016-10-06 00:00:00", DateTime <= "2016-11-15 23:50:00") %>%
  group_by(Day) %>%
  summarise_at(vars(Temp), list(mean_temp = mean), na.rm=TRUE)

# Calculate means by day for Apex data
apex_dayavg <- apex_long %>%
  group_by(Day, treatment) %>%
  summarise_at(vars(Temp), list(mean_temp = mean), na.rm=TRUE)

# Subset apex by treatment, have to cut to 2016-11-15 because that is the limit
# of where we have data for both apex and hobo in the moderate variability treatment,
# but want to consider the same days for all treatments

# Control: 2016-10-06 to 2016-12-11
apex_dayavg_control <-  apex_dayavg %>%
  dplyr::filter(treatment == "ControlTemp") %>%
  dplyr::filter(Day <= "2016-11-15")

# Low Var: 2016-10-06 to 2016-12-11
apex_dayavg_lowvar <-  apex_dayavg %>%
  dplyr::filter(treatment == "LowVarTemp") %>%
  dplyr::filter(Day <= "2016-11-15")

# Mod Var: 2016-10-06 to 2016-11-15
apex_dayavg_modvar <-  apex_dayavg %>%
  dplyr::filter(treatment == "ModVarTemp") %>%
  dplyr::filter(Day <= "2016-11-15")

# High Var: 2016-10-06 to 2016-12-11
apex_dayavg_highvar <-  apex_dayavg %>%
  dplyr::filter(treatment == "HighVarTemp") %>%
  dplyr::filter(Day <= "2016-11-15")

# Compare average temp difference of Hobo to Apex by treatment

mean(apex_dayavg_control$mean_temp - Tank1.1_dayavg$mean_temp) # control
# [1] 0.6860648
mean(apex_dayavg_lowvar$mean_temp - Tank2.1_dayavg$mean_temp) # low variability
# [1] 0.724092
mean(apex_dayavg_modvar$mean_temp - Tank3.1_dayavg$mean_temp) # moderate variability
# [1] 0.7730094
mean(apex_dayavg_highvar$mean_temp - Tank4.2_dayavg$mean_temp) # high variability
# [1] 0.7088613

# Plot the difference over time by treatment
plot(apex_dayavg_control$mean_temp - Tank1.1_dayavg$mean_temp, main = "Control Treatment", ylab = "Apex - Hobo (°C)", ylim = c(0,1.1)) # control
plot(apex_dayavg_lowvar$mean_temp - Tank2.1_dayavg$mean_temp, main = "Low Variability Treatment", ylab = "Apex - Hobo (°C)", ylim = c(0,1.1)) # low variability
plot(apex_dayavg_modvar$mean_temp - Tank3.1_dayavg$mean_temp, main = "Moderate Variability Treatment", ylab = "Apex - Hobo (°C)", ylim = c(0,1.1)) # moderate variability
plot(apex_dayavg_highvar$mean_temp - Tank4.2_dayavg$mean_temp, main = "High Variability Treatment", ylab = "Apex - Hobo (°C)", ylim = c(-0.1,2.2)) # high variability

# There are some strange daily outliers, likely where probes came out of the water for Apex, but
# moving forward with correcting the Hobo logger data (since it is most complete) to the Apex data
# (since Apex was calibrated with a thermometer) based on treatment-specific mean difference

#### Correct Hobos to Apex and trim files ####

## Trim temp files to only include dates in the experiment & separately to only include variability &
## again to include heat stress and recovery

# 9/6/2016 was the last day of recovery, the first day of acclimation started 9/7/2016.
# 9/22/2016 was the first day of TVE treatments.
# 11/10/2016 was the last day of TVE treatments.
# 11/11/2016 - 11/14/2016 ramped up temps from 28 to 32C
# 11/14/2016 - 11/21/2016 held at 32C
# 11/22/2016 - 11/25/2016 ramped down from 32 to 28
# 11/26/2016 - 12/11/2016 recovery at 28 - calling this end of recovery
# The experiment ended on 12/13/2016, everyone frozen

# System 5 didn't get put on Apex until day 10 of acclimation (9/16/16)

Tank1.1_all <- Tank1.1 %>%
  mutate(Temp = Temp + 0.686064) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-12 00:00:00" & DateTime > "2016-09-07 00:00:00") %>%
  mutate(Treatment = "Control")

Tank1.1_var <- Tank1.1 %>%
  mutate(Temp = Temp + 0.686064) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-11-10 18:00:00" & DateTime > "2016-09-22 00:00:00") %>%
  mutate(Treatment = "Control")

Tank1.1_stress <- Tank1.1 %>%
  mutate(Temp = Temp + 0.686064) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-10 18:00:00") %>%
  mutate(Treatment = "Control")

Tank2.1_all <- Tank2.1 %>%
  mutate(Temp = Temp + 0.724092) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-12 00:00:00" & DateTime > "2016-09-07 00:00:00") %>%
  mutate(Treatment = "Low Var")

Tank2.1_var <- Tank2.1 %>%
  mutate(Temp = Temp + 0.724092) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-11-10 18:00:00" & DateTime > "2016-09-22 00:00:00") %>%
  mutate(Treatment = "Low Var")

Tank2.1_stress <- Tank2.1 %>%
  mutate(Temp = Temp + 0.724092) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-10 18:00:00") %>%
  mutate(Treatment = "Low Var")

# Tank 3 doesn't record until the end of the experiment
Tank3.1_all <- Tank3.1 %>%
  mutate(Temp = Temp + 0.7730094) %>% # treatment-specific correction
  dplyr::filter(DateTime > "2016-09-07 00:00:00") %>%
  mutate(Treatment = "Mod Var")

Tank3.1_var <- Tank3.1 %>%
  mutate(Temp = Temp + 0.7730094) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-11-10 18:00:00" & DateTime > "2016-09-22 00:00:00") %>%
  mutate(Treatment = "Mod Var")

Tank3.1_stress <- Tank3.1 %>%
  mutate(Temp = Temp + 0.7730094) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-10 18:00:00") %>%
  mutate(Treatment = "Mod Var")

Tank4.2_all <- Tank4.2 %>%
  mutate(Temp = Temp + 0.7088613) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-12 00:00:00" & DateTime > "2016-09-07 00:00:00") %>%
  mutate(Treatment = "High Var")

Tank4.2_var <- Tank4.2 %>%
  mutate(Temp = Temp + 0.7088613) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-11-10 18:00:00" & DateTime > "2016-09-22 00:00:00") %>%
  mutate(Treatment = "High Var")

Tank4.2_stress <- Tank4.2 %>%
  mutate(Temp = Temp + 0.7088613) %>% # treatment-specific correction
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-10 18:00:00") %>%
  mutate(Treatment = "High Var")

# Tank 5 = Control 2, not using this treatment so no treatment-specific correction for Tank5's
Tank5.3_all <- Tank5.3 %>%
  dplyr::filter(DateTime < "2016-12-12 00:00:00" & DateTime > "2016-09-07 00:00:00") %>%
  mutate(Treatment = "control 2")

Tank5.3_var <- Tank5.3 %>%
  dplyr::filter(DateTime < "2016-11-10 18:00:00" & DateTime > "2016-09-22 00:00:00") %>%
  mutate(Treatment = "control 2")

Tank5.3_stress <- Tank5.3 %>%
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-10 18:00:00") %>%
  mutate(Treatment = "control 2")


#### Statistically compare hobo logger data ####
# combine all variability data into one dataframe, exclude control 2 for now
all.hobo.data <- rbind(Tank1.1_var, Tank2.1_var, Tank3.1_var, Tank4.2_var)

library(car)

# Levene test used to test if samples have equal variance
leveneTest(Temp~Treatment, data = all.hobo.data)
# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)
# group     3  9061.9 < 2.2e-16 ***
#       57304


# now try https://stackoverflow.com/questions/43646987/multiple-comparison-post-hoc-test-for-levenes-test
# first calculate the median by treatment and add as a column
all.hobo.data <- all.hobo.data %>%
  group_by(Treatment) %>%
  mutate(temp.med = ifelse(Temp, median(Temp, na.rm=TRUE), ifelse(Temp==NA, NA)))

# then calculate the residual for each temperature
all.hobo.data$temp.med.res<-abs(all.hobo.data$Temp-all.hobo.data$temp.med)

# Then we run an ANOVA, and post-hoc if necessary:
levene.dat.aov<-aov(temp.med.res~Treatment, all.hobo.data)
summary(levene.dat.aov)
#                 Df Sum Sq Mean Sq F value Pr(>F)
# Treatment       3   9120  3040.1    9062 <2e-16 ***
# Residuals   57304  19224     0.3

TukeyHSD(levene.dat.aov)
# $treat
#                             diff        lwr        upr p adj
# high var-control       1.0918213  1.0742404  1.1094022     0
# low var-control        0.4477593  0.4301785  0.4653402     0
# moderate var-control   0.7264947  0.7089138  0.7440755     0
# low var-high var      -0.6440620 -0.6616429 -0.6264811     0
# moderate var-high var -0.3653267 -0.3829075 -0.3477458     0
# moderate var-low var   0.2787353  0.2611544  0.2963162     0

#### Plot Hobo Loggers Separately to Explore ####
# create separate plots for each treatment, and for the duration of variability only, stitch together with grid.arrange()
#cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

p1 = Tank1.1_all %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = 'darkgrey', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Control") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p1

p1_var = Tank1.1_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = 'darkgrey', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Control") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p1_var

p2 = Tank2.1_all %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#FF9966', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Low Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p2

p2_var = Tank2.1_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#FF9966', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Low Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p2_var

p3 = Tank3.1_all %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#CC3300', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Moderate Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p3

p3_var = Tank3.1_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#CC3300', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Moderate Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p3_var

p4 = Tank4.2_all %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#7f0000', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("High Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p4

p4_var = Tank4.2_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = '#7f0000', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("High Variability") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p4_var

p5 = Tank5.3 %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = 'blue', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Control 2") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p5

p5_var = Tank5.3_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp))+
  geom_line(color = 'blue', lwd=1)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  ggtitle("Control 2") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
p5_var

# ignore control 2, all data
hobo_plot <- grid.arrange(p1, p2, p3, p4, nrow = 2)
ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo.pdf", hobo_plot, width=12, height=6, units=c("in"), useDingbats=FALSE)

# ignore control 2, only variability data
hobo_plot_var <- grid.arrange(p1_var, p2_var, p3_var, p4_var, nrow = 2)
ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_var.pdf", hobo_plot_var, width=8, height=5, units=c("in"), useDingbats=FALSE)

# include control 2, all data
hobo_plot_ctrl2 <- grid.arrange(p1, p2, p3, p4, p5, nrow = 3)
ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_ctrl2.pdf", hobo_plot_ctrl2, width=10, height=6, units=c("in"), useDingbats=FALSE)

# include control 2, only variability data
hobo_plot_ctrl2_var <- grid.arrange(p1_var, p2_var, p3_var, p4_var, p5_var, nrow = 3)
ggsave(file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_ctrl2_var.pdf", hobo_plot_ctrl2_var, width=8, height=6, units=c("in"), useDingbats=FALSE)

#### Combine all hobo variability objects to plot ####
all_temp = rbind(Tank1.1_all, Tank2.1_all, Tank3.1_all, Tank4.2_all)
head(all_temp)
str(all_temp)
all_temp$Treatment = factor(all_temp$Treatment, levels = c("Control","Low Var","Mod Var","High Var"))

all_var = rbind(Tank1.1_var, Tank2.1_var, Tank3.1_var, Tank4.2_var)
head(all_var)
str(all_var)
all_var$Treatment = factor(all_var$Treatment, levels = c("Control","Low Var","Mod Var","High Var"))

all_stress = rbind(Tank1.1_stress, Tank2.1_stress, Tank3.1_stress, Tank4.2_stress)
head(all_stress)
str(all_stress)
all_stress$Treatment = factor(all_stress$Treatment, levels = c("Control","Low Var","Mod Var","High Var"))

cols_treat_reds <- c("Control" = "darkgrey", "Low Var" = "#FF9966", "Mod Var"="#CC3300", "High Var"="#7f0000")

# plot the variability period
str(all_var)
all_var$Treatment = factor(all_var$Treatment, levels = c("High Var", "Mod Var", "Low Var", "Control"))

all_var.plot = all_var %>%
  ggplot(aes(x = DateTime_ct, y = Temp, color = Treatment))+
  geom_line(aes(color = Treatment), lwd=1, alpha = 0.8)+
  scale_color_manual(values = cols_treat_reds)+
  ylim(26,34) +
  geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  #ggtitle("Control 2") +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-22 00:00:00","2016-11-10 18:00:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-22 00:00:00", "2016-11-10 18:00:00")))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
all_var.plot
ggsave(all_var.plot,file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_Lines_Combined.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

# plot the whole time course - figure 1C
all_temp$Treatment = factor(all_temp$Treatment, levels = c("High Var", "Mod Var", "Low Var", "Control"))

all_temp.plot = all_temp %>%
  ggplot(aes(x = DateTime_ct, y = Temp, color = Treatment))+
  annotate("rect", xmin = as.POSIXct("2016-09-07 00:00:00"), xmax = as.POSIXct("2016-09-22 00:00:00"), ymin = - Inf, ymax = Inf, fill = "gray", alpha = 0.3)+
  annotate("rect", xmin = as.POSIXct("2016-09-22 00:00:00"), xmax = as.POSIXct("2016-11-10 18:00:00"), ymin = - Inf, ymax = Inf, fill = "orange", alpha = 0.15)+
  annotate("rect", xmin = as.POSIXct("2016-11-10 18:00:00"), xmax = as.POSIXct("2016-11-25 00:00:00"), ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = as.POSIXct("2016-11-25 00:00:00"), xmax = as.POSIXct("2016-12-11 23:55:00"), ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_line(aes(color = Treatment), lwd=0.75)+
  scale_color_manual(values = cols_treat_reds)+
  #geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=.5) +
  xlab("Day") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-09-07 00:00:00", "2016-09-22 00:00:00","2016-11-10 18:00:00","2016-11-14 00:00:00","2016-11-21 00:00:00","2016-11-25 00:00:00", "2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-09-07 00:00:00", "2016-12-11 23:55:00")))+
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(26,33,1))+
  theme_bw()+
  #theme(axis.text.x = element_blank())
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
all_temp.plot
ggsave(all_temp.plot,file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_AllData_Combined.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

# plot a subset of days to zoom in and illustrate what the profiles looked like
all_temp.plot.subset = all_temp %>%
  subset(DateTime_ct > "2016-10-30" & DateTime_ct <= "2016-11-02") %>%
  ggplot(aes(x = DateTime_ct, y = Temp, color = Treatment))+
  annotate("rect", xmin = as.POSIXct(-Inf), xmax = as.POSIXct(Inf), ymin = - Inf, ymax = Inf, fill = "orange", alpha = 0.15)+
  geom_line(aes(color = Treatment), lwd=1, alpha = 1)+
  scale_color_manual(values = cols_treat_reds)+
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(26,33,1))+
  xlab("Day") +
  theme_bw() +
  theme(axis.text.x = element_blank())
#theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
all_temp.plot.subset
ggsave(all_temp.plot.subset,file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_Subset_Combined.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# plot all stress + recovery periods
all_stress.plot = all_stress %>%
  ggplot(aes(x = DateTime_ct, y = Temp, color = Treatment))+
  geom_line(aes(color = Treatment), lwd=1, alpha = 0.8)+
  scale_color_manual(values = cols_treat_reds)+
  #geom_hline(aes(yintercept = 28.5), colour="grey", linetype="dashed", lwd=1) +
  xlab("Date") +
  scale_x_datetime(breaks = as.POSIXct(c("2016-11-10 18:05:00","2016-11-14", "2016-11-21","2016-11-25", "2016-12-11 23:55:00")),
                   date_labels = "%m/%d/%y",
                   limits = as.POSIXct(c("2016-11-10 18:05:00","2016-12-11 23:55:00")))+
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(26,33,1))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
all_stress.plot
ggsave(all_stress.plot,file="/Users/hannahaichelman/Dropbox/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/TankTempHobo_Lines_Combined_stress+recovery.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

#### Hobo Summary Stats ####
# Some summary stats for hogo logger temp data
describe(Tank1.1_var$Temp) # Control
# vars     n  mean   sd median trimmed  mad   min   max range  skew kurtosis se
# X1    1 14327 29.45 0.14  29.45   29.45 0.11 27.63 29.93   2.3 -1.84    20.35  0

describe(Tank2.1_var$Temp) # Low Var
# vars     n  mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 14327 29.15 0.61  28.99   29.11 0.66 28.06 30.56  2.51 0.43    -1.27 0.01

describe(Tank3.1_var$Temp) # Mod Var
# vars     n mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 14327 29.4 0.97  29.02    29.3 0.81 28.08 31.67  3.59 0.63    -1.04 0.01

describe(Tank4.2_var$Temp) # High Var
# vars     n  mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 14327 29.81 1.39  29.28    29.7 1.14 26.89 34.36  7.47 0.55    -0.85 0.01

# the describe() function looks like it is doing overall range, but want to find the average daily range,
# will use the method from dan barshis:

# Find overall mean during heat challenge and recovery periods:
library(plotrix)

all_stress_2 = all_stress %>%
  #dplyr::filter(DateTime < "2016-11-21 00:00:00" & DateTime > "2016-11-14 00:00:00") %>%# to look at just stress
  dplyr::filter(DateTime < "2016-12-11 23:55:00" & DateTime > "2016-11-25 00:00:00") %>%# to look at just recovery
  summarise_at(vars(Temp), list(mean_temp = mean, sd_temp = sd, se_temp = std.error), na.rm=TRUE)
all_stress_2

# mean_temp stress
# mean_temp   sd_temp    se_temp
# 1  31.81931 0.2255795 0.00282682

# mean_temp recovery
# mean_temp   sd_temp     se_temp
# 1  29.42072 0.4755736 0.003924869

#### DTV Function ####
# use Dan Barshis' function to calculate dtv for these loggers
# but need to drop NA's, because they introduce NA's into daily stats.

Tank1.1_var.nona = Tank1.1_var %>%
  drop_na(Temp)
Tank1.1_var.dailystats<-data.frame("DayRange"=tapply(Tank1.1_var.nona$Temp, Tank1.1_var.nona$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(Tank1.1_var.nona$Temp, Tank1.1_var.nona$Day, min),"DayMax"=tapply(Tank1.1_var.nona$Temp, Tank1.1_var.nona$Day, max), "DayMean"=tapply(Tank1.1_var.nona$Temp, Tank1.1_var.nona$Day, mean))
Tank1.1_var.dailystats <- tibble::rownames_to_column(Tank1.1_var.dailystats, "date")
Tank1.1_var.dailystats$treat <- "Control 1"
head(Tank1.1_var.dailystats)

Tank2.1_var.nona = Tank2.1_var %>%
  drop_na(Temp)
Tank2.1_var.dailystats<-data.frame("DayRange"=tapply(Tank2.1_var.nona$Temp, Tank2.1_var.nona$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(Tank2.1_var.nona$Temp, Tank2.1_var.nona$Day, min),"DayMax"=tapply(Tank2.1_var.nona$Temp, Tank2.1_var.nona$Day, max), "DayMean"=tapply(Tank2.1_var.nona$Temp, Tank2.1_var.nona$Day, mean))
Tank2.1_var.dailystats <- tibble::rownames_to_column(Tank2.1_var.dailystats, "date")
Tank2.1_var.dailystats$treat <- "Low Var"
head(Tank2.1_var.dailystats)

Tank3.1_var.nona = Tank3.1_var %>%
  drop_na(Temp)
Tank3.1_var.dailystats<-data.frame("DayRange"=tapply(Tank3.1_var.nona$Temp, Tank3.1_var.nona$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(Tank3.1_var.nona$Temp, Tank3.1_var.nona$Day, min),"DayMax"=tapply(Tank3.1_var.nona$Temp, Tank3.1_var.nona$Day, max), "DayMean"=tapply(Tank3.1_var.nona$Temp, Tank3.1_var.nona$Day, mean))
Tank3.1_var.dailystats <- tibble::rownames_to_column(Tank3.1_var.dailystats, "date")
Tank3.1_var.dailystats$treat <- "Mod Var"
head(Tank3.1_var.dailystats)

Tank4.2_var.nona = Tank4.2_var %>%
  drop_na(Temp)
Tank4.2_var.dailystats<-data.frame("DayRange"=tapply(Tank4.2_var.nona$Temp, Tank4.2_var.nona$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(Tank4.2_var.nona$Temp, Tank4.2_var.nona$Day, min),"DayMax"=tapply(Tank4.2_var.nona$Temp, Tank4.2_var.nona$Day, max), "DayMean"=tapply(Tank4.2_var.nona$Temp, Tank4.2_var.nona$Day, mean))
Tank4.2_var.dailystats <- tibble::rownames_to_column(Tank4.2_var.dailystats, "date")
Tank4.2_var.dailystats$treat <- "High Var"
head(Tank4.2_var.dailystats)

# combine daily variability stats for all loggers
library(Rmisc)

var.dailystats.all = rbind(Tank1.1_var.dailystats,Tank2.1_var.dailystats,Tank3.1_var.dailystats,Tank4.2_var.dailystats)
str(var.dailystats.all)

var.dailystats.all$treat = factor(var.dailystats.all$treat, levels = c("Control 1","Low Var","Mod Var","High Var"))

summarySE(data = var.dailystats.all, measurevar = "DayRange", groupvar = "treat")
#       treat  N DayRange        sd         se         ci
# 1 Control 1 50  0.43098 0.2996750 0.04238044 0.08516668
# 2   Low Var 50  1.80772 0.1717537 0.02428965 0.04881187
# 3   Mod Var 50  2.88386 0.1295098 0.01831544 0.03680626
# 4  High Var 50  4.09558 0.6225749 0.08804539 0.17693383

summarySE(data = var.dailystats.all, measurevar = "DayMean", groupvar = "treat")
#       treat  N  DayMean         sd         se         ci
# 1 Control 1 50 29.44626 0.08399914 0.01187927 0.02387229
# 2   Low Var 50 29.15265 0.09244450 0.01307363 0.02627244
# 3   Mod Var 50 29.40066 0.13128601 0.01856665 0.03731107
# 4  High Var 50 29.81221 0.28840981 0.04078731 0.08196516

summarySE(data = var.dailystats.all, measurevar = "DayMax", groupvar = "treat")
#       treat  N   DayMax        sd         se         ci
# 1 Control 1 50 29.59380 0.1112273 0.01572991 0.03161044
# 2   Low Var 50 30.20135 0.1932732 0.02733296 0.05492763
# 3   Mod Var 50 31.28463 0.1869496 0.02643867 0.05313049
# 4  High Var 50 32.40988 0.4872109 0.06890202 0.13846380

summarySE(data = var.dailystats.all, measurevar = "DayMin", groupvar = "treat")
#       treat  N   DayMin         sd         se         ci
# 1 Control 1 50 29.16282 0.31183470 0.04410009 0.08862244
# 2   Low Var 50 28.39363 0.09258051 0.01309286 0.02631109
# 3   Mod Var 50 28.40077 0.10868125 0.01536985 0.03088687
# 4  High Var 50 28.31430 0.37830137 0.05349989 0.10751206

#### HOBO DAILY STATS & PLOTS####
# daily range
aov.dtv=aov(DayRange~treat, data=var.dailystats.all)
summary(aov.dtv)
#             Df Sum Sq Mean Sq F value Pr(>F)
# treat         3  365.0  121.67   929.4 <2e-16 ***
# Residuals   196   25.7    0.13

par(mfrow=c(2,2))
plot(aov.dtv)
TukeyHSD(aov.dtv)
#                       diff       lwr      upr p adj
# Low Var-Control 1  1.37674 1.1892259 1.564254     0
# Mod Var-Control 1  2.45288 2.2653659 2.640394     0
# High Var-Control 1 3.66460 3.4770859 3.852114     0
# Mod Var-Low Var    1.07614 0.8886259 1.263654     0
# High Var-Low Var   2.28786 2.1003459 2.475374     0
# High Var-Mod Var   1.21172 1.0242059 1.399234     0

# daily mean
aov.mean=aov(DayMean~treat, data=var.dailystats.all)
summary(aov.mean)
#               Df Sum Sq Mean Sq F value Pr(>F)
# treat         3 11.102   3.701   127.6 <2e-16 ***
# Residuals   196  5.685   0.029

par(mfrow=c(2,2))
plot(aov.mean)
TukeyHSD(aov.mean)
#                           diff        lwr         upr    p adj
# Low Var-Control 1  -0.29360865 -0.3818689 -0.20534842 0.000000
# Mod Var-Control 1  -0.04559441 -0.1338546  0.04266581 0.539521
# High Var-Control 1  0.36595725  0.2776970  0.45421747 0.000000
# Mod Var-Low Var     0.24801424  0.1597540  0.33627446 0.000000
# High Var-Low Var    0.65956590  0.5713057  0.74782612 0.000000
# High Var-Mod Var    0.41155166  0.3232914  0.49981188 0.000000

# daily maximum
aov.max=aov(DayMax~treat, data=var.dailystats.all)
summary(aov.max)
#               Df Sum Sq Mean Sq F value Pr(>F)
# treat         3 230.94   76.98   956.1 <2e-16 ***
# Residuals   196  15.78    0.08

par(mfrow=c(2,2))
plot(aov.max)
TukeyHSD(aov.max)
#                         diff       lwr       upr p adj
# Low Var-Control 1  0.607548 0.4604983 0.7545977     0
# Mod Var-Control 1  1.690825 1.5437757 1.8378751     0
# High Var-Control 1 2.816077 2.6690276 2.9631270     0
# Mod Var-Low Var    1.083277 0.9362277 1.2303271     0
# High Var-Low Var   2.208529 2.0614796 2.3555790     0
# High Var-Mod Var   1.125252 0.9782022 1.2723016     0

# daily minimum
aov.min=aov(DayMin~treat, data=var.dailystats.all)
summary(aov.min)
#               Df Sum Sq Mean Sq F value Pr(>F)
# treat         3  23.83   7.942   121.8 <2e-16 ***
# Residuals   196  12.78   0.065

par(mfrow=c(2,2))
plot(aov.min)
TukeyHSD(aov.min)
#                           diff        lwr         upr     p adj
# Low Var-Control 1  -0.7691920 -0.9015049 -0.63687906 0.0000000
# Mod Var-Control 1  -0.7620546 -0.8943675 -0.62974166 0.0000000
# High Var-Control 1 -0.8485227 -0.9808356 -0.71620976 0.0000000
# Mod Var-Low Var     0.0071374 -0.1251755  0.13945034 0.9990226
# High Var-Low Var   -0.0793307 -0.2116436  0.05298224 0.4076884
# High Var-Mod Var   -0.0864681 -0.2187810  0.04584484 0.3299194

cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

# plot boxplots of  treatments - supplemental fig
tempPlot.dtv <- ggplot(var.dailystats.all, aes(x = treat, y = DayRange)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = treat)) +
  scale_color_manual(values = cols_treat_reds) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.7,
               aes(fill = treat))+
  scale_fill_manual(values = cols_treat_reds) + # for boxplot
  #xlab("Treatment")+
  ylab("Temperature (°C)")+
  #ylim(27,32)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")+
  ggtitle("A. Daily Variability")
tempPlot.dtv

tempPlot.mean <- ggplot(var.dailystats.all, aes(x = treat, y = DayMean)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = treat)) +
  scale_color_manual(values = cols_treat_reds) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.7,
               aes(fill = treat))+
  scale_fill_manual(values = cols_treat_reds) + # for boxplot
  xlab("Treatment")+
  ylab("Temperature (°C)")+
  #ylim(27,32)+
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")+
  ggtitle("B. Daily Mean")
tempPlot.mean

tempPlot.max <- ggplot(var.dailystats.all, aes(x = treat, y = DayMax)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = treat)) +
  scale_color_manual(values = cols_treat_reds) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.7,
               aes(fill = treat))+
  scale_fill_manual(values = cols_treat_reds) + # for boxplot
  xlab("Treatment")+
  ylab("Temperature (°C)")+
  #ylim(27,32)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("C. Daily Maximum")
tempPlot.max

tempPlot.min <- ggplot(var.dailystats.all, aes(x = treat, y = DayMin)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = treat)) +
  scale_color_manual(values = cols_treat_reds) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.7,
               aes(fill = treat))+
  scale_fill_manual(values = cols_treat_reds) + # for boxplot
  xlab("Treatment")+
  ylab("Temperature (°C)")+
  #ylim(27,32)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("D. Daily Minimum")
tempPlot.min


all.temp.plots = ggarrange(tempPlot.dtv,tempPlot.mean,tempPlot.max,tempPlot.min,
                           ncol = 2, nrow = 2)
ggsave(all.temp.plots, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/TankTemps/Tank_Hobo_Loggers/plots/DailyTemperatureBoxplots.pdf", width=6, height=6, units=c("in"), useDingbats=FALSE)

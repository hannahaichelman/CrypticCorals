##### METADATA #####
##  IN SITU TEMPERATURE DATA ##
#   project             PANAMA TVE
#   author              Hannah Aichelman (hannahaichelman@gmail.com)

#   function:
#   Visualize/analyze temperature data from HOBO loggers

#### Useful Notation Notes ####

#   xts
#       sample.xts['2007']            all of 2007
#       sample.xts['2007-03/']        March 2007 to the end of the data set
#       sample.xts['2007-03/2007']    March 2007 to the end of 2007
#       sample.xts['/']               the whole data set
#       sample.xts['/2007']           the beginning of the data through 2007
#       sample.xts['2007-01-03']      just the 3rd of January 2007

## Can skip straight to "Compiling All Logged Data in Longform" section and read in the files highlighted there to skip processing of individual sites


##### Required Packages #####
library(shiny)
library(plotly)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(xts)
library(zoo)
library(TTR)
library(scales)
library(ggpubr)
library(signal)
library(data.table)
library(ggridges)
library(Rmisc)


##### Set Color Palettes ####
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")
its2_cols_greens = c("C1" = "#edf8e9", "C3af" = "#238b45","C3" = "#a1d99b","D1" = "#00441b")

##### Importing and Formatting Data #####

# set wd
setwd("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Trimmed HOBO txt files")

# what's in the wd
list.files("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Trimmed HOBO txt files")

# make each HOBO file an xts/zoo object
# these lines will give you an object where the datetime is the rowname and the temperature value is in a column
# when you set the ISOdate parameters, you are meant to specify the exact start/end times of your dataset. however,
# even though all loggers across sites were logging every 15 minutes, they were not all coordinated to the minute. this makes time of day comparisons challenging.
# you can work around this by telling xts that your dataset starts at midnight and ends at midnight. this will shift all of your values to the same 15 minute intervals.
# important to note that xts recommends you use GMT or UTC when possible, as dealing with timezones is difficult;
# additionally, there is no simple way to account for DST. this means that after spring forward, data are consistently misaligned by one hour.
# not a problem for large scale patterns but be mindful of the impact on more granular analyses

# main loggers (one logger per site)
Cayo.OR3.arr3 <- xts(zoo(read.table("Cayo.OR3.arr3.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
Drago.OR4 <- xts(zoo(read.table("Drago.OR4.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
Punta.IR1 <- xts(zoo(read.table("Punta.IR1.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2 <- xts(zoo(read.table("STRI.IR2.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
Cristo.IR3.arr1 <- xts(zoo(read.table("Cristo.IR3.arr1.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
# fix cristobal temp to degrees C because the logger was deployed in F
Cristo.IR3.arr1 <- ((Cristo.IR3.arr1 - 32) * (5/9))

# STRI array loggers
STRI.IR2.arr4 <- xts(zoo(read.table("STRI.IR2.arr4.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr5 <- xts(zoo(read.table("STRI.IR2.arr5.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr6 <- xts(zoo(read.table("STRI.IR2.arr6.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr7 <- xts(zoo(read.table("STRI.IR2.arr7.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr10 <- xts(zoo(read.table("STRI.IR2.arr10.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr12 <- xts(zoo(read.table("STRI.IR2.arr12.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,8,14,0,0,0), "15 min", tz="GMT")))


##### Monthly Max, Min, and Mean #####

# i. maximum temperature by month #

# read in the below csv to skip this section of code
monthly.mmm <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/MonthlyTempParams.csv")
###

# (step 1) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.mmax <- apply.monthly(Cayo.OR3.arr3['/'], FUN = max)
Drago.OR4.mmax <- apply.monthly(Drago.OR4['/'], FUN = max)
Punta.IR1.mmax <- apply.monthly(Punta.IR1['/'], FUN = max)
STRI.IR2.mmax <- apply.monthly(STRI.IR2['/'], FUN = max)
Cristo.IR3.arr1.mmax <- apply.monthly(Cristo.IR3.arr1['/'], FUN = max)

# (step 2) merge the sites together
#          can only merge two objects at a time
merge <- merge(Cayo.OR3.arr3.mmax,Drago.OR4.mmax)
merge <- merge(merge, Punta.IR1.mmax)
merge <- merge(merge,STRI.IR2.mmax)
merged <- merge(merge, Cristo.IR3.arr1.mmax)

# (step 3) rename the columns to keep identifying info
colnames(merged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 4) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
merged <- as.data.frame(merged)
merged <- setDT(merged, keep.rownames = T)
colnames(merged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) melt to put the data in longform
melt.mmax <- melt(merged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.mmax) <- c('datetime', 'logger', 'mmax')
head(melt.mmax)

# repeat steps 1-5 for all parameters of interest...

# ii. minimum temperature by month #
Cayo.OR3.arr3.mmin<- apply.monthly(Cayo.OR3.arr3['/'], FUN = min)
Drago.OR4.mmin <- apply.monthly(Drago.OR4['/'], FUN = min)
Punta.IR1.mmin <- apply.monthly(Punta.IR1['/'], FUN = min)
STRI.IR2.mmin <- apply.monthly(STRI.IR2['/'], FUN = min)
Cristo.IR3.arr1.mmin <- apply.monthly(Cristo.IR3.arr1['/'], FUN = min)
merge <- merge(Cayo.OR3.arr3.mmin,Drago.OR4.mmin)
merge <- merge(merge, Punta.IR1.mmin)
merge <- merge(merge,STRI.IR2.mmin)
merged <- merge(merge, Cristo.IR3.arr1.mmin)
colnames(merged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged <- as.data.frame(merged)
merged <- setDT(merged, keep.rownames = T)
colnames(merged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.mmin <- melt(merged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.mmin) <- c('datetime', 'logger', 'mmin')
head(melt.mmin)

# iii. mean temperature by month #
Cayo.OR3.arr3.mmean<- apply.monthly(Cayo.OR3.arr3['/'], FUN = mean)
Drago.OR4.mmean <- apply.monthly(Drago.OR4['/'], FUN = mean)
Punta.IR1.mmean <- apply.monthly(Punta.IR1['/'], FUN = mean)
STRI.IR2.mmean <- apply.monthly(STRI.IR2['/'], FUN = mean)
Cristo.IR3.arr1.mmean <- apply.monthly(Cristo.IR3.arr1['/'], FUN = mean)
merge <- merge(Cayo.OR3.arr3.mmean,Drago.OR4.mmean)
merge <- merge(merge, Punta.IR1.mmean)
merge <- merge(merge,STRI.IR2.mmean)
merged <- merge(merge, Cristo.IR3.arr1.mmean)
colnames(merged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged <- as.data.frame(merged)
merged <- setDT(merged, keep.rownames = T)
colnames(merged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.mmean <- melt(merged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.mmean) <- c('datetime', 'logger', 'mmean')
head(melt.mmean)

# merge the max, min, and mean and melt
m.maxmin <- merge(melt.mmax, melt.mmin)
m.maxminmean <- merge(m.maxmin, melt.mmean)
melt.mmmm <- melt(m.maxminmean, id.vars = c('datetime','logger'), measure.vars = c('mmax','mmin','mmean'))
head(melt.mmmm)

# it's easier to do use setDT to get "typical day" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of day
monthly.mmm<-setDT(melt.mmmm)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
monthly.mmm$range <- (monthly.mmm$max - monthly.mmm$min)
monthly.mmm$logger <- factor(monthly.mmm$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(monthly.mmm)

#add in separate date column, since time doesn't matter for weekly data
monthly.mmm$date <- as.Date(monthly.mmm$datetime,"%Y-%m-%d" )

# add in separate month column
head(monthly.mmm)

monthly.mmm$month <- ifelse(monthly.mmm$date == '2015-06-30', 'jun',
                     ifelse(monthly.mmm$date == '2015-07-31', 'jul',
                     ifelse(monthly.mmm$date == '2015-08-31', 'aug',
                     ifelse(monthly.mmm$date == '2015-09-30', 'sep',
                     ifelse(monthly.mmm$date == '2015-10-31', 'oct',
                     ifelse(monthly.mmm$date == '2015-11-30', 'nov',
                     ifelse(monthly.mmm$date == '2015-12-31', 'dec',
                     ifelse(monthly.mmm$date == '2016-01-31', 'jan',
                     ifelse(monthly.mmm$date == '2016-02-29', 'feb',
                     ifelse(monthly.mmm$date == '2016-03-31', 'mar',
                     ifelse(monthly.mmm$date == '2016-04-30', 'apr',
                     ifelse(monthly.mmm$date == '2016-05-31', 'may',
                     ifelse(monthly.mmm$date == '2016-06-30', 'jun',
                     ifelse(monthly.mmm$date == '2016-07-31', 'jul',
                     ifelse(monthly.mmm$date == '2016-08-12', 'aug',
                     'NA')))))))))))))))

#write.csv(monthly.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/MonthlyTempParams.csv", row.names = FALSE)

aggregate(range ~ logger, data = monthly.mmm, FUN = "mean")
# logger    range
# 1       Punta.IR1 2.422933
# 2        STRI.IR2 1.932067
# 3 Cristo.IR3.arr1 3.525704
# 4   Cayo.OR3.arr3 2.916000
# 5       Drago.OR4 2.631533


##### Weekly Max, Min, Mean, and Range #####

# i. maximum temperature by week #

# read in this csv to skip running this section of code:
weekly.mmm <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/WeeklyTempParams.csv")
###

# (step 1) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.wmax <- apply.weekly(Cayo.OR3.arr3['/'], FUN = max)
Drago.OR4.wmax <- apply.weekly(Drago.OR4['/'], FUN = max)
Punta.IR1.wmax <- apply.weekly(Punta.IR1['/'], FUN = max)
STRI.IR2.wmax <- apply.weekly(STRI.IR2['/'], FUN = max)
Cristo.IR3.arr1.wmax <- apply.weekly(Cristo.IR3.arr1['/'], FUN = max)

# (step 2) merge the sites together
#          can only merge two objects at a time
weekmerge <- merge(Cayo.OR3.arr3.wmax,Drago.OR4.wmax)
weekmerge <- merge(weekmerge, Punta.IR1.wmax)
weekmerge <- merge(weekmerge,STRI.IR2.wmax)
weekmerged <- merge(weekmerge, Cristo.IR3.arr1.wmax)

# (step 3) rename the columns to keep identifying info
colnames(weekmerged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 4) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
weekmerged <- as.data.frame(weekmerged)
weekmerged <- setDT(weekmerged, keep.rownames = T)
colnames(weekmerged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) melt to put the data in longform
melt.wmax <- melt(weekmerged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmax) <- c('datetime', 'logger', 'wmax')
head(melt.wmax)

# repeat steps 1-5 for all parameters of interest...

# ii. minimum temperature by week #
Cayo.OR3.arr3.wmin<- apply.weekly(Cayo.OR3.arr3['/'], FUN = min)
Drago.OR4.wmin <- apply.weekly(Drago.OR4['/'], FUN = min)
Punta.IR1.wmin <- apply.weekly(Punta.IR1['/'], FUN = min)
STRI.IR2.wmin <- apply.weekly(STRI.IR2['/'], FUN = min)
Cristo.IR3.arr1.wmin <- apply.weekly(Cristo.IR3.arr1['/'], FUN = min)
weekmerge <- merge(Cayo.OR3.arr3.wmin,Drago.OR4.wmin)
weekmerge <- merge(weekmerge, Punta.IR1.wmin)
weekmerge <- merge(weekmerge,STRI.IR2.wmin)
weekmerged <- merge(weekmerge, Cristo.IR3.arr1.wmin)
colnames(weekmerged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
weekmerged <- as.data.frame(weekmerged)
weekmerged <- setDT(weekmerged, keep.rownames = T)
colnames(weekmerged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.wmin <- melt(weekmerged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmin) <- c('datetime', 'logger', 'wmin')
head(melt.wmin)

# iii. mean temperature by week #
Cayo.OR3.arr3.wmean<- apply.weekly(Cayo.OR3.arr3['/'], FUN = mean)
Drago.OR4.wmean <- apply.weekly(Drago.OR4['/'], FUN = mean)
Punta.IR1.wmean <- apply.weekly(Punta.IR1['/'], FUN = mean)
STRI.IR2.wmean <- apply.weekly(STRI.IR2['/'], FUN = mean)
Cristo.IR3.arr1.wmean <- apply.weekly(Cristo.IR3.arr1['/'], FUN = mean)
weekmerge <- merge(Cayo.OR3.arr3.wmean,Drago.OR4.wmean)
weekmerge <- merge(weekmerge, Punta.IR1.wmean)
weekmerge <- merge(weekmerge,STRI.IR2.wmean)
weekmerged <- merge(weekmerge, Cristo.IR3.arr1.wmean)
colnames(weekmerged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
weekmerged <- as.data.frame(weekmerged)
weekmerged <- setDT(weekmerged, keep.rownames = T)
colnames(weekmerged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.wmean <- melt(weekmerged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmean) <- c('datetime', 'logger', 'wmean')
head(melt.wmean)

# merge the max, min, and mean and melt
w.maxmin <- merge(melt.wmax, melt.wmin)
w.maxminmean <- merge(w.maxmin, melt.wmean)
melt.wmmm <- melt(w.maxminmean, id.vars = c('datetime','logger'), measure.vars = c('wmax','wmin','wmean'))
head(melt.wmmm)

# it's easier to do use setDT to get "typical day" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of day
weekly.mmm<-setDT(melt.wmmm)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
weekly.mmm$range <- (weekly.mmm$max - weekly.mmm$min)
weekly.mmm$logger <- factor(weekly.mmm$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(weekly.mmm)

#add in separate date column, since time doesn't matter for weekly data
weekly.mmm$date <- as.Date(weekly.mmm$datetime,"%Y-%m-%d" )

# add in separate month column
head(weekly.mmm)

#write.csv(weekly.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/WeeklyTempParams.csv", row.names = FALSE)

aggregate(range ~ logger, data = weekly.mmm, FUN = "mean")
# logger    range
# 1       Punta.IR1 1.485613
# 2        STRI.IR2 1.100710
# 3 Cristo.IR3.arr1 2.372599
# 4   Cayo.OR3.arr3 1.932871
# 5       Drago.OR4 1.683935


#### Compiling All Logged Data in Longform ####

# read in this csv file to skip running this code for all temp data
master.melt <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllTempData.csv")

###

# merge the sites together
# can only merge two objects at a time
merge <- merge(Cayo.OR3.arr3,Drago.OR4)
merge <- merge(merge, Punta.IR1)
merge <- merge(merge,STRI.IR2)
merged <- merge(merge, Cristo.IR3.arr1)
head(merged)

# rename the columns to keep identifying info
colnames(merged) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged <- as.data.frame(merged)
# merged is still an xts/zoo object. make it a data.table and keep the datetimes as rownames

merged <- setDT(merged, keep.rownames = T)
colnames(merged) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
head(merged)

# melt, format, add in separate date and time columns
master.melt<- melt(merged, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(master.melt) <- c('datetime','logger','temp')
head(master.melt)
# master.melt is now every measurement you have for each site's "main logger" from the period specified in the initial import as zoo/xts files

# when you converted master.melt to a data.table, the datetime became a character string. this is useful because we can pull date and time out individually.
master.melt$time <- substr(master.melt$datetime, 12, 19)
master.melt$date <- as.Date(master.melt$datetime,"%Y-%m-%d" )

# but now put them back in POSIX format so R knows they're dates
# note that when you do this, the time column will be arbitrarily assigned all the same date (current date)
master.melt$datetime <- as.POSIXct(strptime(master.melt$datetime, "%Y-%m-%d %H:%M:%S"))
master.melt$time <- as.POSIXct(strptime(master.melt$time, "%H:%M:%S"))
head(master.melt)
tail(master.melt)

#write.csv(master.melt, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllTempData.csv", row.names = FALSE)

## Summary Stats

aggregate(temp ~ logger, data = master.melt, FUN = "min")
# logger   temp
# 1   Cayo.OR3.arr3 27.038
# 2       Drago.OR4 26.769
# 3       Punta.IR1 27.038
# 4        STRI.IR2 27.038
# 5 Cristo.IR3.arr1 27.235

aggregate(temp ~ logger, data = master.melt, FUN = "max")
# logger     temp
# 1   Cayo.OR3.arr3 32.43300
# 2       Drago.OR4 31.66300
# 3       Punta.IR1 32.15000
# 4        STRI.IR2 31.58600
# 5 Cristo.IR3.arr1 33.54778

aggregate(temp ~ logger, data = master.melt, FUN = "mean")
# logger     temp
# 1   Cayo.OR3.arr3 29.45465
# 2       Drago.OR4 29.13399
# 3       Punta.IR1 29.70488
# 4        STRI.IR2 29.31704
# 5 Cristo.IR3.arr1 29.86366

#### Seasonal Parameters ####

# read in these csv files to skip running this section of code for all data, daytime data only, and nighttime data only, can skip straight to plots

dailys <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/SeasonalDailyTempData.csv")

###

# separate into three-month chunks
# in Panama, the rainy season typically extends from about May to December. Jan, Feb, and March are relatively dry; rain starts to pick up in April and May, is in full swing June, July, August, September, October, and November; starts to taper in December
# so in this blocking, Jan-Mar = dry, Apr-Jun = ramp up into rain, Jul-Sep = full on rain, and Oct-Dec = ramp down to dry
# also nice because August was my colony collection month, so I can get the month preceding, month of, and month after


Jul.Aug.Sep <-
  master.melt %>%
  dplyr::filter(datetime >= '2015-07-1 00:00:00', datetime <= '2015-09-30 11:59:59' | datetime >= '2016-07-1 00:00:00', datetime <= '2016-09-30 11:59:59')
Jul.Aug.Sep$season <- 'julaugsep'

Oct.Nov.Dec <-
  master.melt %>%
  dplyr::filter(datetime >= '2015-10-1 00:00:00', datetime <= '2015-12-31 11:59:59')
Oct.Nov.Dec$season <- 'octnovdec'

Jan.Feb.Mar <-
  master.melt %>%
  dplyr::filter(datetime >= '2016-01-1 00:00:00', datetime <= '2016-03-31 11:59:59')
Jan.Feb.Mar$season <- 'janfebmar'

Apr.May.Jun <-
  master.melt %>%
  dplyr::filter(datetime >= '2015-06-1 00:00:00', datetime <= '2015-06-30 11:59:59' | datetime >= '2016-04-1 00:00:00', datetime <= '2016-06-30 11:59:59')
Apr.May.Jun$season <- 'aprmayjun'

seas <- rbind(Oct.Nov.Dec,Jan.Feb.Mar,Apr.May.Jun,Jul.Aug.Sep)
head(seas)


dailys<-setDT(master.melt)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
dailys$range <- (dailys$max - dailys$min)
dailys$season <- ifelse(dailys$date >= '2015-10-1' & dailys$date <= '2015-12-31', 'octnovdec',
                 ifelse(dailys$date >= '2016-01-1' & dailys$date <= '2016-03-31', 'janfebmar',
                 ifelse(dailys$date >= '2016-04-1' & dailys$date <= '2016-06-30', 'aprmayjun',
                 ifelse(dailys$date >= '2015-06-1' & dailys$date <= '2015-06-30', 'aprmayjun',
                  'julaugsep'))))

dailys$season <- factor(dailys$season, levels = c('janfebmar','aprmayjun','julaugsep','octnovdec'))
dailys$logger <- factor(dailys$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))

head(dailys)

#write.csv(dailys, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/SeasonalDailyTempData.csv", row.names = FALSE)

# rank sites by daily variability by season
aggregate(range ~ logger + season, data = dailys, FUN = "mean")
# Order for all seasons : Cristobal, Cayo de Agua, Drago Mar, Punta Donato, STRI


#PLOT seasonal DTV and mean temperature data
str(dailys)
dailys$season = as.factor(dailys$season)
dailys$season <- factor(dailys$season, levels = c('janfebmar','aprmayjun','julaugsep','octnovdec'))
levels(dailys$season) <- c('JanFebMar','AprMayJun','JulAugSep','OctNovDec')

dailys$logger = as.factor(dailys$logger)
levels(dailys$logger) <- c("PD", "SP", "CI", "CA", "Drago.OR4")

dailys_plot = dailys %>%
  dplyr::filter(logger!="Drago.OR4")

#plot mean temperature by season
mean_boxplot_season <- ggplot(dailys_plot, aes(x=logger, y=mean)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Mean Temperature (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~season, nrow = 1, ncol = 4)
mean_boxplot_season

#plot mean temperature by season
dtv_boxplot_season <- ggplot(dailys_plot, aes(x=logger, y=range)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Daily Temperature Range (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~season, nrow = 1, ncol = 4)
dtv_boxplot_season


boxplots_season = ggarrange(mean_boxplot_season, dtv_boxplot_season, nrow = 2, ncol = 1)

ggsave(boxplots_season, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/SeasonTemps_boxplot.pdf", width=8, height=6, units=c("in"), useDingbats=FALSE)

# rank sites by daily variability by season
aggregate(range ~ logger + season, data = dailys_plot, FUN = "mean")
# regardless of season, the order of sites is:
#CI > CA > PD > SP

# logger    season     range
# 1      PD JanFebMar 0.6784725
# 2      SP JanFebMar 0.6061538
# 3      CI JanFebMar 1.2898107
# 4      CA JanFebMar 0.9969670
# 5      PD AprMayJun 0.6939375
# 6      SP AprMayJun 0.5237411
# 7      CI AprMayJun 1.2052083
# 8      CA AprMayJun 0.9999732
# 9      PD JulAugSep 0.7270073
# 10     SP JulAugSep 0.5229635
# 11     CI JulAugSep 1.0572749
# 12     CA JulAugSep 0.9770438
# 13     PD OctNovDec 0.7241957
# 14     SP OctNovDec 0.5255652
# 15     CI OctNovDec 1.1569143
# 16     CA OctNovDec 0.9303804

# rank sites by mean by season
aggregate(mean ~ logger + season, data = dailys_plot, FUN = "mean")
# regardless of season,
# janfebmar: CI > PD > CA > SP
# aprmayjun: CI > PD > CA > SP
# julaugsep: CI > PD > CA > SP
# octnovdec: CI > PD > CA > SP

m1 = lm(mean ~ season + logger, data = dailys_plot)
summary(m1)
anova(m1)
# Response: range
# Df  Sum Sq Mean Sq  F value  Pr(>F)
# season       3   1.226   0.409   2.9021 0.03376 *
# logger       3 100.021  33.340 236.8259 < 2e-16 ***
# Residuals 1721 242.282   0.141

# Response: mean
# Df Sum Sq Mean Sq F value    Pr(>F)
# season       3 411.88 137.292 325.092 < 2.2e-16 ***
# logger       3  78.15  26.050  61.684 < 2.2e-16 ***
# Residuals 1721 726.81   0.422

# try looking at just the month before collection
collection <- master.melt %>%
  dplyr::filter(datetime >= '2016-07-12 00:00:00', datetime <= '2016-08-12 00:00:00')

collections<-setDT(collection)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
collections$range <- (collections$max - collections$min)
head(collections)

aggregate(mean ~ logger, data = collections, FUN = "mean")
# logger     mean
# 1   Cayo.OR3.arr3 29.12143
# 2       Drago.OR4 28.89048
# 3       Punta.IR1 29.41228
# 4        STRI.IR2 28.99112
# 5 Cristo.IR3.arr1 29.32083


aggregate(range ~ logger, data = collections, FUN = "mean")
# logger     range
# 1   Cayo.OR3.arr3 0.9967813
# 2       Drago.OR4 0.7754687
# 3       Punta.IR1 0.7657813
# 4        STRI.IR2 0.4335625
# 5 Cristo.IR3.arr1 1.0529167

str(collections)
collections$logger = as.factor(collections$logger)
levels(collections$logger) <- c("CA", "Drago.OR4", "PD", "SP", "CI")

collections_plot = collections %>%
  dplyr::filter(logger!="Drago.OR4")

#plot mean temperature before collection
mean_boxplot_collection <- ggplot(collections_plot, aes(x=logger, y=mean)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Mean Temperature (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none")
mean_boxplot_collection

#plot mean temperature by season
dtv_boxplot_collection <- ggplot(collections_plot, aes(x=logger, y=range)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Daily Temperature Range (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none")
dtv_boxplot_collection


boxplots_collection = ggarrange(mean_boxplot_collection, dtv_boxplot_collection, nrow = 1, ncol = 2)

ggsave(boxplots_collection, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/PreCollectionTemps_boxplot.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


#### Daily Temperature Variation ####

# read in this csv to skip running this section of code
daily.mmm <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/DailyTempRangeData.csv")
###

# it's easier to do use setDT to get "typical day" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of day
daily.mmm<-setDT(master.melt)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
# add in range
daily.mmm$range <- (daily.mmm$max - daily.mmm$min)
daily.mmm$logger <- factor(daily.mmm$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(daily.mmm)

#remove 0's from short last day of deployment
daily.mmm = daily.mmm %>%
  dplyr::filter(range>0.01)

#write.csv(daily.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/DailyTempRangeData.csv", row.names = FALSE)

daily.mmm = daily.mmm %>%
  dplyr::filter(logger!="Drago.OR4")

# rank sites by daily variability overall
aggregate(range ~ logger, data = daily.mmm, FUN = "mean")

# logger     range
# 1       Punta.IR1 0.7092529
# 2        STRI.IR2 0.5424988
# 3 Cristo.IR3.arr1 1.1685357
# 4   Cayo.OR3.arr3 0.9795151

aggregate(range ~ logger, data = daily.mmm, FUN = "max")

# logger    range
# 1       Punta.IR1 1.399000
# 2        STRI.IR2 1.263000
# 3 Cristo.IR3.arr1 3.167222
# 4   Cayo.OR3.arr3 2.411000

summarySE(data = daily.mmm, groupvar = "logger", measurevar = "mean")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "min")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "max")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "range")



#### 90% Daily Thermal Range ####

# Calculate 90% quantile of the daily temperature range (Kenkel et al., Ecology)
# an estimate of the magnitude of high-frequency temp fluctuations
# we can use the quantile() function to calculate this

# daily variability logger data
daily.mmm %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
# Punta.IR1           1.02
# STRI.IR2            0.822
# Cristo.IR3.arr1     1.93
# Cayo.OR3.arr3       1.49

# weekly variability logger data, all data points included to calculate weekly variability
head(weekly.mmm)

weekly.mmm %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
# Punta.IR1            1.90
# STRI.IR2             1.41
# Cristo.IR3.arr1      3.14
# Cayo.OR3.arr3        2.56
# Drago.OR4            2.27

#### Figure 1 Boxplots ####

daily.mmm.plot = daily.mmm %>%
  dplyr::filter(logger != "Drago.OR4")

daily.mmm.plot$logger = droplevels(daily.mmm.plot$logger)

str(daily.mmm.plot)

daily.mmm.plot$logger = as.factor(daily.mmm.plot$logger)
levels(daily.mmm.plot$logger) <- c("PD", "SP", "CI", "CA") # double check this
daily.mmm.plot$logger <- factor(daily.mmm.plot$logger, levels = c("PD", "SP", "CI", "CA"))

# daily mean as boxplot

mean_boxplot <- ggplot(daily.mmm.plot, aes(x=logger, y=mean)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Daily Mean Temperature (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none")
mean_boxplot

ggsave(mean_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/DailyMeanTemp_boxplot.pdf", width=2.5, height=3.7, units=c("in"), useDingbats=FALSE)

# daily range as boxplot
dtv_boxplot <- ggplot(daily.mmm.plot, aes(x=logger, y=range)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.99,
              aes(color = logger)) +
  scale_color_manual(values = cols_site) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = logger))+
  scale_fill_manual(values = cols_site) + # for boxplot
  ylab("Daily Temperature Range (°C)") +
  xlab("Site") +
  #ylim(0,3.2) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none")
dtv_boxplot

ggsave(dtv_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/DailyTempRange_boxplot.pdf", width=2.5, height=3.7, units=c("in"), useDingbats=FALSE)

# stats by logger for mean and dtv
aov1 = aov(range ~ logger, data = daily.mmm.plot)
summary(aov1)
# Df Sum Sq Mean Sq F value Pr(>F)
# logger         3  100.2   33.42   239.1 <2e-16 ***
# Residuals   1720  240.4    0.14

TukeyHSD(aov1)
# $logger
# diff        lwr        upr p adj
# SP-PD -0.1667541 -0.2322429 -0.1012652     0
# CI-PD  0.4592828  0.3937940  0.5247716     0
# CA-PD  0.2702622  0.2047734  0.3357510     0
# CI-SP  0.6260369  0.5605481  0.6915257     0
# CA-SP  0.4370162  0.3715274  0.5025051     0
# CA-CI -0.1890206 -0.2545094 -0.1235318     0

aov2 = aov(mean ~ logger, data = daily.mmm.plot)
summary(aov2)
# Df Sum Sq Mean Sq F value Pr(>F)
# logger         3   77.9  25.977   39.26 <2e-16 ***
# Residuals   1720 1138.0   0.662

TukeyHSD(aov2)
# $logger
# diff          lwr        upr     p adj
# SP-PD -0.3878341 -0.530322481 -0.2453458 0.0000000
# CI-PD  0.1587846  0.016296206  0.3012729 0.0218856
# CA-PD -0.2502358 -0.392724163 -0.1077475 0.0000398
# CI-SP  0.5466187  0.404130336  0.6891070 0.0000000
# CA-SP  0.1375983 -0.004890032  0.2800867 0.0628787
# CA-CI -0.4090204 -0.551508719 -0.2665320 0.0000000

# find average parameters for each site
daily.mmm.plot %>%
  group_by(logger) %>%
  summarise_each(funs(mean))

# logger date         max   min  mean range
# PD     2016-01-11  30.1  29.4  29.7 0.709
# SP     2016-01-11  29.6  29.1  29.3 0.542
# CI     2016-01-11  30.5  29.4  29.9 1.17
# CA     2016-01-11  30.0  29.0  29.5 0.980

# compare to collections
collections_plot %>%
  group_by(logger) %>%
  summarise_each(funs(mean))

# logger date         max   min  mean range
# PD     2016-07-27  29.8  29.0  29.4 0.766
# SP     2016-07-27  29.2  28.8  29.0 0.434
# CI     2016-07-27  29.9  28.9  29.3 1.05
# CA     2016-07-27  29.6  28.6  29.1 0.997

#### Extra Plots ####

# i. daily range distribution by site

str(daily.mmm)
daily.mmm$date = as.Date(daily.mmm$date, format = "%Y-%m-%d")

daily.mmm$logger = as.factor(daily.mmm$logger)
daily.mmm$logger <- factor(daily.mmm$logger, levels = c("STRI.IR2", "Punta.IR1", "Drago.OR4",
                                                        "Cayo.OR3.arr3","Cristo.IR3.arr1"))
levels(daily.mmm$logger) <- c("SP", "PD", "Drago.OR4", "CA","CI")

cols_site = c("CA" = "#80cdc1", "CI"="#543005", "Drago.OR4"="black",
              "PD"="#bf812d", "SP"="#dfc27d")


range.plot = daily.mmm %>%
  dplyr::filter(logger != "Drago.OR4") %>%
  ggplot(aes(x = range, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.5, rel_min_height = .5, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_manual(values= cols_site)+
  stat_density_ridges(quantile_lines = T, scale = 1.5)+
  scale_x_continuous(breaks = seq(0,3,.5))+
  theme_ridges(center = T)+
  theme(legend.position = "none")
range.plot
ggsave(range.plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/DensityPlot_range.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

daily.mmm %>%
  ggplot(aes(x = min, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.5, rel_min_height = .5, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_manual(values= cols_site)+
  stat_density_ridges(quantile_lines = T, scale = 1.5)+
  scale_x_continuous(breaks = seq(26,34,1))+
  theme_ridges(center = T)

daily.mmm %>%
  ggplot(aes(x = max, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.5, rel_min_height = .5, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_manual(values= cols_site)+
  stat_density_ridges(quantile_lines = T, scale = 1.5)+
  scale_x_continuous(breaks = seq(26,34,1))+
  theme_ridges(center = T)


#### Compare with external loggers around STRI ####

# Data from Maggie Johnson
mj.temp<-read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/ExternalLabs_Hobo_Loggers/MaggieJohnson/FieldParameters_CRandHP_17-18.csv")
mj.temp$Date<-strptime(mj.temp$Date, format="%m/%d/%y %H:%M")
colnames(mj.temp)[2] <- "datetime"
mj.temp$datetime_ct <- as.POSIXct(mj.temp$datetime, format="%Y-%m-%dT%H:%M:%S")
mj.temp$Day<-format(mj.temp$datetime,"%D")
#mj.temp$Day<-as.POSIXct(mj.temp$Day, format="%m/%d/%y")
mj.temp$Site <- factor(mj.temp$Site)

str(mj.temp)

# plot temperature data through time

mj.temp %>%
  ggplot(aes(x = datetime_ct, y = Temp, color = Site, group = Site))+
  geom_line(lwd = 1)+
  scale_color_manual(name = "Site",
                     labels = c("Cayo_Roldan","Hospital_Point"),
                     values = c("coral1","darkslategray4")) +
  ggtitle("Maggie Johnson Data") +
  theme_classic()

# Calculate daily ranges for Maggie's sites
mj.temp.cayo = mj.temp %>%
  subset(Site=="Cayo")

mj.temp.hosp = mj.temp %>%
  subset(Site=="Hosp")

mj.temp.cayo.daily<-data.frame("DayRange"=tapply(mj.temp.cayo$Temp, mj.temp.cayo$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(mj.temp.cayo$Day, mj.temp.cayo$Day, min),"DayMax"=tapply(mj.temp.cayo$Temp, mj.temp.cayo$Day, max), "DayMean"=tapply(mj.temp.cayo$Temp, mj.temp.cayo$Day, mean))
mj.temp.cayo.daily <- tibble::rownames_to_column(mj.temp.cayo.daily, "date")
mj.temp.cayo.daily$site <- "Cayo_Roldan"

mj.temp.hosp.daily<-data.frame("DayRange"=tapply(mj.temp.hosp$Temp, mj.temp.hosp$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(mj.temp.hosp$Day, mj.temp.hosp$Day, min),"DayMax"=tapply(mj.temp.hosp$Temp, mj.temp.hosp$Day, max), "DayMean"=tapply(mj.temp.hosp$Temp, mj.temp.hosp$Day, mean))
mj.temp.hosp.daily <- tibble::rownames_to_column(mj.temp.hosp.daily, "date")
mj.temp.hosp.daily$site <- "Hospital_Point"

# Data from Noelle Lucey
nl.temp<-read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/ExternalLabs_Hobo_Loggers/NoelleLucey/BTD_Temp_cleaned.csv")
nl.temp$datetime<-paste(nl.temp$date,nl.temp$time,sep=" ")
nl.temp$datetime<-strptime(nl.temp$datetime, format="%m/%d/%y %H:%M:%S")
nl.temp$datetime_ct <- as.POSIXct(nl.temp$datetime, format="%Y-%m-%dT%H:%M:%S")

nl.temp$Day<-format(nl.temp$datetime,"%D")
nl.temp$date<-as.POSIXct(nl.temp$date, format="%m/%d/%y")

str(nl.temp)

nl.temp$site <- factor(nl.temp$site)

# noelle says Tranquilo Bay is new data, so remove any values of 36 or 23
nl.temp = nl.temp %>%
  dplyr::filter(temp <= 36, temp >= 23)
summary(nl.temp)

nl.temp %>%
  ggplot(aes(x = datetime_ct, y = temp, color = site, group = site))+
  geom_line()+
  scale_color_manual(name = "Site",
                     labels = c("Bella_Vista","Cayo_Wilson","Punta_Caracol","Tranquilo_Bay"),
                     values = c("cornflowerblue","indianred1","darkgray","seagreen")) +
  ggtitle("Noelle Lucey Data") +

  theme_classic()

# Other data from Noelle Lucey for Crawl Cay
# Add in Crawl Cay data from Noelle to see how it fits:
library(readxl)
crawl_cay = read_excel("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/ExternalLabs_Hobo_Loggers/NoelleLucey/temp_data.xlsx",
                       sheet = "Crawl-Cay")
head(crawl_cay)
crawl_cay$DateTime<-strptime(crawl_cay$DateTime, format="%Y-%m-%d %H:%M:%S")

crawl_cay$Day<-format(crawl_cay$DateTime,"%D")

str(crawl_cay)

# Calculate daily ranges for Noelle's sites
nl.temp.vista = nl.temp %>%
  subset(site=="Bella_Vista")

nl.temp.wilson = nl.temp %>%
  subset(site=="Cayo_Wilson")

nl.temp.caracol = nl.temp %>%
  subset(site=="Punta_Caracol")

nl.temp.bay.spring = nl.temp %>%
  subset(site=="Tranquilo_Bay") %>%
  dplyr::filter(Day > "04/03/19" & Day < "06/10/19")

nl.temp.bay.summer = nl.temp %>%
  subset(site=="Tranquilo_Bay") %>%
  dplyr::filter(Day > "06/20/19" & Day < "09/25/19")

p1 = nl.temp.vista %>%
  ggplot(aes(x = datetime_ct, y = temp))+
  geom_line(color = 'coral', lwd=1.5) +
  ggtitle("Noelle Lucey - Bella Vista") +
  ylim(26,34) +
  theme_classic()

p2 = nl.temp.wilson %>%
  ggplot(aes(x = datetime_ct, y = temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Noelle Lucey - Cayo Wilson") +
  ylim(26,34) +
  theme_classic()

p3 = nl.temp.caracol %>%
  ggplot(aes(x = datetime_ct, y = temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Noelle Lucey - Punta Caracol") +
  ylim(26,34) +
  theme_classic()

p4 = nl.temp.bay.spring %>%
  ggplot(aes(x = datetime_ct, y = temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Noelle Lucey - Tranquilo Bay - Spring") +
  ylim(26,34) +
  theme_classic()

p5 = nl.temp.bay.summer %>%
  ggplot(aes(x = datetime_ct, y = temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Noelle Lucey - Tranquilo Bay - Summer") +
  theme_classic()

p6 = mj.temp.cayo %>%
  ggplot(aes(x = datetime_ct, y = Temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Maggie Johnson - Cayo Roldan") +
  ylim(26,34) +
  theme_classic()

p7 = mj.temp.hosp %>%
  ggplot(aes(x = datetime_ct, y = Temp))+
  geom_line(color = 'coral', lwd=1.5)+
  ggtitle("Maggie Johnson - Hospital Point") +
  ylim(26,34) +
  theme_classic()


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4)

# use Dan Barshis' function to calculate dtv for these loggers
nl.temp.vista.daily<-data.frame("DayRange"=tapply(nl.temp.vista$temp, nl.temp.vista$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nl.temp.vista$temp, nl.temp.vista$Day, min),"DayMax"=tapply(nl.temp.vista$temp, nl.temp.vista$Day, max), "DayMean"=tapply(nl.temp.vista$temp, nl.temp.vista$Day, mean))
nl.temp.vista.daily <- tibble::rownames_to_column(nl.temp.vista.daily, "date")
nl.temp.vista.daily$site <- "Bella_Vista"

nl.temp.wilson.daily<-data.frame("DayRange"=tapply(nl.temp.wilson$temp, nl.temp.wilson$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nl.temp.wilson$temp, nl.temp.wilson$Day, min),"DayMax"=tapply(nl.temp.wilson$temp, nl.temp.wilson$Day, max), "DayMean"=tapply(nl.temp.wilson$temp, nl.temp.wilson$Day, mean))
nl.temp.wilson.daily <- tibble::rownames_to_column(nl.temp.wilson.daily, "date")
nl.temp.wilson.daily$site <- "Cayo_Wilson"

nl.temp.caracol.daily<-data.frame("DayRange"=tapply(nl.temp.caracol$temp, nl.temp.caracol$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nl.temp.caracol$temp, nl.temp.caracol$Day, min),"DayMax"=tapply(nl.temp.caracol$temp, nl.temp.caracol$Day, max), "DayMean"=tapply(nl.temp.caracol$temp, nl.temp.caracol$Day, mean))
nl.temp.caracol.daily <- tibble::rownames_to_column(nl.temp.caracol.daily, "date")
nl.temp.caracol.daily$site <- "Punta_Caracol"

nl.temp.bay.spring.daily<-data.frame("DayRange"=tapply(nl.temp.bay.spring$temp, nl.temp.bay.spring$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nl.temp.bay.spring$temp, nl.temp.bay.spring$Day, min),"DayMax"=tapply(nl.temp.bay.spring$temp, nl.temp.bay.spring$Day, max), "DayMean"=tapply(nl.temp.bay.spring$temp, nl.temp.bay.spring$Day, mean))
nl.temp.bay.spring.daily <- tibble::rownames_to_column(nl.temp.bay.spring.daily, "date")
nl.temp.bay.spring.daily$site <- "Tranquilo_Bay_Spring"

nl.temp.bay.summer.daily<-data.frame("DayRange"=tapply(nl.temp.bay.summer$temp, nl.temp.bay.summer$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nl.temp.bay.summer$temp, nl.temp.bay.summer$Day, min),"DayMax"=tapply(nl.temp.bay.summer$temp, nl.temp.bay.summer$Day, max), "DayMean"=tapply(nl.temp.bay.summer$temp, nl.temp.bay.summer$Day, mean))
nl.temp.bay.summer.daily <- tibble::rownames_to_column(nl.temp.bay.summer.daily, "date")
nl.temp.bay.summer.daily$site <- "Tranquilo_Bay_Summer"

crawl.cay.daily<-data.frame("DayRange"=tapply(crawl_cay$Temp, crawl_cay$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(crawl_cay$Temp, crawl_cay$Day, min),"DayMax"=tapply(crawl_cay$Temp, crawl_cay$Day, max), "DayMean"=tapply(crawl_cay$Temp, crawl_cay$Day, mean))
crawl.cay.daily <- tibble::rownames_to_column(crawl.cay.daily, "date")
crawl.cay.daily$site <- "Crawl_Cay"
#write.csv(crawl.cay.daily, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/CrawlCay_Daily.csv", row.names = FALSE)


# combine all of Maggie and Noelle's data together:
all.loggers <- rbind(mj.temp.cayo.daily, mj.temp.hosp.daily, nl.temp.vista.daily, nl.temp.wilson.daily, nl.temp.caracol.daily, nl.temp.bay.spring.daily, nl.temp.bay.summer.daily)
str(all.loggers)
all.loggers$site <- factor(all.loggers$site)
head(all.loggers)


#write.csv(all.loggers, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllLoggers.csv", row.names = FALSE)

# plot boxplot of daily ranges
ggplot(all.loggers, aes(x=site, y=DayRange, fill=site)) +
  geom_boxplot()+
  scale_fill_brewer(palette="RdBu") +
  geom_jitter(shape=16, position=position_jitter(0.2), color="gray", alpha=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6))

summarySE(all.loggers, measurevar="DayRange", groupvars = "site")
#                   site   N  DayRange        sd         se         ci
# 1          Bella_Vista  60 1.0827500 0.6900857 0.08908968 0.17826804
# 2          Cayo_Roldan 258 0.6356434 0.2612238 0.01626309 0.03202588
# 3          Cayo_Wilson  65 0.9188462 0.3278141 0.04066034 0.08122836
# 4       Hospital_Point 300 0.9822233 0.3880608 0.02240470 0.04409088
# 5        Punta_Caracol  63 0.2689683 0.1142346 0.01439221 0.02876962
# 6 Tranquilo_Bay_Spring  38 1.4920526 0.6484901 0.10519898 0.21315339
# 7 Tranquilo_Bay_Summer  60 5.6892000 2.9233902 0.37740806 0.75519177


summarySE(all.loggers, measurevar="DayMean", groupvars = "site")
#                   site   N  DayMean        sd         se         ci
# 1          Bella_Vista  60 30.49130 0.8029545 0.10366098 0.20742515
# 2          Cayo_Roldan 258 29.42425 0.5992965 0.03731057 0.07347338
# 3          Cayo_Wilson  65 29.19707 0.6489998 0.08049852 0.16081428
# 4       Hospital_Point 300 29.21883 1.0832968 0.06254417 0.12308252
# 5        Punta_Caracol  63 26.94833 0.3031607 0.03819466 0.07635003
# 6 Tranquilo_Bay_Spring  38 29.28549 0.6696942 0.10863875 0.22012302
# 7 Tranquilo_Bay_Summer  60 27.40121 1.1669257 0.15064946 0.30144888




#### SHINY ####
# no mods required here - skip
ui <- fluidPage(
  radioButtons("plotType", "Plot Type:", choices = c("ggplotly", "plotly")),
  plotlyOutput("plot"),
  verbatimTextOutput("hover"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brush"),
  verbatimTextOutput("zoom"))

# THIS SECTION TO BE MODIFIED
server <- function(input, output, session) {
  output$plot <- renderPlotly({
    # USE THE KEY AESTHETIC ARGUJMENT TO HELP UNIQUE IDENTIFY SELECTED OBSERVATIONS
    key <- master.melt$logger
    if (identical(input$plotType, "ggplotly")) {
      # YOUR PLOT CODE HERE (STORE AS p)
      p <- master.melt %>%
        ggplot(aes(x = time, y = temp, group = logger,color = logger))+
        geom_smooth()+
        scale_color_manual(values = palsite)+
        xlab('Time of Day')+
        ylab(expression(paste("Temp (",degree,"C)")))+
        scale_x_datetime(
          breaks = seq(as.POSIXct("2019-03-04 00:00:00"),as.POSIXct("2019-03-05 00:00:00 "), "4 hours"), labels = c('00:00', '04:00','08:00','12:00','16:00','20:00','00:00'))+
        facet_wrap(~logger, 1,5)
      ggplotly(p) %>% layout(dragmode = "select")
    } else {
      # PLOTY VERSION OPTION
      plot_ly(melt.mmmm, x = ~datetime, y = ~value, key = ~key) %>%
        layout(dragmode = "select")
    }
  })
  # NO FURTHER MODS REQUIRED
  output$hover <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  })
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click events appear here (double-click to clear)" else d
  })

  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })

  output$zoom <- renderPrint({
    d <- event_data("plotly_relayout")
    if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
  })

}

# OPEN THE APP
shinyApp(ui, server)

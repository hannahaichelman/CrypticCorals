########## 
  ##  TRIMMING HOBO DATA BY DATE ##
#   project             PANAMA TVE
#   author              Brooke Benson
#   updated             Thu Feb 21 2019

#   function:
#   Trims HOBO data from .csv files to reflect logger deployment and recovery dates and export as txt for use in xts/zoo

##########
# required packages #

library(tidyverse)
library(dplyr)

##########
# importing, trimming, and exporting files #

# set to wd for csv files
setwd("/Users/bebenson/Desktop/Panama Temperature Data/Raw HOBO Outputs")


# read in the csv
# change the datetime format
# trim using tidyverse to reflect logger deployment an recovery dates
# export as txt file
Cayo.OR3.arr3 <-read.csv("CayodeAgua-OR3_array3.csv")
Cayo.OR3.arr3$datetimew <- as.POSIXct(strptime(Cayo.OR3.arr3$datetime, format = "%m/%d/%y %H:%M"))
Cayo.OR3.arr3.trmd <- Cayo.OR3.arr3 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(Cayo.OR3.arr3.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/Cayo.OR3.arr3.trmd.txt', sep = '\t')

Cayo.OR3.arr10 <-read.csv("CayodeAgua-OR3_array10.csv")
Cayo.OR3.arr10$datetimew <- as.POSIXct(strptime(Cayo.OR3.arr10$datetime, format = "%m/%d/%y %H:%M"))
Cayo.OR3.arr10.trmd <- Cayo.OR3.arr10 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(Cayo.OR3.arr10.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/Cayo.OR3.arr10.trmd.txt', sep = '\t')

Drago.OR4 <- read.csv("DragoMar-OR4.csv")
Drago.OR4$datetimew <- as.POSIXct(strptime(Drago.OR4$datetime, format = "%m/%d/%y %H:%M"))
Drago.OR4.trmd <- Drago.OR4 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(Drago.OR4.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/Drago.OR4.trmd.txt', sep = '\t')

Punta.IR1 <- read.csv("PuntaDonato-IR1.csv")
Punta.IR1$datetimew <- as.POSIXct(strptime(Punta.IR1$datetime, format = "%m/%d/%y %H:%M"))
Punta.IR1.trmd <- Punta.IR1 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(Punta.IR1.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/Punta.IR1.trmd.txt', sep = '\t')

STRI.IR2 <- read.csv("STRIPoint-IR2.csv")
STRI.IR2$datetimew <- as.POSIXct(strptime(STRI.IR2$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.trmd <- STRI.IR2 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.trmd.txt', sep = '\t')

STRI.IR2.arr4 <- read.csv("STRIPoint-IR2_array4.csv")
STRI.IR2.arr4$datetimew <- as.POSIXct(strptime(STRI.IR2.arr4$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr4.trmd <- STRI.IR2.arr4 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr4.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr4.trmd.txt', sep = '\t')

STRI.IR2.arr5 <- read.csv("STRIPoint-IR2_array5.csv")
STRI.IR2.arr5$datetimew <- as.POSIXct(strptime(STRI.IR2.arr5$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr5.trmd <- STRI.IR2.arr5 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr5.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr5.trmd.txt', sep = '\t')

STRI.IR2.arr6 <- read.csv("STRIPoint-IR2_array6.csv")
STRI.IR2.arr6$datetimew <- as.POSIXct(strptime(STRI.IR2.arr6$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr6.trmd <- STRI.IR2.arr6 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr6.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr6.trmd.txt', sep = '\t')

STRI.IR2.arr7 <- read.csv("STRIPoint-IR2_array7.csv")
STRI.IR2.arr7$datetimew <- as.POSIXct(strptime(STRI.IR2.arr7$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr7.trmd <- STRI.IR2.arr7 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr7.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr7.trmd.txt', sep = '\t')

STRI.IR2.arr10 <- read.csv("STRIPoint-IR2_array10.csv")
STRI.IR2.arr10$datetimew <- as.POSIXct(strptime(STRI.IR2.arr10$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr10.trmd <- Cayo.OR3.arr3 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr10.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr10.trmd.txt', sep = '\t')

STRI.IR2.arr12 <- read.csv("STRIPoint-IR2_array12.csv")
STRI.IR2.arr12$datetimew <- as.POSIXct(strptime(STRI.IR2.arr12$datetime, format = "%m/%d/%y %H:%M"))
STRI.IR2.arr12.trmd <- STRI.IR2.arr12 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(STRI.IR2.arr12.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/STRI.IR2.arr12.trmd.txt', sep = '\t')

Cristo.IR3.arr1 <- read.csv("Cristobal-IR3_array1.csv")
Cristo.IR3.arr1$datetimew <- as.POSIXct(strptime(Cristo.IR3.arr1$datetime, format = "%m/%d/%y %H:%M"))
Cristo.IR3.arr1.trmd <- Cristo.IR3.arr1 %>%
  dplyr::filter(datetimew >= '2015-06-10 00:00:00' & datetimew <= '2016-08-14 11:59:59')
write.table(Cristo.IR3.arr1.trmd,'/Users/bebenson/Desktop/Panama Temperature Data/HOBO txt Files/Cristo.IR3.arr1.trmd.txt', sep = '\t')





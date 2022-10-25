##### METADATA #####
##  IN SITU TEMPERATURE DATA ##
#   project             PANAMA TVE
#   author              Brooke Benson
#   updated             Tue Mar 5 2019

#   function:
#   Visualize/analyze temperature data from HOBO loggers

#   TOC by line:
#       Shiny                                       19
#       Useful Notation Notes                       74
#       Required Packages                           85
#       Frequent Color Palettes                     103
#       Importing and Formatting Data               107
#       Monthly Max, Min, and Mean                  136
#       Compiling All Logged Data in Longform       219


#### SHINY ####
# no mods required here
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


#### Useful Notation Notes ####

#   xts
#       sample.xts['2007']            all of 2007
#       sample.xts['2007-03/']        March 2007 to the end of the data set
#       sample.xts['2007-03/2007']    March 2007 to the end of 2007
#       sample.xts['/']               the whole data set
#       sample.xts['/2007']           the beginning of the data through 2007
#       sample.xts['2007-01-03']      just the 3rd of January 2007


##### I. Required Packages #####
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


##### II: Frequent Color Palettes ####
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_lineage <- c("L1" = "#bcbddc", "L2" = "#756bb1")

##### III: Importing and Formatting Data #####

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
Cayo.OR3.arr3 <- xts(zoo(read.table("Cayo.OR3.arr3.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,7,1,0,0,0),ISOdate(2016,6,30,0,0,0), "15 min", tz="GMT")))
Drago.OR4 <- xts(zoo(read.table("Drago.OR4.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,7,1,0,0,0),ISOdate(2016,6,30,0,0,0), "15 min", tz="GMT")))
Punta.IR1 <- xts(zoo(read.table("Punta.IR1.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,7,1,0,0,0),ISOdate(2016,6,30,0,0,0), "15 min", tz="GMT")))
STRI.IR2 <- xts(zoo(read.table("STRI.IR2.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,7,1,0,0,0),ISOdate(2016,6,30,0,0,0), "15 min", tz="GMT")))
Cristo.IR3.arr1 <- xts(zoo(read.table("Cristo.IR3.arr1.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,7,1,0,0,0),ISOdate(2016,6,30,0,0,0), "15 min", tz="GMT")))
# fix cristobal temp to degrees C because the logger was deployed in F
Cristo.IR3.arr1 <- ((Cristo.IR3.arr1 - 32) * (5/9))

# STRI array loggers
STRI.IR2.arr4 <- xts(zoo(read.table("STRI.IR2.arr4.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr5 <- xts(zoo(read.table("STRI.IR2.arr5.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr6 <- xts(zoo(read.table("STRI.IR2.arr6.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr7 <- xts(zoo(read.table("STRI.IR2.arr7.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr10 <- xts(zoo(read.table("STRI.IR2.arr10.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))
STRI.IR2.arr12 <- xts(zoo(read.table("STRI.IR2.arr12.trmd.txt", header=TRUE, sep = '\t')$temp,seq.POSIXt(ISOdate(2015,6,10,0,0,0),ISOdate(2016,7,31,0,0,0), "15 min", tz="GMT")))


##### IV: Monthly Max, Min, and Mean #####

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

monthly.mmm$month <- ifelse(monthly.mmm$date == '2015-07-31', 'jul',
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
                                                                                                  ifelse(monthly.mmm$date == '2016-06-29', 'jun',
                                                                                                         'NA'))))))))))))

monthly.mmm$season <- ifelse(monthly.mmm$date  >= '2015-07-31' & monthly.mmm$date <= '2015-09-30', 'julaugsep',
                             ifelse(monthly.mmm$date >= '2015-10-31' & monthly.mmm$date <= '2015-12-31', 'octnovdec',
                                    ifelse(monthly.mmm$date >= '2016-01-31' & monthly.mmm$date <= '2016-03-31', 'janfebmar',
                                           ifelse(monthly.mmm$date >= '2016-04-30' & monthly.mmm$date <= '2016-06-29', 'aprmayjun', 'NA'))))


#write.csv(monthly.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/MonthlyTempParams.csv", row.names = FALSE)

aggregate(range ~ logger, data = monthly.mmm, FUN = "mean")

##### IV: Monthly Max, Min, and Mean Daytime Data Only #####


# (step 1) subset loggers to include daylight hours only
# https://www.worlddata.info/america/panama/sunset.php
# According to above website, seems like a reasonable average day is between 7am and 6pm (18:00)
Cayo.OR3.arr3.day <- Cayo.OR3.arr3["T07/T17"]
Drago.OR4.day <- Drago.OR4["T07/T17"]
Punta.IR1.day <- Punta.IR1["T07/T17"]
STRI.IR2.day <- STRI.IR2["T07/T17"]
Cristo.IR3.arr1.day <- Cristo.IR3.arr1["T07/T17"]

# i. maximum temperature by month #

# read in this csv file to skip running this section of code
monthly.mmm.day <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/MonthlyTempParams_daytime.csv")
###

# (step 1) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.day.mmax <- apply.monthly(Cayo.OR3.arr3.day['/'], FUN = max)
Drago.OR4.day.mmax <- apply.monthly(Drago.OR4.day['/'], FUN = max)
Punta.IR1.day.mmax <- apply.monthly(Punta.IR1.day['/'], FUN = max)
STRI.IR2.day.mmax <- apply.monthly(STRI.IR2.day['/'], FUN = max)
Cristo.IR3.arr1.day.mmax <- apply.monthly(Cristo.IR3.arr1.day['/'], FUN = max)

# (step 2) merge the sites together
#          can only merge two objects at a time
merge.day <- merge(Cayo.OR3.arr3.day.mmax,Drago.OR4.day.mmax)
merge.day <- merge(merge.day, Punta.IR1.day.mmax)
merge.day <- merge(merge.day,STRI.IR2.day.mmax)
merged.day <- merge(merge.day, Cristo.IR3.arr1.day.mmax)

# (step 3) rename the columns to keep identifying info
colnames(merged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 4) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
merged.day <- as.data.frame(merged.day)
merged.day <- setDT(merged.day, keep.rownames = T)
colnames(merged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) melt to put the data in longform
melt.day.mmax <- melt(merged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.day.mmax) <- c('datetime', 'logger', 'mmax')
head(melt.day.mmax)

# repeat steps 1-5 for all parameters of interest...

# ii. minimum temperature by month #
Cayo.OR3.arr3.day.mmin<- apply.monthly(Cayo.OR3.arr3.day['/'], FUN = min)
Drago.OR4.day.mmin <- apply.monthly(Drago.OR4.day['/'], FUN = min)
Punta.IR1.day.mmin <- apply.monthly(Punta.IR1.day['/'], FUN = min)
STRI.IR2.day.mmin <- apply.monthly(STRI.IR2.day['/'], FUN = min)
Cristo.IR3.arr1.day.mmin <- apply.monthly(Cristo.IR3.arr1.day['/'], FUN = min)
merge.day <- merge(Cayo.OR3.arr3.day.mmin,Drago.OR4.day.mmin)
merge.day <- merge(merge.day, Punta.IR1.day.mmin)
merge.day <- merge(merge.day,STRI.IR2.day.mmin)
merged.day <- merge(merge.day, Cristo.IR3.arr1.day.mmin)
colnames(merged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.day <- as.data.frame(merged.day)
merged.day <- setDT(merged.day, keep.rownames = T)
colnames(merged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.day.mmin <- melt(merged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.day.mmin) <- c('datetime', 'logger', 'mmin')
head(melt.day.mmin)

# iii. mean temperature by month #
Cayo.OR3.arr3.day.mmean<- apply.monthly(Cayo.OR3.arr3.day['/'], FUN = mean)
Drago.OR4.day.mmean <- apply.monthly(Drago.OR4.day['/'], FUN = mean)
Punta.IR1.day.mmean <- apply.monthly(Punta.IR1.day['/'], FUN = mean)
STRI.IR2.day.mmean <- apply.monthly(STRI.IR2.day['/'], FUN = mean)
Cristo.IR3.arr1.day.mmean <- apply.monthly(Cristo.IR3.arr1.day['/'], FUN = mean)
merge.day <- merge(Cayo.OR3.arr3.day.mmean,Drago.OR4.day.mmean)
merge.day <- merge(merge.day, Punta.IR1.day.mmean)
merge.day <- merge(merge.day,STRI.IR2.day.mmean)
merged.day <- merge(merge.day, Cristo.IR3.arr1.day.mmean)
colnames(merged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.day <- as.data.frame(merged.day)
merged.day <- setDT(merged.day, keep.rownames = T)
colnames(merged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.day.mmean <- melt(merged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.day.mmean) <- c('datetime', 'logger', 'mmean')
head(melt.day.mmean)

# merge the max, min, and mean and melt
m.maxmin.day <- merge(melt.day.mmax, melt.day.mmin)
m.maxminmean.day <- merge(m.maxmin.day, melt.day.mmean)
melt.mmmm.day <- melt(m.maxminmean.day, id.vars = c('datetime','logger'), measure.vars = c('mmax','mmin','mmean'))
head(melt.mmmm.day)

# it's easier to do use setDT to get "typical day" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of day
monthly.mmm.day<-setDT(melt.mmmm.day)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
monthly.mmm.day$range <- (monthly.mmm.day$max - monthly.mmm.day$min)
monthly.mmm.day$logger <- factor(monthly.mmm.day$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(monthly.mmm.day)

#add in separate date column, since time doesn't matter for weekly data
monthly.mmm.day$date <- as.Date(monthly.mmm.day$datetime,"%Y-%m-%d" )

# add in separate month column
head(monthly.mmm.day)

monthly.mmm.day$month <- ifelse(monthly.mmm.day$date == '2015-07-31', 'jul',
                                ifelse(monthly.mmm.day$date == '2015-08-31', 'aug',
                                       ifelse(monthly.mmm.day$date == '2015-09-30', 'sep',
                                              ifelse(monthly.mmm.day$date == '2015-10-31', 'oct',
                                                     ifelse(monthly.mmm.day$date == '2015-11-30', 'nov',
                                                            ifelse(monthly.mmm.day$date == '2015-12-31', 'dec',
                                                                   ifelse(monthly.mmm.day$date == '2016-01-31', 'jan',
                                                                          ifelse(monthly.mmm.day$date == '2016-02-29', 'feb',
                                                                                 ifelse(monthly.mmm.day$date == '2016-03-31', 'mar',
                                                                                        ifelse(monthly.mmm.day$date == '2016-04-30', 'apr',
                                                                                               ifelse(monthly.mmm.day$date == '2016-05-31', 'may',
                                                                                                      ifelse(monthly.mmm.day$date == '2016-06-29', 'jun',
                                                                                                             'NA'))))))))))))

monthly.mmm.day$season <- ifelse(monthly.mmm.day$date  >= '2015-07-31' & monthly.mmm.day$date <= '2015-09-30', 'julaugsep',
                                 ifelse(monthly.mmm.day$date >= '2015-10-31' & monthly.mmm.day$date <= '2015-12-31', 'octnovdec',
                                        ifelse(monthly.mmm.day$date >= '2016-01-31' & monthly.mmm.day$date <= '2016-03-31', 'janfebmar',
                                               ifelse(monthly.mmm.day$date >= '2016-04-30' & monthly.mmm.day$date <= '2016-06-29', 'aprmayjun', 'NA'))))

#write.csv(monthly.mmm.day, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/MonthlyTempParams_daytime.csv", row.names = FALSE)

aggregate(range ~ logger, data = monthly.mmm.day, FUN = "mean")

##### IV: Monthly Max, Min, and Mean Nighttime Data Only #####

# (step 1) subset loggers to include nighttime hours only
# https://www.worlddata.info/america/panama/sunset.php
# According to above website, seems like a reasonable average day is between 7am and 6pm (18:00)
# Therefore, we'll call average night between 6pm (18:00) and 7am (07:00)
Cayo.OR3.arr3.night <- Cayo.OR3.arr3["T18/T06"]
Drago.OR4.night <- Drago.OR4["T18/T06"]
Punta.IR1.night <- Punta.IR1["T18/T06"]
STRI.IR2.night <- STRI.IR2["T18/T06"]
Cristo.IR3.arr1.night <- Cristo.IR3.arr1["T18/T06"]

# i. maximum temperature by month #

# read in this csv to skip running this section of code
monthly.mmm.night <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/MonthlyTempParams_nighttime.csv")
###

# (step 1) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.night.mmax <- apply.monthly(Cayo.OR3.arr3.night['/'], FUN = max)
Drago.OR4.night.mmax <- apply.monthly(Drago.OR4.night['/'], FUN = max)
Punta.IR1.night.mmax <- apply.monthly(Punta.IR1.night['/'], FUN = max)
STRI.IR2.night.mmax <- apply.monthly(STRI.IR2.night['/'], FUN = max)
Cristo.IR3.arr1.night.mmax <- apply.monthly(Cristo.IR3.arr1.night['/'], FUN = max)

# (step 2) merge the sites together
#          can only merge two objects at a time
merge.night <- merge(Cayo.OR3.arr3.night.mmax,Drago.OR4.night.mmax)
merge.night <- merge(merge.night, Punta.IR1.night.mmax)
merge.night <- merge(merge.night,STRI.IR2.night.mmax)
merged.night <- merge(merge.night, Cristo.IR3.arr1.night.mmax)

# (step 3) rename the columns to keep identifying info
colnames(merged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 4) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
merged.night <- as.data.frame(merged.night)
merged.night <- setDT(merged.night, keep.rownames = T)
colnames(merged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) melt to put the data in longform
melt.night.mmax <- melt(merged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.night.mmax) <- c('datetime', 'logger', 'mmax')
head(melt.night.mmax)

# repeat steps 1-5 for all parameters of interest...

# ii. minimum temperature by month #
Cayo.OR3.arr3.night.mmin<- apply.monthly(Cayo.OR3.arr3.night['/'], FUN = min)
Drago.OR4.night.mmin <- apply.monthly(Drago.OR4.night['/'], FUN = min)
Punta.IR1.night.mmin <- apply.monthly(Punta.IR1.night['/'], FUN = min)
STRI.IR2.night.mmin <- apply.monthly(STRI.IR2.night['/'], FUN = min)
Cristo.IR3.arr1.night.mmin <- apply.monthly(Cristo.IR3.arr1.night['/'], FUN = min)
merge.night <- merge(Cayo.OR3.arr3.night.mmin,Drago.OR4.night.mmin)
merge.night <- merge(merge.night, Punta.IR1.night.mmin)
merge.night <- merge(merge.night,STRI.IR2.night.mmin)
merged.night <- merge(merge.night, Cristo.IR3.arr1.night.mmin)
colnames(merged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.night <- as.data.frame(merged.night)
merged.night <- setDT(merged.night, keep.rownames = T)
colnames(merged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.night.mmin <- melt(merged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.night.mmin) <- c('datetime', 'logger', 'mmin')
head(melt.night.mmin)

# iii. mean temperature by month #
Cayo.OR3.arr3.night.mmean<- apply.monthly(Cayo.OR3.arr3.night['/'], FUN = mean)
Drago.OR4.night.mmean <- apply.monthly(Drago.OR4.night['/'], FUN = mean)
Punta.IR1.night.mmean <- apply.monthly(Punta.IR1.night['/'], FUN = mean)
STRI.IR2.night.mmean <- apply.monthly(STRI.IR2.night['/'], FUN = mean)
Cristo.IR3.arr1.night.mmean <- apply.monthly(Cristo.IR3.arr1.night['/'], FUN = mean)
merge.night <- merge(Cayo.OR3.arr3.night.mmean,Drago.OR4.night.mmean)
merge.night <- merge(merge.night, Punta.IR1.night.mmean)
merge.night <- merge(merge.night,STRI.IR2.night.mmean)
merged.night <- merge(merge.night, Cristo.IR3.arr1.night.mmean)
colnames(merged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.night <- as.data.frame(merged.night)
merged.night <- setDT(merged.night, keep.rownames = T)
colnames(merged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
melt.night.mmean <- melt(merged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.night.mmean) <- c('datetime', 'logger', 'mmean')
head(melt.night.mmean)

# merge the max, min, and mean and melt
m.maxmin.night <- merge(melt.night.mmax, melt.night.mmin)
m.maxminmean.night <- merge(m.maxmin.night, melt.night.mmean)
melt.mmmm.night <- melt(m.maxminmean.night, id.vars = c('datetime','logger'), measure.vars = c('mmax','mmin','mmean'))
head(melt.mmmm.night)

# it's easier to do use setDT to get "typical night" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of night
monthly.mmm.night<-setDT(melt.mmmm.night)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
monthly.mmm.night$range <- (monthly.mmm.night$max - monthly.mmm.night$min)
monthly.mmm.night$logger <- factor(monthly.mmm.night$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(monthly.mmm.night)

#add in separate date column, since time doesn't matter for weekly data
monthly.mmm.night$date <- as.Date(monthly.mmm.night$datetime,"%Y-%m-%d" )

# add in separate month column
head(monthly.mmm.night)

monthly.mmm.night$month <- ifelse(monthly.mmm.night$date == '2015-07-31', 'jul',
                                  ifelse(monthly.mmm.night$date == '2015-08-31', 'aug',
                                         ifelse(monthly.mmm.night$date == '2015-09-30', 'sep',
                                                ifelse(monthly.mmm.night$date == '2015-10-31', 'oct',
                                                       ifelse(monthly.mmm.night$date == '2015-11-30', 'nov',
                                                              ifelse(monthly.mmm.night$date == '2015-12-31', 'dec',
                                                                     ifelse(monthly.mmm.night$date == '2016-01-31', 'jan',
                                                                            ifelse(monthly.mmm.night$date == '2016-02-29', 'feb',
                                                                                   ifelse(monthly.mmm.night$date == '2016-03-31', 'mar',
                                                                                          ifelse(monthly.mmm.night$date == '2016-04-30', 'apr',
                                                                                                 ifelse(monthly.mmm.night$date == '2016-05-31', 'may',
                                                                                                        ifelse(monthly.mmm.night$date == '2016-06-29', 'jun',
                                                                                                               'NA'))))))))))))

monthly.mmm.night$season <- ifelse(monthly.mmm.night$date  >= '2015-07-31' & monthly.mmm.night$date <= '2015-09-30', 'julaugsep',
                                   ifelse(monthly.mmm.night$date >= '2015-10-31' & monthly.mmm.night$date <= '2015-12-31', 'octnovdec',
                                          ifelse(monthly.mmm.night$date >= '2016-01-31' & monthly.mmm.night$date <= '2016-03-31', 'janfebmar',
                                                 ifelse(monthly.mmm.night$date >= '2016-04-30' & monthly.mmm.night$date <= '2016-06-29', 'aprmayjun', 'NA'))))

#write.csv(monthly.mmm.night, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/MonthlyTempParams_nighttime.csv", row.names = FALSE)

aggregate(range ~ logger, data = monthly.mmm.night, FUN = "mean")

##### V: Weekly Max, Min, Mean, and Range #####

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

weekly.mmm$month <- ifelse(weekly.mmm$date  >= '2015-07-05' & weekly.mmm$date <= '2015-07-26', 'jul',
                           ifelse(weekly.mmm$date >= '2015-08-02' & weekly.mmm$date <= '2015-08-30', 'aug',
                                  ifelse(weekly.mmm$date >= '2015-09-06' & weekly.mmm$date <= '2015-09-27', 'sep',
                                         ifelse(weekly.mmm$date >= '2015-10-04' & weekly.mmm$date <= '2015-10-25', 'oct',
                                                ifelse(weekly.mmm$date >= '2015-11-01' & weekly.mmm$date <= '2015-11-29', 'nov',
                                                       ifelse(weekly.mmm$date >= '2015-12-06' & weekly.mmm$date <= '2015-12-27', 'dec',
                                                              ifelse(weekly.mmm$date >= '2016-01-03' & weekly.mmm$date <= '2016-01-31', 'jan',
                                                                     ifelse(weekly.mmm$date >= '2016-02-07' & weekly.mmm$date <= '2016-02-28', 'feb',
                                                                            ifelse(weekly.mmm$date >= '2016-03-06' & weekly.mmm$date <= '2016-03-27', 'mar',
                                                                                   ifelse(weekly.mmm$date >= '2016-04-03' & weekly.mmm$date <= '2016-04-24', 'apr',
                                                                                          ifelse(weekly.mmm$date >= '2016-05-01' & weekly.mmm$date <= '2016-05-29', 'may',
                                                                                                 ifelse(weekly.mmm$date >= '2016-06-05' & weekly.mmm$date <= '2016-06-29', 'jun',
                                                                                                        'NA'))))))))))))

weekly.mmm$season <- ifelse(weekly.mmm$date  >= '2015-07-05' & weekly.mmm$date <= '2015-09-27', 'julaugsep',
                            ifelse(weekly.mmm$date >= '2015-10-04' & weekly.mmm$date <= '2015-12-27', 'octnovdec',
                                   ifelse(weekly.mmm$date >= '2016-01-03' & weekly.mmm$date <= '2016-03-27', 'janfebmar',
                                          ifelse(weekly.mmm$date >= '2016-04-03' & weekly.mmm$date <= '2016-06-29', 'aprmayjun', 'NA'))))

#write.csv(weekly.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams.csv", row.names = FALSE)

aggregate(range ~ logger, data = weekly.mmm, FUN = "mean")

##### VI: Weekly Max, Min, Mean, and Range Daylight Hours Only #####

# i. maximum temperature by week, daylight hours #

# read in this csv to skip running this section of code:
weekly.mmm.day <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/WeeklyTempParams_daytime.csv")
###

# (step 2) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.day.wmax <- apply.weekly(Cayo.OR3.arr3.day['/'], FUN = max)
Drago.OR4.day.wmax <- apply.weekly(Drago.OR4.day['/'], FUN = max)
Punta.IR1.day.wmax <- apply.weekly(Punta.IR1.day['/'], FUN = max)
Punta.IR1.day.wmax <- apply.weekly(STRI.IR2.day['/'], FUN = max)
Cristo.IR3.arr1.day.wmax <- apply.weekly(Cristo.IR3.arr1.day['/'], FUN = max)

# (step 3) merge the sites together
#          can only merge two objects at a time
weekmerge.day <- merge(Cayo.OR3.arr3.day.wmax,Drago.OR4.day.wmax)
weekmerge.day <- merge(weekmerge.day, Punta.IR1.day.wmax)
weekmerge.day <- merge(weekmerge.day,Punta.IR1.day.wmax)
weekmerged.day <- merge(weekmerge.day, Cristo.IR3.arr1.day.wmax)

# (step 4) rename the columns to keep identifying info
colnames(weekmerged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
weekmerged.day <- as.data.frame(weekmerged.day)
weekmerged.day <- setDT(weekmerged.day, keep.rownames = T)
colnames(weekmerged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 6) melt to put the data in longform
melt.wmax.day <- melt(weekmerged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmax.day) <- c('datetime', 'logger', 'wmax.day')
head(melt.wmax.day)

# repeat steps 1-6 for all parameters of interest...

# ii. minimum temperature by week #
Cayo.OR3.arr3.day.wmin<- apply.weekly(Cayo.OR3.arr3.day['/'], FUN = min)
Drago.OR4.day.wmin <- apply.weekly(Drago.OR4.day['/'], FUN = min)
Punta.IR1.day.wmin <- apply.weekly(Punta.IR1.day['/'], FUN = min)
STRI.IR2.day.wmin <- apply.weekly(STRI.IR2.day['/'], FUN = min)
Cristo.IR3.arr1.day.wmin <- apply.weekly(Cristo.IR3.arr1.day['/'], FUN = min)

weekmerge.day <- merge(Cayo.OR3.arr3.day.wmin,Drago.OR4.day.wmin)
weekmerge.day <- merge(weekmerge.day, Punta.IR1.day.wmin)
weekmerge.day <- merge(weekmerge.day,STRI.IR2.day.wmin)
weekmerged.day <- merge(weekmerge.day, Cristo.IR3.arr1.day.wmin)
colnames(weekmerged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

weekmerged.day <- as.data.frame(weekmerged.day)
weekmerged.day <- setDT(weekmerged.day, keep.rownames = T)
colnames(weekmerged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

melt.wmin.day <- melt(weekmerged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmin.day) <- c('datetime', 'logger', 'wmin.day')
head(melt.wmin.day)

# iii. mean temperature by week #
Cayo.OR3.arr3.day.wmean<- apply.weekly(Cayo.OR3.arr3.day['/'], FUN = mean)
Drago.OR4.day.wmean <- apply.weekly(Drago.OR4.day['/'], FUN = mean)
Punta.IR1.day.wmean <- apply.weekly(Punta.IR1.day['/'], FUN = mean)
STRI.IR2.day.wmean <- apply.weekly(STRI.IR2.day['/'], FUN = mean)
Cristo.IR3.arr1.day.wmean <- apply.weekly(Cristo.IR3.arr1.day['/'], FUN = mean)

weekmerge.day <- merge(Cayo.OR3.arr3.day.wmean,Drago.OR4.day.wmean)
weekmerge.day <- merge(weekmerge.day, Punta.IR1.day.wmean)
weekmerge.day <- merge(weekmerge.day,STRI.IR2.day.wmean)
weekmerged.day <- merge(weekmerge.day, Cristo.IR3.arr1.day.wmean)
colnames(weekmerged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

weekmerged.day <- as.data.frame(weekmerged.day)
weekmerged.day <- setDT(weekmerged.day, keep.rownames = T)
colnames(weekmerged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

melt.wmean.day <- melt(weekmerged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmean.day) <- c('datetime', 'logger', 'wmean.day')
head(melt.wmean.day)

# merge the max, min, and mean and melt
w.maxmin.day <- merge(melt.wmax.day, melt.wmin.day)
w.maxminmean.day <- merge(w.maxmin.day, melt.wmean.day)
melt.wmmm.day <- melt(w.maxminmean.day, id.vars = c('datetime','logger'), measure.vars = c('wmax.day','wmin.day','wmean.day'))
head(melt.wmmm.day)

# it's easier to do use setDT to get "typical day" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of day
weekly.mmm.day<-setDT(melt.wmmm.day)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
weekly.mmm.day$range <- (weekly.mmm.day$max - weekly.mmm.day$min)
weekly.mmm.day$logger <- factor(weekly.mmm.day$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(weekly.mmm.day)

#add in separate date column, since time doesn't matter for weekly data
weekly.mmm.day$date <- as.Date(weekly.mmm.day$datetime,"%Y-%m-%d" )

# add in month indication
head(weekly.mmm.day)

weekly.mmm.day$month <- ifelse(weekly.mmm.day$date  >= '2015-07-05' & weekly.mmm.day$date <= '2015-07-26', 'jul',
                               ifelse(weekly.mmm.day$date >= '2015-08-02' & weekly.mmm.day$date <= '2015-08-30', 'aug',
                                      ifelse(weekly.mmm.day$date >= '2015-09-06' & weekly.mmm.day$date <= '2015-09-27', 'sep',
                                             ifelse(weekly.mmm.day$date >= '2015-10-04' & weekly.mmm.day$date <= '2015-10-25', 'oct',
                                                    ifelse(weekly.mmm.day$date >= '2015-11-01' & weekly.mmm.day$date <= '2015-11-29', 'nov',
                                                           ifelse(weekly.mmm.day$date >= '2015-12-06' & weekly.mmm.day$date <= '2015-12-27', 'dec',
                                                                  ifelse(weekly.mmm.day$date >= '2016-01-03' & weekly.mmm.day$date <= '2016-01-31', 'jan',
                                                                         ifelse(weekly.mmm.day$date >= '2016-02-07' & weekly.mmm.day$date <= '2016-02-28', 'feb',
                                                                                ifelse(weekly.mmm.day$date >= '2016-03-06' & weekly.mmm.day$date <= '2016-03-27', 'mar',
                                                                                       ifelse(weekly.mmm.day$date >= '2016-04-03' & weekly.mmm.day$date <= '2016-04-24', 'apr',
                                                                                              ifelse(weekly.mmm.day$date >= '2016-05-01' & weekly.mmm.day$date <= '2016-05-29', 'may',
                                                                                                     ifelse(weekly.mmm.day$date >= '2016-06-05' & weekly.mmm.day$date <= '2016-06-29', 'jun',
                                                                                                            'NA'))))))))))))

# add in season info
weekly.mmm.day$season <- ifelse(weekly.mmm.day$date  >= '2015-07-05' & weekly.mmm.day$date <= '2015-09-27', 'julaugsep',
                                ifelse(weekly.mmm.day$date >= '2015-10-04' & weekly.mmm.day$date <= '2015-12-27', 'octnovdec',
                                       ifelse(weekly.mmm.day$date >= '2016-01-03' & weekly.mmm.day$date <= '2016-03-27', 'janfebmar',
                                              ifelse(weekly.mmm.day$date >= '2016-04-03' & weekly.mmm.day$date <= '2016-06-29', 'aprmayjun', 'NA'))))

#write.csv(weekly.mmm.day, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams_daytime.csv", row.names = FALSE)

head(weekly.mmm.day)

aggregate(range ~ logger, data = weekly.mmm.day, FUN = "mean")

##### VI: Weekly Max, Min, Mean, and Range Nighttime Hours Only #####

# i. maximum temperature by week, nighttime hours #

# read in this csv to skip running this section of code:
weekly.mmm.night <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/WeeklyTempParams_nighttime.csv")
###

# (step 2) use xts/zoo to get monthly max temperatures by site as xts/zoo objects

Cayo.OR3.arr3.night.wmax <- apply.weekly(Cayo.OR3.arr3.night['/'], FUN = max)
Drago.OR4.night.wmax <- apply.weekly(Drago.OR4.night['/'], FUN = max)
Punta.IR1.night.wmax <- apply.weekly(Punta.IR1.night['/'], FUN = max)
Punta.IR1.night.wmax <- apply.weekly(STRI.IR2.night['/'], FUN = max)
Cristo.IR3.arr1.night.wmax <- apply.weekly(Cristo.IR3.arr1.night['/'], FUN = max)

# (step 3) merge the sites together
#          can only merge two objects at a time
weekmerge.night <- merge(Cayo.OR3.arr3.night.wmax,Drago.OR4.night.wmax)
weekmerge.night <- merge(weekmerge.night, Punta.IR1.night.wmax)
weekmerge.night <- merge(weekmerge.night,Punta.IR1.night.wmax)
weekmerged.night <- merge(weekmerge.night, Cristo.IR3.arr1.night.wmax)

# (step 4) rename the columns to keep identifying info
colnames(weekmerged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 5) make it a data.frame and carry over the datetime info as a variable via keep row names, rename column
weekmerged.night <- as.data.frame(weekmerged.night)
weekmerged.night <- setDT(weekmerged.night, keep.rownames = T)
colnames(weekmerged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

# (step 6) melt to put the data in longform
melt.wmax.night <- melt(weekmerged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmax.night) <- c('datetime', 'logger', 'wmax.night')
head(melt.wmax.night)

# repeat steps 1-6 for all parameters of interest...

# ii. minimum temperature by week #
Cayo.OR3.arr3.night.wmin<- apply.weekly(Cayo.OR3.arr3.night['/'], FUN = min)
Drago.OR4.night.wmin <- apply.weekly(Drago.OR4.night['/'], FUN = min)
Punta.IR1.night.wmin <- apply.weekly(Punta.IR1.night['/'], FUN = min)
STRI.IR2.night.wmin <- apply.weekly(STRI.IR2.night['/'], FUN = min)
Cristo.IR3.arr1.night.wmin <- apply.weekly(Cristo.IR3.arr1.night['/'], FUN = min)

weekmerge.night <- merge(Cayo.OR3.arr3.night.wmin,Drago.OR4.night.wmin)
weekmerge.night <- merge(weekmerge.night, Punta.IR1.night.wmin)
weekmerge.night <- merge(weekmerge.night,STRI.IR2.night.wmin)
weekmerged.night <- merge(weekmerge.night, Cristo.IR3.arr1.night.wmin)
colnames(weekmerged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

weekmerged.night <- as.data.frame(weekmerged.night)
weekmerged.night <- setDT(weekmerged.night, keep.rownames = T)
colnames(weekmerged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

melt.wmin.night <- melt(weekmerged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmin.night) <- c('datetime', 'logger', 'wmin.night')
head(melt.wmin.night)

# iii. mean temperature by week #
Cayo.OR3.arr3.night.wmean<- apply.weekly(Cayo.OR3.arr3.night['/'], FUN = mean)
Drago.OR4.night.wmean <- apply.weekly(Drago.OR4.night['/'], FUN = mean)
Punta.IR1.night.wmean <- apply.weekly(Punta.IR1.night['/'], FUN = mean)
STRI.IR2.night.wmean <- apply.weekly(STRI.IR2.night['/'], FUN = mean)
Cristo.IR3.arr1.night.wmean <- apply.weekly(Cristo.IR3.arr1.night['/'], FUN = mean)

weekmerge.night <- merge(Cayo.OR3.arr3.night.wmean,Drago.OR4.night.wmean)
weekmerge.night <- merge(weekmerge.night, Punta.IR1.night.wmean)
weekmerge.night <- merge(weekmerge.night,STRI.IR2.night.wmean)
weekmerged.night <- merge(weekmerge.night, Cristo.IR3.arr1.night.wmean)
colnames(weekmerged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

weekmerged.night <- as.data.frame(weekmerged.night)
weekmerged.night <- setDT(weekmerged.night, keep.rownames = T)
colnames(weekmerged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')

melt.wmean.night <- melt(weekmerged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(melt.wmean.night) <- c('datetime', 'logger', 'wmean.night')
head(melt.wmean.night)

# merge the max, min, and mean and melt
w.maxmin.night <- merge(melt.wmax.night, melt.wmin.night)
w.maxminmean.night <- merge(w.maxmin.night, melt.wmean.night)
melt.wmmm.night <- melt(w.maxminmean.night, id.vars = c('datetime','logger'), measure.vars = c('wmax.night','wmin.night','wmean.night'))
head(melt.wmmm.night)

# it's easier to do use setDT to get "typical night" patterns than xtx/zoo because the unique datetime values in the xts/zoo preclude looking at things by time of night
weekly.mmm.night<-setDT(melt.wmmm.night)[, list(max = max(value), min = min(value), mean=mean(value)), by=list(logger,datetime)]
# add in range
weekly.mmm.night$range <- (weekly.mmm.night$max - weekly.mmm.night$min)
weekly.mmm.night$logger <- factor(weekly.mmm.night$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))
head(weekly.mmm.night)

#add in separate date column, since time doesn't matter for weekly data
weekly.mmm.night$date <- as.Date(weekly.mmm.night$datetime,"%Y-%m-%d" )

# add in month indication
head(weekly.mmm.night)

weekly.mmm.night$month <- ifelse(weekly.mmm.night$date  >= '2015-07-05' & weekly.mmm.night$date <= '2015-07-26', 'jul',
                                 ifelse(weekly.mmm.night$date >= '2015-08-02' & weekly.mmm.night$date <= '2015-08-30', 'aug',
                                        ifelse(weekly.mmm.night$date >= '2015-09-06' & weekly.mmm.night$date <= '2015-09-27', 'sep',
                                               ifelse(weekly.mmm.night$date >= '2015-10-04' & weekly.mmm.night$date <= '2015-10-25', 'oct',
                                                      ifelse(weekly.mmm.night$date >= '2015-11-01' & weekly.mmm.night$date <= '2015-11-29', 'nov',
                                                             ifelse(weekly.mmm.night$date >= '2015-12-06' & weekly.mmm.night$date <= '2015-12-27', 'dec',
                                                                    ifelse(weekly.mmm.night$date >= '2016-01-03' & weekly.mmm.night$date <= '2016-01-31', 'jan',
                                                                           ifelse(weekly.mmm.night$date >= '2016-02-07' & weekly.mmm.night$date <= '2016-02-28', 'feb',
                                                                                  ifelse(weekly.mmm.night$date >= '2016-03-06' & weekly.mmm.night$date <= '2016-03-27', 'mar',
                                                                                         ifelse(weekly.mmm.night$date >= '2016-04-03' & weekly.mmm.night$date <= '2016-04-24', 'apr',
                                                                                                ifelse(weekly.mmm.night$date >= '2016-05-01' & weekly.mmm.night$date <= '2016-05-29', 'may',
                                                                                                       ifelse(weekly.mmm.night$date >= '2016-06-05' & weekly.mmm.night$date <= '2016-06-29', 'jun',
                                                                                                              'NA'))))))))))))

# add in season info
weekly.mmm.night$season <- ifelse(weekly.mmm.night$date  >= '2015-07-05' & weekly.mmm.night$date <= '2015-09-27', 'julaugsep',
                                  ifelse(weekly.mmm.night$date >= '2015-10-04' & weekly.mmm.night$date <= '2015-12-27', 'octnovdec',
                                         ifelse(weekly.mmm.night$date >= '2016-01-03' & weekly.mmm.night$date <= '2016-03-27', 'janfebmar',
                                                ifelse(weekly.mmm.night$date >= '2016-04-03' & weekly.mmm.night$date <= '2016-06-29', 'aprmayjun', 'NA'))))

#write.csv(weekly.mmm.night, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams_nighttime.csv", row.names = FALSE)

head(weekly.mmm.night)

aggregate(range ~ logger, data = weekly.mmm.night, FUN = "mean")

#### Absolute Max/Min/Mean from Loggers and Satellite ####

# get the absolute max, min, mean etc for the deployment year
# deployment year lasts 2015-07-01 00:00:00 to 2016-06-30 00:00:00
pd <- data.frame(site = 'Punta.IR1', max = max(Punta.IR1), min = min(Punta.IR1), mean = mean(Punta.IR1), sd = sd(Punta.IR1), se = (sd(Punta.IR1)/sqrt(nrow(Punta.IR1))))
sp <- data.frame(site = 'STRI.IR2', max = max(STRI.IR2), min = min(STRI.IR2), mean = mean(STRI.IR2), sd = sd(STRI.IR2), se = (sd(STRI.IR2)/sqrt(nrow(STRI.IR2))))
ci <-data.frame(site = 'Cristo.IR3', max = max(Cristo.IR3.arr1), min = min(Cristo.IR3.arr1), mean = mean(Cristo.IR3.arr1), sd = sd(Cristo.IR3.arr1), se = (sd(Cristo.IR3.arr1)/sqrt(nrow(Cristo.IR3.arr1))))
ca <- data.frame(site = 'Cayo.OR3', max = max(Cayo.OR3.arr3), min = min(Cayo.OR3.arr3), mean = mean(Cayo.OR3.arr3), sd = sd(Cayo.OR3.arr3), se = (sd(Cayo.OR3.arr3)/sqrt(nrow(Cayo.OR3.arr3))))
dm <- data.frame(site = 'Drago.OR4', max = max(Drago.OR4), min = min(Drago.OR4), mean = mean(Drago.OR4), sd = sd(Drago.OR4), se = (sd(Drago.OR4)/sqrt(nrow(Drago.OR4))))

abspar <-rbind(pd, sp, ci, ca, dm)
abspar$ann.range <- (abspar$max - abspar$min)

# site      max    min     mean        sd          se ann.range
# Punta.IR1 32.15000 27.038 29.70215 0.8545586 0.004565135  5.112000
# STRI.IR2 31.58600 27.038 29.31217 0.7784966 0.004158804  4.548000
# Cristo.IR3 33.54778 27.235 29.88817 0.9600723 0.005128799  6.312778
# Cayo.OR3 32.43300 27.038 29.47132 0.9496586 0.005073168  5.395000
# Drago.OR4 31.66300 26.769 29.12715 0.8223259 0.004392944  4.894000


# now do this with daytime values only
# date range includes 2015-07-01 07:00:00	to 2016-06-29 17:45:00 (daytime only)
pd.day <- data.frame(site = 'Punta.IR1.day', max = max(Punta.IR1.day), min = min(Punta.IR1.day), mean = mean(Punta.IR1.day), sd = sd(Punta.IR1.day), se = (sd(Punta.IR1.day)/sqrt(nrow(Punta.IR1.day))))
sp.day <- data.frame(site = 'STRI.IR2.day', max = max(STRI.IR2.day), min = min(STRI.IR2.day), mean = mean(STRI.IR2.day), sd = sd(STRI.IR2.day), se = (sd(STRI.IR2.day)/sqrt(nrow(STRI.IR2.day))))
ci.day <-data.frame(site = 'Cristo.IR3.day', max = max(Cristo.IR3.arr1.day), min = min(Cristo.IR3.arr1.day), mean = mean(Cristo.IR3.arr1.day), sd = sd(Cristo.IR3.arr1.day), se = (sd(Cristo.IR3.arr1.day)/sqrt(nrow(Cristo.IR3.arr1.day))))
ca.day <- data.frame(site = 'Cayo.OR3.day', max = max(Cayo.OR3.arr3.day), min = min(Cayo.OR3.arr3.day), mean = mean(Cayo.OR3.arr3.day), sd = sd(Cayo.OR3.arr3.day), se = (sd(Cayo.OR3.arr3.day)/sqrt(nrow(Cayo.OR3.arr3.day))))
dm.day <- data.frame(site = 'Drago.OR4.day', max = max(Drago.OR4.day), min = min(Drago.OR4.day), mean = mean(Drago.OR4.day), sd = sd(Drago.OR4.day), se = (sd(Drago.OR4.day)/sqrt(nrow(Drago.OR4.day))))

abspar.day <-rbind(pd.day, sp.day, ci.day, ca.day, dm.day)
abspar.day$ann.range <- (abspar.day$max - abspar.day$min)

# site      max    min     mean        sd          se ann.range
# Punta.IR1.day 32.15000 27.038 29.75048 0.8805864 0.006948630  5.112000
# STRI.IR2.day 31.58600 27.038 29.37357 0.8015326 0.006324823  4.548000
# Cristo.IR3.day 33.54778 27.235 29.88817 0.9600723 0.005128799  6.312778
# Cayo.OR3.day 32.43300 27.038 29.56584 0.9798479 0.007731894  5.395000
# Drago.OR4.day 31.66300 26.769 29.16504 0.8567839 0.006760807  4.894000

pd.night <- data.frame(site = 'Punta.IR1.night', max = max(Punta.IR1.night), min = min(Punta.IR1.night), mean = mean(Punta.IR1.night), sd = sd(Punta.IR1.night), se = (sd(Punta.IR1.night)/sqrt(nrow(Punta.IR1.night))))
sp.night <- data.frame(site = 'STRI.IR2.night', max = max(STRI.IR2.night), min = min(STRI.IR2.night), mean = mean(STRI.IR2.night), sd = sd(STRI.IR2.night), se = (sd(STRI.IR2.night)/sqrt(nrow(STRI.IR2.night))))
ci.night <-data.frame(site = 'Cristo.IR3.night', max = max(Cristo.IR3.arr1.night), min = min(Cristo.IR3.arr1.night), mean = mean(Cristo.IR3.arr1.night), sd = sd(Cristo.IR3.arr1.night), se = (sd(Cristo.IR3.arr1.night)/sqrt(nrow(Cristo.IR3.arr1.night))))
ca.night <- data.frame(site = 'Cayo.OR3.night', max = max(Cayo.OR3.arr3.night), min = min(Cayo.OR3.arr3.night), mean = mean(Cayo.OR3.arr3.night), sd = sd(Cayo.OR3.arr3.night), se = (sd(Cayo.OR3.arr3.night)/sqrt(nrow(Cayo.OR3.arr3.night))))
dm.night <- data.frame(site = 'Drago.OR4.night', max = max(Drago.OR4.night), min = min(Drago.OR4.night), mean = mean(Drago.OR4.night), sd = sd(Drago.OR4.night), se = (sd(Drago.OR4.night)/sqrt(nrow(Drago.OR4.night))))

abspar.night <-rbind(pd.night, sp.night, ci.night, ca.night, dm.night)
abspar.night$ann.range <- (abspar.night$max - abspar.night$min)

# site    max      min     mean        sd          se ann.range
# Punta.IR1.night 32.047 27.03800 29.66126 0.8297277 0.006022490  5.009000
# STRI.IR2.night 31.306 27.03800 29.26022 0.7545885 0.005477099  4.268000
# Cristo.IR3.night 32.975 27.30778 29.79554 0.9021838 0.006548404  5.667222
# Cayo.OR3.night 32.150 27.03800 29.39135 0.9157836 0.006647117  5.112000
# Drago.OR4.night 31.306 26.79300 29.09509 0.7906063 0.005738531  4.513000


# looks like the logger absolute data is the same regardless of whether we use daytime only data or not...

# Now look at the satellite absolute data
# this includes dates from 2015-06-30 to 2016-06-30
abspar.sat <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Satellite_Temps/Satellite_All_Data.csv")
str(abspar.sat)
head(abspar.sat)

abspar.sat = abspar.sat %>%
  rename(wirci=X9.27_.82.24,wirpd=X9.37_.82.36,wirpl=X9.13_.82.12,wirsp=X9.35_.82.26,worbs=X9.29_.82.09,
         worca=X9.19_.82.05,wordm=X9.42_.82.32,worbn=X9.35_.82.18) %>%
  slice(2:366)

ci.sat <-data.frame(site = 'wirci', max = max(abspar.sat$wirci), min = min(abspar.sat$wirci), mean = mean(abspar.sat$wirci), sd = sd(abspar.sat$wirci), se = (sd(abspar.sat$wirci)/sqrt(nrow(abspar.sat))))
pd.sat <- data.frame(site = 'wirpd', max = max(abspar.sat$wirpd), min = min(abspar.sat$wirpd), mean = mean(abspar.sat$wirpd), sd = sd(abspar.sat$wirpd), se = (sd(abspar.sat$wirpd)/sqrt(nrow(abspar.sat))))
sp.sat <- data.frame(site = 'wirsp', max = max(abspar.sat$wirsp), min = min(abspar.sat$wirsp), mean = mean(abspar.sat$wirsp), sd = sd(abspar.sat$wirsp), se = (sd(abspar.sat$wirsp)/sqrt(nrow(abspar.sat))))
ca.sat <- data.frame(site = 'worca', max = max(abspar.sat$worca), min = min(abspar.sat$worca), mean = mean(abspar.sat$worca), sd = sd(abspar.sat$worca), se = (sd(abspar.sat$worca)/sqrt(nrow(abspar.sat))))
dm.sat <- data.frame(site = 'wordm', max = max(abspar.sat$wordm), min = min(abspar.sat$wordm), mean = mean(abspar.sat$wordm), sd = sd(abspar.sat$wordm), se = (sd(abspar.sat$wordm)/sqrt(nrow(abspar.sat))))

abspar.sat.all <-rbind(ci.sat, pd.sat, sp.sat, ca.sat, dm.sat)
abspar.sat.all$ann.range <- (abspar.sat.all$max - abspar.sat.all$min)

# site    max    min     mean        sd         se ann.range
# 1 wirci 30.832 27.188 28.93167 0.7733489 0.04047893  3.644001
# 2 wirpd 30.665 26.636 28.82865 0.7728303 0.04045179  4.029001
# 3 wirsp 30.695 26.752 28.83524 0.7799372 0.04082378  3.942999
# 4 worca 30.671 26.861 28.78550 0.7982675 0.04178323  3.809999
# 5 wordm 30.479 26.529 28.72862 0.7829606 0.04098203  3.950001

#### Compiling All Logged Data in Longform ####

# read in this csv file to skip running this code for all temp data
master.melt <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllTempData.csv")

# and for daytime data only
master.melt.day <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllTempData_daytime.csv")

# and for nighttime data only
master.melt.night <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllTempData_nighttime.csv")
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

# merge the sites together for daytime data only
# can only merge two objects at a time
merge.day <- merge(Cayo.OR3.arr3.day,Drago.OR4.day)
merge.day <- merge(merge.day, Punta.IR1.day)
merge.day <- merge(merge.day,STRI.IR2.day)
merged.day <- merge(merge.day, Cristo.IR3.arr1.day)
head(merged.day)

# rename the columns to keep identifying info
colnames(merged.day) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.day <- as.data.frame(merged.day)
# merged is still an xts/zoo object. make it a data.table and keep the datetimes as rownames

merged.day <- setDT(merged.day, keep.rownames = T)
colnames(merged.day) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
head(merged.day)
# melt, format, add in separate date and time columns
master.melt.day<- melt(merged.day, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(master.melt.day) <- c('datetime','logger','temp')
head(master.melt.day)

# master.melt.day is now every daytime measurement you have for each site's "main logger" from the period specified in the initial import as zoo/xts files
# when you converted master.melt to a data.table, the datetime became a character string. this is useful because we can pull date and time out individually.
master.melt.day$time <- substr(master.melt.day$datetime, 12, 19)
master.melt.day$date <- as.Date(master.melt.day$datetime,"%Y-%m-%d" )

# but now put them back in POSIX format so R knows they're dates
# note that when you do this, the time column will be arbitrarily assigned all the same date (current date)
master.melt.day$datetime <- as.POSIXct(strptime(master.melt.day$datetime, "%Y-%m-%d %H:%M:%S"))
master.melt.day$time <- as.POSIXct(strptime(master.melt.day$time, "%H:%M:%S"))
head(master.melt.day)
tail(master.melt.day)

#write.csv(master.melt.day, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/AllTempData_daytime.csv", row.names = FALSE)

# merge the sites together for nighttime data only
# can only merge two objects at a time
merge.night <- merge(Cayo.OR3.arr3.night,Drago.OR4.night)
merge.night <- merge(merge.night, Punta.IR1.night)
merge.night <- merge(merge.night,STRI.IR2.night)
merged.night <- merge(merge.night, Cristo.IR3.arr1.night)
head(merged.night)

# rename the columns to keep identifying info
colnames(merged.night) <- c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
merged.night <- as.data.frame(merged.night)
# merged is still an xts/zoo object. make it a data.table and keep the datetimes as rownames

merged.night <- setDT(merged.night, keep.rownames = T)
colnames(merged.night) <- c('datetime','Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1')
head(merged.night)
# melt, format, add in separate date and time columns
master.melt.night<- melt(merged.night, id.vars = 'datetime', measure.vars = c('Cayo.OR3.arr3','Drago.OR4','Punta.IR1','STRI.IR2','Cristo.IR3.arr1'))
colnames(master.melt.night) <- c('datetime','logger','temp')
head(master.melt.night)

# master.melt.night is now every nighttime measurement you have for each site's "main logger" from the period specified in the initial import as zoo/xts files
# when you converted master.melt to a data.table, the datetime became a character string. this is useful because we can pull date and time out individually.
master.melt.night$time <- substr(master.melt.night$datetime, 12, 19)
master.melt.night$date <- as.Date(master.melt.night$datetime,"%Y-%m-%d" )

# but now put them back in POSIX format so R knows they're dates
# note that when you do this, the time column will be arbitrarily assigned all the same date (current date)
master.melt.night$datetime <- as.POSIXct(strptime(master.melt.night$datetime, "%Y-%m-%d %H:%M:%S"))
master.melt.night$time <- as.POSIXct(strptime(master.melt.night$time, "%H:%M:%S"))
head(master.melt.night)
tail(master.melt.night)

#write.csv(master.melt.night, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/AllTempData_nighttime.csv", row.names = FALSE)


## Summary Stats

aggregate(temp ~ logger, data = master.melt, FUN = "min")

aggregate(temp ~ logger, data = master.melt, FUN = "max")

aggregate(temp ~ logger, data = master.melt, FUN = "mean")

#### Seasonal Parameters ####

# read in these csv files to skip running this section of code for all data, daytime data only, and nighttime data only

dailys <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/SeasonalDailyTempData.csv")
dailys.day <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/SeasonalDailyTempData_daytime.csv")
dailys.night <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/SeasonalDailyTempData_nighttime.csv")
###

# separate into three-month chunks
# in Panama, the rainy season typically extends from about May to December. Jan, Feb, and March are relatively dry; rain starts to pick up in April and May, is in full swing June, July, August, September, October, and November; starts to taper in December
# so in this blocking, Jan-Mar = dry, Apr-Jun = ramp up into rain, Jul-Sep = full on rain, and Oct-Dec = ramp down to dry
# also nice because August was my colony collection month, so I can get the month preceding, month of, and month after


Jul.Aug.Sep <-
  master.melt %>%
  dplyr::filter(datetime >= '2015-07-1 00:00:00', datetime <= '2015-09-30 11:59:59')
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
  dplyr::filter(datetime >= '2016-04-1 00:00:00', datetime <= '2016-06-30 11:59:59')
Apr.May.Jun$season <- 'aprmayjun'

seas <- rbind(Oct.Nov.Dec,Jan.Feb.Mar,Apr.May.Jun,Jul.Aug.Sep)
head(seas)


dailys<-setDT(master.melt)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
dailys$range <- (dailys$max - dailys$min)
dailys$season <- ifelse(dailys$date  >= '2015-07-1' & dailys$date <= '2015-09-30', 'julaugsep',
                        ifelse(dailys$date >= '2015-10-1' & dailys$date <= '2015-12-31', 'octnovdec',
                               ifelse(dailys$date >= '2016-01-1' & dailys$date <= '2016-03-31', 'janfebmar',
                                      ifelse(dailys$date >= '2016-04-1' & dailys$date <= '2016-06-30', 'aprmayjun', 'NA'))))

dailys$season <- factor(dailys$season, levels = c('janfebmar','aprmayjun','julaugsep','octnovdec'))
dailys$logger <- factor(dailys$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))

head(dailys)

#write.csv(dailys, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/SeasonalDailyTempData.csv", row.names = FALSE)

# rank sites by daily variability by season
aggregate(range ~ logger + season, data = dailys, FUN = "mean")
# Order for all seasons : Cristobal, Cayo de Agua, Drago Mar, Punta Donato, STRI

# do this again but for daytime data only
Jul.Aug.Sep.day <-
  master.melt.day %>%
  dplyr::filter(datetime >= '2015-07-1 00:00:00', datetime <= '2015-09-30 11:59:59')
Jul.Aug.Sep.day$season <- 'julaugsep'

Oct.Nov.Dec.day <-
  master.melt.day %>%
  dplyr::filter(datetime >= '2015-10-1 00:00:00', datetime <= '2015-12-31 11:59:59')
Oct.Nov.Dec.day$season <- 'octnovdec'

Jan.Feb.Mar.day <-
  master.melt.day %>%
  dplyr::filter(datetime >= '2016-01-1 00:00:00', datetime <= '2016-03-31 11:59:59')
Jan.Feb.Mar.day$season <- 'janfebmar'

Apr.May.Jun.day <-
  master.melt.day %>%
  dplyr::filter(datetime >= '2016-04-1 00:00:00', datetime <= '2016-06-30 11:59:59')
Apr.May.Jun.day$season <- 'aprmayjun'

seas.day <- rbind(Oct.Nov.Dec.day,Jan.Feb.Mar.day,Apr.May.Jun.day,Jul.Aug.Sep.day)
head(seas.day)


dailys.day<-setDT(master.melt.day)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
dailys.day$range <- (dailys.day$max - dailys.day$min)
dailys.day$season <- ifelse(dailys.day$date  >= '2015-07-1' & dailys.day$date <= '2015-09-30', 'julaugsep',
                            ifelse(dailys.day$date >= '2015-10-1' & dailys.day$date <= '2015-12-31', 'octnovdec',
                                   ifelse(dailys.day$date >= '2016-01-1' & dailys.day$date <= '2016-03-31', 'janfebmar',
                                          ifelse(dailys.day$date >= '2016-04-1' & dailys.day$date <= '2016-06-30', 'aprmayjun', 'NA'))))

dailys.day$season <- factor(dailys.day$season, levels = c('janfebmar','aprmayjun','julaugsep','octnovdec'))
dailys.day$logger <- factor(dailys.day$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))

head(dailys.day)

#write.csv(dailys.day, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/SeasonalDailyTempData_daytime.csv", row.names = FALSE)

# do this again but for nighttime data only
Jul.Aug.Sep.night <-
  master.melt.night %>%
  dplyr::filter(datetime >= '2015-07-1 00:00:00', datetime <= '2015-09-30 11:59:59')
Jul.Aug.Sep.night$season <- 'julaugsep'

Oct.Nov.Dec.night <-
  master.melt.night %>%
  dplyr::filter(datetime >= '2015-10-1 00:00:00', datetime <= '2015-12-31 11:59:59')
Oct.Nov.Dec.night$season <- 'octnovdec'

Jan.Feb.Mar.night <-
  master.melt.night %>%
  dplyr::filter(datetime >= '2016-01-1 00:00:00', datetime <= '2016-03-31 11:59:59')
Jan.Feb.Mar.night$season <- 'janfebmar'

Apr.May.Jun.night <-
  master.melt.night %>%
  dplyr::filter(datetime >= '2016-04-1 00:00:00', datetime <= '2016-06-30 11:59:59')
Apr.May.Jun.night$season <- 'aprmayjun'

seas.night <- rbind(Oct.Nov.Dec.night,Jan.Feb.Mar.night,Apr.May.Jun.night,Jul.Aug.Sep.night)
head(seas.night)


dailys.night<-setDT(master.melt.night)[, list(max = max(temp), min = min(temp), mean=mean(temp)), by=list(logger,date)]
dailys.night$range <- (dailys.night$max - dailys.night$min)
dailys.night$season <- ifelse(dailys.night$date  >= '2015-07-1' & dailys.night$date <= '2015-09-30', 'julaugsep',
                              ifelse(dailys.night$date >= '2015-10-1' & dailys.night$date <= '2015-12-31', 'octnovdec',
                                     ifelse(dailys.night$date >= '2016-01-1' & dailys.night$date <= '2016-03-31', 'janfebmar',
                                            ifelse(dailys.night$date >= '2016-04-1' & dailys.night$date <= '2016-06-30', 'aprmayjun', 'NA'))))

dailys.night$season <- factor(dailys.night$season, levels = c('janfebmar','aprmayjun','julaugsep','octnovdec'))
dailys.night$logger <- factor(dailys.night$logger, levels = c('Punta.IR1', 'STRI.IR2', 'Cristo.IR3.arr1', 'Cayo.OR3.arr3', 'Drago.OR4'))

head(dailys.night)

#write.csv(dailys.night, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/SeasonalDailyTempData_nighttime.csv", row.names = FALSE)

# rank sites by daily variability by season
aggregate(range ~ logger + season, data = dailys.night, FUN = "mean")
# Order Jan Feb Mar & Apr May Jun: Cristobal, Cayo de Agua, Drago Mar, Punta Donato, STRI
# Order Jul Aug Sep & Oct Nov Dec: Cayo de Agua, Cristobal, Drago Mar, Punta Donato, STRI

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
#write.csv(daily.mmm, file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/DailyTempRangeData.csv", row.names = FALSE)

daily.mmm = daily.mmm %>%
  dplyr::filter(logger!="Drago.OR4")

# rank sites by daily variability overall
aggregate(range ~ logger, data = daily.mmm, FUN = "mean")

#       Punta.IR1 0.7041393 (4)
#        STRI.IR2 0.5603306 (5)
# Cristo.IR3.arr1 1.1835231 (1)
#   Cayo.OR3.arr3 0.9772268 (2)

aggregate(range ~ logger, data = daily.mmm, FUN = "max")

#   Cayo.OR3.arr3 2.411000
# Cristo.IR3.arr1 2.856111
#       Punta.IR1 1.315000
#        STRI.IR2 1.263000

summarySE(data = daily.mmm, groupvar = "logger", measurevar = "mean")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "min")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "max")
summarySE(data = daily.mmm, groupvar = "logger", measurevar = "range")


#### Correlate logger weekly temperature range with satellite weekly range ####

satellite_weekly <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Satellite_Temps/Satellite_Weekly_Range.csv")
head(satellite_weekly) # dates range 2015-06-14 to 2016-08-07

# make sure R knows this is a date
satellite_weekly$date <- as.POSIXct(strptime(satellite_weekly$date, "%Y-%m-%d"))

# re-name columns
satellite_weekly = satellite_weekly %>%
  rename(wirci=X9.27_.82.24,wirpd=X9.37_.82.36,wirpl=X9.13_.82.12,wirsp=X9.35_.82.26,worbs=X9.29_.82.09,
         worca=X9.19_.82.05,wordm=X9.42_.82.32,worbn=X9.35_.82.18)

head(satellite_weekly)

# filter dates to match loggers
satellite_weekly <- satellite_weekly[4:56, ]
head(satellite_weekly)

satellite.weekly.melt<- melt(satellite_weekly, id.vars = 'date', measure.vars = c("wirci","wirpd","wirpl","wirsp","worbs","worca","wordm","worbn"))
head(satellite.weekly.melt)

colnames(satellite.weekly.melt) <- c('date','site','range')
head(satellite.weekly.melt)
str(satellite.weekly.melt)

satellite.weekly.melt$month <- ifelse(satellite.weekly.melt$date  >= '2015-07-05' & satellite.weekly.melt$date <= '2015-07-26', 'jul',
                                      ifelse(satellite.weekly.melt$date >= '2015-08-02' & satellite.weekly.melt$date <= '2015-08-30', 'aug',
                                             ifelse(satellite.weekly.melt$date >= '2015-09-06' & satellite.weekly.melt$date <= '2015-09-27', 'sep',
                                                    ifelse(satellite.weekly.melt$date >= '2015-10-04' & satellite.weekly.melt$date <= '2015-10-25', 'oct',
                                                           ifelse(satellite.weekly.melt$date >= '2015-11-01' & satellite.weekly.melt$date <= '2015-11-29', 'nov',
                                                                  ifelse(satellite.weekly.melt$date >= '2015-12-06' & satellite.weekly.melt$date <= '2015-12-27', 'dec',
                                                                         ifelse(satellite.weekly.melt$date >= '2016-01-03' & satellite.weekly.melt$date <= '2016-01-31', 'jan',
                                                                                ifelse(satellite.weekly.melt$date >= '2016-02-07' & satellite.weekly.melt$date <= '2016-02-28', 'feb',
                                                                                       ifelse(satellite.weekly.melt$date >= '2016-03-06' & satellite.weekly.melt$date <= '2016-03-27', 'mar',
                                                                                              ifelse(satellite.weekly.melt$date >= '2016-04-03' & satellite.weekly.melt$date <= '2016-04-24', 'apr',
                                                                                                     ifelse(satellite.weekly.melt$date >= '2016-05-01' & satellite.weekly.melt$date <= '2016-05-29', 'may',
                                                                                                            ifelse(satellite.weekly.melt$date >= '2016-06-05' & satellite.weekly.melt$date <= '2016-06-26', 'jun',
                                                                                                                   ifelse(satellite.weekly.melt$date == '2016-07-03', 'jul',
                                                                                                                          'NA')))))))))))))

satellite.weekly.melt$season <- ifelse(satellite.weekly.melt$date >= '2015-07-05' & satellite.weekly.melt$date <= '2015-09-27', 'julaugsep',
                                       ifelse(satellite.weekly.melt$date >= '2015-10-04' & satellite.weekly.melt$date <= '2015-12-27', 'octnovdec',
                                              ifelse(satellite.weekly.melt$date >= '2016-01-03' & satellite.weekly.melt$date <= '2016-03-27', 'janfebmar',
                                                     ifelse(satellite.weekly.melt$date >= '2016-04-03' & satellite.weekly.melt$date <= '2016-06-26', 'aprmayjun',
                                                            ifelse(satellite.weekly.melt$date == '2016-07-03', 'julaugsep',
                                                                   'NA')))))

satellite.weekly.melt$season <- as.factor(satellite.weekly.melt$season)
satellite.weekly.melt$month <- factor(satellite.weekly.melt$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

aggregate(range ~ site, data = satellite.weekly.melt, FUN = "mean")

# now logger data, either subset by daytime temps only or not
logger_weekly <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams.csv")
logger_weekly_daytime <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams_daytime.csv")
logger_weekly_nighttime <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/WeeklyTempParams_nighttime.csv")

logger_weekly$season <- as.factor(logger_weekly$season)
logger_weekly$logger <- as.factor(logger_weekly$logger)
logger_weekly$month <- factor(logger_weekly$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

logger_weekly_daytime$season <- as.factor(logger_weekly_daytime$season)
logger_weekly_daytime$logger <- as.factor(logger_weekly_daytime$logger)
logger_weekly_daytime$month <- factor(logger_weekly_daytime$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

logger_weekly_nighttime$season <- as.factor(logger_weekly_nighttime$season)
logger_weekly_nighttime$logger <- as.factor(logger_weekly_nighttime$logger)
logger_weekly_nighttime$month <- factor(logger_weekly_nighttime$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

str(logger_weekly) # dates range 2015-07-05 to 2016-06-30
str(logger_weekly_daytime)
str(logger_weekly_nighttime)
str(satellite.weekly.melt)

# remove  satellite info that we don't have logger info for
satellite_weekly_tocombine =
  satellite.weekly.melt %>%
  dplyr::filter(site != "wirpl") %>% # remove Punta Laurel
  dplyr::filter(site != "worbs") %>% # remove Bastimentos South
  dplyr::filter(site != "worbn") %>% # remove Bastimentos North
  mutate(type = "satellite")
str(satellite_weekly_tocombine)
#write.csv(satellite_weekly_tocombine, "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/satellite_weekly_tocombine.csv")

# rename loggers so we can combine the data
logger_weekly_tocombine =
  logger_weekly %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_weekly_tocombine)
#write.csv(logger_weekly_tocombine, "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/logger_weekly_tocombine_alldata.csv")

logger_weekly_daytime_tocombine =
  logger_weekly_daytime %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_weekly_daytime_tocombine)

logger_weekly_nighttime_tocombine =
  logger_weekly_nighttime %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_weekly_nighttime_tocombine)

# subset weekly logger data by logger name to facilitate correlations
Cayo_logger <- subset(logger_weekly_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger <- subset(logger_weekly_tocombine, logger=="Punta.IR1")
STRI_logger <- subset(logger_weekly_tocombine, logger=="STRI.IR2")
Cristo_logger <- subset(logger_weekly_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger <- subset(logger_weekly_tocombine, logger=="Drago.OR4")

Cayo_logger_daytime <- subset(logger_weekly_daytime_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger_daytime <- subset(logger_weekly_daytime_tocombine, logger=="Punta.IR1")
STRI_logger_daytime <- subset(logger_weekly_daytime_tocombine, logger=="STRI.IR2")
Cristo_logger_daytime <- subset(logger_weekly_daytime_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger_daytime <- subset(logger_weekly_daytime_tocombine, logger=="Drago.OR4")

Cayo_logger_nighttime <- subset(logger_weekly_nighttime_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger_nighttime <- subset(logger_weekly_nighttime_tocombine, logger=="Punta.IR1")
STRI_logger_nighttime <- subset(logger_weekly_nighttime_tocombine, logger=="STRI.IR2")
Cristo_logger_nighttime <- subset(logger_weekly_nighttime_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger_nighttime <- subset(logger_weekly_nighttime_tocombine, logger=="Drago.OR4")

Cayo_satellite <- subset(satellite_weekly_tocombine, site=="worca")
Punta_satellite <- subset(satellite_weekly_tocombine, site=="wirpd")
STRI_satellite <- subset(satellite_weekly_tocombine, site=="wirsp")
Cristo_satellite <- subset(satellite_weekly_tocombine, site=="wirci")
Drago_satelite <- subset(satellite_weekly_tocombine, site=="wordm")

# Satellite data names: wirci wirpd wirpl wirsp worbs worca wordm worbn
# Logger data names: Punta.IR1 STRI.IR2 Cristo.IR3.arr1 Cayo.OR3.arr3 Drago.OR4

cor.test(Cayo_logger$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger$range, Punta_satellite$range, method = "pearson") # Punta Donato
cor.test(STRI_logger$range, STRI_satellite$range, method = "pearson") #STRI
cor.test(Cristo_logger$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger$range, Drago_satelite$range, method = "pearson") #Drago Mar

cor.test(Cayo_logger_nighttime$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger_nighttime$range, Punta_satellite$range, method = "pearson") # Punta Donato, sig
cor.test(STRI_logger_nighttime$range, STRI_satellite$range, method = "pearson") #STRI, sig
cor.test(Cristo_logger_nighttime$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger_nighttime$range, Drago_satelite$range, method = "pearson") #Drago Mar

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/WeeklyTempCorrelations_nighttime.pdf",
    width = 6,
    height = 6)
par(mfrow=c(3,2))
plot(Cayo_logger_nighttime$range, Cayo_satellite$range,
     main = "Cayo de Agua",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, col="black", frame = TRUE) + # Cayo de Agua
  abline(lm(Cayo_satellite$range~Cayo_logger_nighttime$range))

plot(Punta_logger_nighttime$range, Punta_satellite$range, main = "Punta Donato, *p=0.01*",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Punta_satellite$range~Punta_logger_nighttime$range)) # Punta Donato

plot(STRI_logger_nighttime$range, STRI_satellite$range, main = "STRI, *p<0.01*",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(STRI_satellite$range~STRI_logger_nighttime$range)) #STRI

plot(Cristo_logger_nighttime$range, Cristo_satellite$range, main = "Cristobal Island",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Cristo_satellite$range~Cristo_logger_nighttime$range)) #Cristobal

plot(Drago_logger_nighttime$range, Drago_satelite$range, main = "Drago Mar",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) + # Cayo de Agua
  abline(lm(Drago_satelite$range~Drago_logger_nighttime$range)) #Drago Mar
dev.off()


#### Order sites by logger and satellite weekly range ####
# this is looking at all weekly range data, not subset by month or season, for loggers (both daytime only and all data) and satellite
mean(Cayo_logger$range) #1.927038 (2)
mean(Punta_logger$range) #1.477019 (4)
mean(STRI_logger$range) #1.132623 (5)
mean(Cristo_logger$range) #2.380881 (1)
mean(Drago_logger$range) #1.688396 (3)

mean(Cayo_logger_daytime$range) #1.869604 (2)
mean(Punta_logger_daytime$range) #0.9109582 (5)
mean(STRI_logger_daytime$range) #1.105906 (4)
mean(Cristo_logger_daytime$range) #2.303218 (1)
mean(Drago_logger_daytime$range) #1.562075 (3)

mean(Cayo_logger_nighttime$range) #1.742302 (2)
mean(Punta_logger_nighttime$range) #0.7261651 (5)
mean(STRI_logger_nighttime$range) #0.8523774 (4)
mean(Cristo_logger_nighttime$range) #1.913606 (1)
mean(Drago_logger_nighttime$range) #1.473075 (3)

mean(Cayo_satellite$range) #0.7945094 (2)
mean(Punta_satellite$range) #0.7899812 (4)
mean(STRI_satellite$range) #0.79383 (3)
mean(Cristo_satellite$range) #0.816 (1)
mean(Drago_satelite$range) #0.7449434 (5)

# so there are hardly any differences in mean range across the satellite, but the logger is capturing this much better

# sarah wants to correlate the sites by mean weekly ranges
# put all of the ranges above into a spreadsheet by hand because i am lazy.
weekranges <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/AllWeekRanges.csv")
str(weekranges)

weekranges$site <- as.factor(weekranges$site)
weekranges$type <- as.factor(weekranges$type)

satellite <- subset(weekranges, type=="sat" & site!="drago")
logger <- subset(weekranges, type=="logger"  & site!="drago")
logger_day <- subset(weekranges, type=="logger_day"  & site!="drago")
logger_night <- subset(weekranges, type=="logger_night"  & site!="drago")

cor.test(satellite$avg_range, logger$avg_range, method = "pearson")
cor.test(satellite$avg_range, logger_day$avg_range, method = "pearson")
cor.test(satellite$avg_range, logger_night$avg_range, method = "pearson")

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/AvgWeeklyTempCorrelations.pdf",
    width = 10,
    height = 4)
par(mfrow=c(1,3))

plot(logger$avg_range, satellite$avg_range,
     main = "Average Weekly Range \n Satellite vs. Logger",
     xlab = "Logger Weekly Range", ylab = "Satellite Weekly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger$avg_range)) +
  text(satellite$avg_range~logger$avg_range, labels=logger$site, cex=0.9, font=2)


plot(logger_day$avg_range, satellite$avg_range, main = "Average Weekly Range \n Satellite vs. Daytime Logger",
     xlab = "Logger Weekly Range, Daytime", ylab = "Satellite Weekly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger_day$avg_range)) +
  text(satellite$avg_range~logger_day$avg_range, labels=logger_day$site, cex=0.9, font=2)

plot(logger_night$avg_range, satellite$avg_range, main = "Average Weekly Range \n Satellite vs. Nighttime Logger",
     xlab = "Logger Weekly Range, Nighttime", ylab = "Satellite Weekly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger_night$avg_range)) +
  text(satellite$avg_range~logger_night$avg_range, labels=logger_night$site, cex=0.9, font=2)

dev.off()

#### Correlate monthly logger and monthly satellite data ####

satellite_monthly <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Satellite_Temps/Satellite_Monthly_Range.csv")
head(satellite_monthly) # dates range 2015-06-30 to 2016-08-31
# loggers range 2015-07-31 to 2016-06-29
# make sure R knows this is a date
satellite_monthly$date <- as.POSIXct(strptime(satellite_monthly$date, "%Y-%m-%d"))

# re-name columns
satellite_monthly = satellite_monthly %>%
  rename(wirci=X9.27_.82.24,wirpd=X9.37_.82.36,wirpl=X9.13_.82.12,wirsp=X9.35_.82.26,worbs=X9.29_.82.09,
         worca=X9.19_.82.05,wordm=X9.42_.82.32,worbn=X9.35_.82.18) %>%
  slice(2:13) # filter dates to match loggers


head(satellite_monthly)

satellite.monthly.melt<- melt(satellite_monthly, id.vars = 'date', measure.vars = c("wirci","wirpd","wirpl","wirsp","worbs","worca","wordm","worbn"))
colnames(satellite.monthly.melt) <- c('date','site','range')

# now logger data, either subset by daytime temps only or not
logger_monthly <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/MonthlyTempParams.csv")
logger_monthly_daytime <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/MonthlyTempParams_daytime.csv")
logger_monthly_nighttime <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/MonthlyTempParams_nighttime.csv")

logger_monthly$season <- as.factor(logger_monthly$season)
logger_monthly$logger <- as.factor(logger_monthly$logger)
logger_monthly$month <- factor(logger_monthly$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

logger_monthly_daytime$season <- as.factor(logger_monthly_daytime$season)
logger_monthly_daytime$logger <- as.factor(logger_monthly_daytime$logger)
logger_monthly_daytime$month <- factor(logger_monthly_daytime$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

logger_monthly_nighttime$season <- as.factor(logger_monthly_nighttime$season)
logger_monthly_nighttime$logger <- as.factor(logger_monthly_nighttime$logger)
logger_monthly_nighttime$month <- factor(logger_monthly_nighttime$month, levels=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))

str(logger_monthly) # dates range 2015-07-05 to 2016-06-30
str(logger_monthly_daytime)
str(logger_monthly_nighttime)
str(satellite.monthly.melt)


# remove  satellite info that we don't have logger info for
satellite_monthly_tocombine =
  satellite.monthly.melt %>%
  dplyr::filter(site != "wirpl") %>% # remove Punta Laurel
  dplyr::filter(site != "worbs") %>% # remove Bastimentos South
  dplyr::filter(site != "worbn") %>% # remove Bastimentos North
  mutate(type = "satellite")
str(satellite_monthly_tocombine)

# rename loggers so we can combine the data
logger_monthly_tocombine =
  logger_monthly %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_monthly_tocombine)

logger_monthly_daytime_tocombine =
  logger_monthly_daytime %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_monthly_daytime_tocombine)

logger_monthly_nighttime_tocombine =
  logger_monthly_nighttime %>%
  mutate(site = recode_factor(logger, "Cayo.OR3.arr3"="wirca", "Cristo.IR3.arr1"="wirci", "Drago.OR4"="wordm", "Punta.IR1"="wirpd", "STRI.IR2"="wirsp")) %>%
  mutate(type = "logger")
str(logger_monthly_nighttime_tocombine)

# subset monthly logger data by logger name to facilitate correlations
Cayo_logger <- subset(logger_monthly_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger <- subset(logger_monthly_tocombine, logger=="Punta.IR1")
STRI_logger <- subset(logger_monthly_tocombine, logger=="STRI.IR2")
Cristo_logger <- subset(logger_monthly_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger <- subset(logger_monthly_tocombine, logger=="Drago.OR4")

Cayo_logger_daytime <- subset(logger_monthly_daytime_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger_daytime <- subset(logger_monthly_daytime_tocombine, logger=="Punta.IR1")
STRI_logger_daytime <- subset(logger_monthly_daytime_tocombine, logger=="STRI.IR2")
Cristo_logger_daytime <- subset(logger_monthly_daytime_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger_daytime <- subset(logger_monthly_daytime_tocombine, logger=="Drago.OR4")

Cayo_logger_nighttime <- subset(logger_monthly_nighttime_tocombine, logger=="Cayo.OR3.arr3")
Punta_logger_nighttime <- subset(logger_monthly_nighttime_tocombine, logger=="Punta.IR1")
STRI_logger_nighttime <- subset(logger_monthly_nighttime_tocombine, logger=="STRI.IR2")
Cristo_logger_nighttime <- subset(logger_monthly_nighttime_tocombine, logger=="Cristo.IR3.arr1")
Drago_logger_nighttime <- subset(logger_monthly_nighttime_tocombine, logger=="Drago.OR4")

Cayo_satellite <- subset(satellite_monthly_tocombine, site=="worca")
Punta_satellite <- subset(satellite_monthly_tocombine, site=="wirpd")
STRI_satellite <- subset(satellite_monthly_tocombine, site=="wirsp")
Cristo_satellite <- subset(satellite_monthly_tocombine, site=="wirci")
Drago_satelite <- subset(satellite_monthly_tocombine, site=="wordm")

# Satellite data names: wirci wirpd wirpl wirsp worbs worca wordm worbn
# Logger data names: Punta.IR1 STRI.IR2 Cristo.IR3.arr1 Cayo.OR3.arr3 Drago.OR4

cor.test(Cayo_logger$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger$range, Punta_satellite$range, method = "pearson") # Punta Donato
cor.test(STRI_logger$range, STRI_satellite$range, method = "pearson") #STRI
cor.test(Cristo_logger$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger$range, Drago_satelite$range, method = "pearson") #Drago Mar

cor.test(Cayo_logger_daytime$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger_daytime$range, Punta_satellite$range, method = "pearson") # Punta Donato
cor.test(STRI_logger_daytime$range, STRI_satellite$range, method = "pearson") #STRI
cor.test(Cristo_logger_daytime$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger_daytime$range, Drago_satelite$range, method = "pearson") #Drago Mar

cor.test(Cayo_logger_nighttime$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger_nighttime$range, Punta_satellite$range, method = "pearson") # Punta Donato
cor.test(STRI_logger_nighttime$range, STRI_satellite$range, method = "pearson") #STRI
cor.test(Cristo_logger_nighttime$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger_nighttime$range, Drago_satelite$range, method = "pearson") #Drago Mar

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/MonthlyTempCorrelations_nighttime.pdf",
    width = 6,
    height = 6)
par(mfrow=c(3,2))
plot(Cayo_logger_nighttime$range, Cayo_satellite$range, main = "Cayo de Agua",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, frame = TRUE) + # Cayo de Agua
  abline(lm(Cayo_satellite$range~Cayo_logger_nighttime$range))

plot(Punta_logger_nighttime$range, Punta_satellite$range, main = "Punta Donato",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Punta_satellite$range~Punta_logger_nighttime$range)) # Punta Donato

plot(STRI_logger_nighttime$range, STRI_satellite$range, main = "STRI",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, frame = TRUE) +
  abline(lm(STRI_satellite$range~STRI_logger_nighttime$range)) #STRI

plot(Cristo_logger_nighttime$range, Cristo_satellite$range, main = "Cristobal Island",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Cristo_satellite$range~Cristo_logger_nighttime$range)) #Cristobal

plot(Drago_logger_nighttime$range, Drago_satelite$range, main = "Drago Mar",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, frame = TRUE) + # Cayo de Agua
  abline(lm(Drago_satelite$range~Drago_logger_nighttime$range)) #Drago Mar
dev.off()

# compare mean monthly temperatures across satellites and loggers
# rank sites by daily variability by season
aggregate(range ~ logger, data = logger_monthly, FUN = "mean")
aggregate(range ~ logger, data = logger_monthly_daytime, FUN = "mean")
aggregate(range ~ logger, data = logger_monthly_nighttime, FUN = "mean")
aggregate(range ~ site, data = satellite.monthly.melt, FUN = "mean")

# sarah wants to correlate the sites by mean monthly ranges
# put all of the ranges above into a spreadsheet by hand because i am lazy.
monthranges <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/AllMonthRanges.csv")
str(monthranges)

monthranges$site <- as.factor(monthranges$site)
monthranges$type <- as.factor(monthranges$type)

satellite <- subset(monthranges, type=="sat" & site!="drago")
logger <- subset(monthranges, type=="logger"  & site!="drago")
logger_day <- subset(monthranges, type=="logger_day"  & site!="drago")
logger_night <- subset(monthranges, type=="logger_night"  & site!="drago")

cor.test(satellite$avg_range, logger$avg_range, method = "pearson")
cor.test(satellite$avg_range, logger_day$avg_range, method = "pearson")
cor.test(satellite$avg_range, logger_night$avg_range, method = "pearson")

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/AvgMonthlyTempCorrelations.pdf",
    width = 10,
    height = 4)
par(mfrow=c(1,3))

plot(logger$avg_range, satellite$avg_range,
     main = "Average Monthly Range \n Satellite vs. Logger",
     xlab = "Logger Monthly Range", ylab = "Satellite Monthly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger$avg_range)) +
  text(satellite$avg_range~logger$avg_range, labels=logger$site, cex=0.9, font=2)

plot(logger_day$avg_range, satellite$avg_range, main = "Average Monthly Range \n Satellite vs. Daytime Logger",
     xlab = "Logger Monthly Range, Daytime", ylab = "Satellite Monthly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger_day$avg_range)) +
  text(satellite$avg_range~logger_day$avg_range, labels=logger_day$site, cex=0.9, font=2)

plot(logger_night$avg_range, satellite$avg_range, main = "Average Monthly Range \n Satellite vs. Nighttime Logger",
     xlab = "Logger Monthly Range, Nighttime", ylab = "Satellite Monthly Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$avg_range~logger_night$avg_range)) +
  text(satellite$avg_range~logger_night$avg_range, labels=logger_night$site, cex=0.9, font=2)

dev.off()

#### Correlate annual ranges from satellite and loggers ####

absranges <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Hobo_Loggers/data_sheets/AbsoluteTempParams.csv")
str(absranges)

absranges$site <- as.factor(absranges$site)
absranges$type <- as.factor(absranges$type)

satellite <- subset(absranges, type=="sat" & site!="drago")
logger <- subset(absranges, type=="logger"  & site!="drago")
logger_day <- subset(absranges, type=="logger_day"  & site!="drago")
logger_night <- subset(absranges, type=="logger_night"  & site!="drago")

# values for daytime and all data loggers are the same - so only doing the entire logger data for correlations and plots
cor.test(satellite$ann.range, logger$ann.range, method = "pearson")
cor.test(satellite$max, logger$max, method = "pearson")
cor.test(satellite$ann.range, logger_night$ann.range, method = "pearson")
cor.test(satellite$max, logger_night$max, method = "pearson")

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/AbsoluteTempCorrelations.pdf",
    width = 7,
    height = 7)
par(mfrow=c(2,2))

plot(logger$ann.range, satellite$ann.range,
     main = "Absolute Range \n Satellite vs. Logger",
     xlab = "Logger Range", ylab = "Satellite Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$ann.range~logger$ann.range)) +
  text(satellite$ann.range~logger$ann.range, labels=logger$site, cex=0.9, font=2)

plot(logger$max, satellite$max,
     main = "Absolute Max Temp \n Satellite vs. Logger",
     xlab = "Logger Max Temp", ylab = "Satellite Max Temp",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$max~logger$max)) +
  text(satellite$max~logger$max, labels=logger$site, cex=0.9, font=2)

plot(logger_night$ann.range, satellite$ann.range,
     main = "Absolute Range \n Satellite vs. Nighttime Logger",
     xlab = "Logger Range, Nighttime", ylab = "Satellite Range",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$ann.range~logger_night$ann.range)) +
  text(satellite$ann.range~logger_night$ann.range, labels=logger_night$site, cex=0.9, font=2)

plot(logger_night$max, satellite$max,
     main = "Absolute Max Temp \n Satellite vs. Nighttime Logger",
     xlab = "Logger Max Temp, Nighttime", ylab = "Satellite Max Temp",
     pch = 19, col="lightblue", frame = TRUE) +
  abline(lm(satellite$max~logger_night$max)) +
  text(satellite$max~logger_night$max, labels=logger_night$site, cex=0.9, font=2)

dev.off()

#### Correlate logger daily temperature range with satellite weekly temperature range ####
Cayo_logger_daily <- subset(daily.mmm, logger=="Cayo.OR3.arr3")
Punta_logger_daily <- subset(daily.mmm, logger=="Punta.IR1")
STRI_logger_daily <- subset(daily.mmm, logger=="STRI.IR2")
Cristo_logger_daily <- subset(daily.mmm, logger=="Cristo.IR3.arr1")
Drago_logger_daily <- subset(daily.mmm, logger=="Drago.OR4")

# need to reduce the daily logger data to be the same size as the satellite data, 53 rows
Cayo_logger_daily_reduced <- Cayo_logger_daily[sample(nrow(Cayo_logger_daily), 53), ]
Punta_logger_daily_reduced <- Punta_logger_daily[sample(nrow(Punta_logger_daily), 53), ]
STRI_logger_daily_reduced <- STRI_logger_daily[sample(nrow(STRI_logger_daily), 53), ]
Cristo_logger_daily_reduced <- Cristo_logger_daily[sample(nrow(Cristo_logger_daily), 53), ]
Drago_logger_daily_reduced <- Drago_logger_daily[sample(nrow(Drago_logger_daily), 53), ]

cor.test(Cayo_logger_daily_reduced$range, Cayo_satellite$range, method = "pearson")
cor.test(Punta_logger_daily_reduced$range, Punta_satellite$range, method = "pearson") # Punta Donato
cor.test(STRI_logger_daily_reduced$range, STRI_satellite$range, method = "pearson") #STRI
cor.test(Cristo_logger_daily_reduced$range, Cristo_satellite$range, method = "pearson") #Cristobal
cor.test(Drago_logger_daily_reduced$range, Drago_satelite$range, method = "pearson") #Drago Mar

plot(Cayo_logger_daily_reduced$range, Cayo_satellite$range)
plot(Punta_logger_daily_reduced$range, Punta_satellite$range)
plot(STRI_logger_daily_reduced$range, STRI_satellite$range)
plot(Cristo_logger_daily_reduced$range, Cristo_satellite$range)
plot(Drago_logger_daily_reduced$range, Drago_satelite$range)

pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Logger_Sat_Comparisons/WeeklySatvsDailyLoggerTempCorrelations_all.pdf",
    width = 6,
    height = 6)
par(mfrow=c(3,2))
plot(Cayo_logger_daily_reduced$range, Cayo_satellite$range, main = "Cayo de Agua",
     xlab = "Logger Daily Range", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) + # Cayo de Agua
  abline(lm(Cayo_satellite$range~Cayo_logger_daily_reduced$range))

plot(Punta_logger_daily_reduced$range, Punta_satellite$range, main = "Punta Donato",
     xlab = "Logger Daily Range", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Punta_satellite$range~Punta_logger_daily_reduced$range)) # Punta Donato

plot(STRI_logger_daily_reduced$range, STRI_satellite$range, main = "STRI",
     xlab = "Logger Daily Range", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(STRI_satellite$range~STRI_logger_daily_reduced$range)) #STRI

plot(Cristo_logger_daily_reduced$range, Cristo_satellite$range, main = "Cristobal Island",
     xlab = "Logger Daily Range", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) +
  abline(lm(Cristo_satellite$range~Cristo_logger_daily_reduced$range)) #Cristobal

plot(Drago_logger_daily_reduced$range, Drago_satelite$range, main = "Drago Mar",
     xlab = "Logger Daily Range", ylab = "Satellite Weekly Range",
     pch = 19, frame = TRUE) + # Cayo de Agua
  abline(lm(Drago_satelite$range~Drago_logger_daily_reduced$range)) #Drago Mar
dev.off()

#### 90% Daily Thermal Range ####

# Calculate 90% quantile of the daily temperature range (Kenkel et al., Ecology)
# an estimate of the magnitude of high-frequency temp fluctuations
# we can use the quantile() function to calculate this

# daily variability logger data
daily.mmm %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
#  Punta.IR1           0.999
#  STRI.IR2            0.833
#  Cristo.IR3.arr1     1.92
#  Cayo.OR3.arr3       1.49
#  Drago.OR4           1.24


# weekly variability logger data, all data points included to calculate weekly variability
head(weekly.mmm)

weekly.mmm %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
#  Punta.IR1            1.90
#  STRI.IR2             1.49
#  Cristo.IR3.arr1      3.14
#  Cayo.OR3.arr3        2.55
#  Drago.OR4            2.27

# weekly variability logger data, only daylight data points included too calculate weekly variability
head(weekly.mmm.day)

weekly.mmm.day %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
#  Punta.IR1            1.28
#  STRI.IR2             1.48
#  Cristo.IR3.arr1      3.09
#  Cayo.OR3.arr3        2.45
#  Drago.OR4            2.07

weekly.mmm.night %>%
  group_by(logger) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# logger          percent90
#  Punta.IR1            1.08
#  STRI.IR2             1.24
#  Cristo.IR3.arr1      2.64
#  Cayo.OR3.arr3        2.34
#  Drago.OR4            2.07

head(satellite.weekly.melt)

satellite.weekly.melt %>%
  group_by(site) %>%
  summarise(percent90 = quantile(range, probs = 0.9))

# site  percent90
#  wirci      1.40 - Cristobal Island (CI)
#  worca      1.27 - Cayo de Agua (CA)
#  wirsp      1.34 - STRI Point (SP)
#  wirpd      1.27 - Punta Donato (PD)

#  wirpl      1.32 - Punta Laurel (NA)
#  worbs      1.31 - Bastimentos South (NA)
#  wordm      1.27 - Drago Mar (DM)
#  worbn      1.23 - Bastimentos North (NA)

#### Plots ####

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

# present this data in a slightly different way
daily.mmm.plot = daily.mmm %>%
  dplyr::filter(logger != "Drago.OR4")

str(daily.mmm.plot)

levels(daily.mmm.plot$logger)

#### make figures for paper ####
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

ggsave(mean_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/DailyMeanTemp_boxplot.pdf", width=3.5, height=4, units=c("in"), useDingbats=FALSE)

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
  ylim(0,3) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_bw() +
  theme(legend.position = "none")
dtv_boxplot

ggsave(dtv_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/DailyTempRange_boxplot.pdf", width=3.5, height=4, units=c("in"), useDingbats=FALSE)

# stats by logger for mean and dtv
aov1 = aov(range ~ logger, data = daily.mmm.plot)
summary(aov1)
#               Df Sum Sq Mean Sq F value Pr(>F)
# logger         3  85.08  28.359   204.6 <2e-16 ***
# Residuals   1460 202.35   0.139

TukeyHSD(aov1)
# $logger
#           diff       lwr       upr   p adj
# PD-SP 0.1438087 0.0730274 0.2145901 1.2e-06
# CA-SP 0.4168962 0.3461148 0.4876775 0.0e+00
# CI-SP 0.6231925 0.5524111 0.6939738 0.0e+00
# CA-PD 0.2730874 0.2023061 0.3438688 0.0e+00
# CI-PD 0.4793837 0.4086024 0.5501651 0.0e+00
# CI-CA 0.2062963 0.1355150 0.2770776 0.0e+00

aov2 = aov(mean ~ logger, data = daily.mmm.plot)
summary(aov2)
#               Df Sum Sq Mean Sq F value Pr(>F)
# logger         3   70.8  23.613   32.83 <2e-16 ***
# Residuals   1460 1050.1   0.719

TukeyHSD(aov2)
# $logger
#             diff          lwr         upr     p adj
# PD-SP  0.3901668  0.228923474  0.55141020 0.0000000
# CA-SP  0.1596180 -0.001625365  0.32086136 0.0535289
# CI-SP  0.5775439  0.416300570  0.73878730 0.0000000
# CA-PD -0.2305488 -0.391792202 -0.06930548 0.0013927
# CI-PD  0.1873771  0.026133733  0.34862046 0.0150881
# CI-CA  0.4179259  0.256682572  0.57916930 0.0000000

# find maximum DTV for each site
daily.mmm.plot %>%
  group_by(logger) %>%
  summarise_each(funs(mean))

## now dot plots of the same data

#SummarySE to format data for plotting
min_means <- summarySE(daily.mmm.plot, measurevar="min", groupvars=c("logger"))

# plot, minimum temp
min_temp_plot <- ggplot(min_means,aes(x = logger, y = min, color = logger, pch = logger))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = logger, ymax = min+se, ymin = min-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values = c("lightblue","red4","indianred3","mistyrose3"))+
  scale_shape_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values=c(19,19,19,17))+
  xlab("Site Name")+
  ylab("Average Minimum Temperature (°C)")+
  ylim(28.5,31) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
min_temp_plot

ggsave(min_temp_plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/mintemp_tvesites.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
max_means <- summarySE(daily.mmm.plot, measurevar="max", groupvars=c("logger"))

# plot, minimum temp
max_temp_plot <- ggplot(max_means,aes(x = logger, y = max, color = logger, pch = logger))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = logger, ymax = max+se, ymin = max-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values = c("lightblue","red4","indianred3","mistyrose3"))+
  scale_shape_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values=c(19,19,19,17))+
  xlab("Site Name")+
  ylab("Average Maximum Temperature (°C)")+
  ylim(28.5,31) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
max_temp_plot

ggsave(max_temp_plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/maxtemp_tvesites.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting
mean_means <- summarySE(daily.mmm.plot, measurevar="mean", groupvars=c("logger"))

# plot, minimum temp
mean_temp_plot <- ggplot(mean_means,aes(x = logger, y = mean, color = logger, pch = logger))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = logger, ymax = mean+se, ymin = mean-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values = c("lightblue","red4","indianred3","mistyrose3"))+
  scale_shape_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values=c(19,19,19,17))+
  xlab("Site Name")+
  ylab("Average Mean Temperature (°C)")+
  ylim(28.5,31) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
mean_temp_plot

ggsave(mean_temp_plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/meantemp_tvesites.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
range_means <- summarySE(daily.mmm.plot, measurevar="range", groupvars=c("logger"))

# plot, minimum temp
range_temp_plot <- ggplot(range_means,aes(x = logger, y = range, color = logger, pch = logger))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = logger, ymax = range+se, ymin = range-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values = c("lightblue","red4","indianred3","mistyrose3"))+
  scale_shape_manual(name = "Site",
                     labels = c("CA","CI","PD","SP"),
                     values=c(19,19,19,17))+
  xlab("Site Name")+
  ylab("Average DTV (°C)")+
  #ylim(28.5,31) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
range_temp_plot

ggsave(range_temp_plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/Figures/rangetemp_tvesites.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)

##

master.melt %>%
  ggplot(aes(x = time, y = temp, group = logger,color = logger))+
  geom_smooth()+
  scale_color_manual(values = palsite)+
  xlab('Time of Day')+
  ylab(expression(paste("Temp (",degree,"C)")))+
  scale_x_datetime(
    breaks = seq(as.POSIXct("2019-03-04 00:00:00"),as.POSIXct("2019-03-05 00:00:00 "), "4 hours"), labels = c('00:00', '04:00','08:00','12:00','16:00','20:00','00:00'))+
  facet_wrap(~logger, 1,5)


master.melt %>%
  ggplot(aes(x = time, y = temp, group = logger,color = logger))+
  geom_smooth()+
  scale_color_manual(values = palsite)+
  xlab('Time of Day')+
  ylab(expression(paste("Temp (",degree,"C)")))+
  #scale_x_datetime(
  #  breaks = seq(as.POSIXct("2019-03-04 00:00:00"),as.POSIXct("2019-03-05 00:00:00 "), "6 hours"), labels = c('00:00', '06:00','12:00','18:00','00:00'))+
  #facet_wrap(~logger, 1,5)+
  scale_y_continuous(breaks = seq(28.0,31.0,0.25))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))


# plot all data points as line plot, but only for TVE sites
master.melt.plot = master.melt %>%
  dplyr::filter(logger != "Drago.OR4")

quartz()
ggplot(master.melt.plot, aes(x = datetime, y = temp, color = logger, group = logger))+
  #geom_line()+
  geom_line(aes(group = logger), size = 0.5, alpha = 0.8)+
  scale_color_manual(name = "Site",
                     labels = c("CA","PD","SP","CI"),
                     values = c("lightblue","indianred3","mistyrose3","red4")) +
  ylab("Temperature (°C)")+
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))


# make a figure for Laura with her site color codes
# purple (#6157A2) for Inshore and green (#68B08A) for offshore
# levels of site: Cayo.OR3.arr3 Drago.OR4 Punta.IR1 STRI.IR2 Cristo.IR3.arr1
laurapalsite <- c('#68B08A','#68B08A','#68B08A','#6157A2','#6157A2')


# make a figure of mean temperature by month
head(melt.mmean)
melt.mmean %>%
  ggplot(aes(x = datetime, y = mmean, color = logger, group = logger, shape = logger)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = laurapalsite) +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme_minimal_grid(12) + # the number here sets the font size
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Date", labels = c("2015-07-31 23:45:00" = "2015-07-31",
                                             "2015-08-31 23:45:00" = "2015-08-31",
                                             "2015-09-30 23:45:00" = "2015-09-30",
                                             "2015-10-31 23:45:00" = "2015-10-31",
                                             "2015-11-30 23:45:00" = "2015-11-30",
                                             "2015-12-31 23:45:00" = "2015-12-31",
                                             "2016-01-31 23:45:00" = "2016-01-31",
                                             "2016-02-29 23:45:00" = "2016-02-29",
                                             "2016-03-31 23:45:00" = "2016-03-31",
                                             "2016-04-30 23:45:00" = "2016-04-30",
                                             "2016-05-31 23:45:00" = "2016-05-31",
                                             "2016-06-30 00:00:00" = "2016-06-30")) +
  ylab(expression(paste("Mean Monthly Temp (",degree,"C)")))


# plot mean temperatures and facet by season
head(dailys)
range_sum <- summarySE(data = dailys, measurevar = 'range', groupvars = c('season', 'logger'), na.rm = T)
mean_sum <- summarySE(data = dailys, measurevar = 'mean', groupvars = c('logger', 'season'), na.rm = T)
max_sum <- summarySE(data = dailys, measurevar = 'max', groupvars = c('logger', 'season'), na.rm = T)
min_sum <- summarySE(data = dailys, measurevar = 'min', groupvars = c('logger', 'season'), na.rm = T)

# plot daily range faceted by season
range_sum %>%
  ggplot(aes(x = logger, y = range, color = logger, group = logger, shape = logger))+
  geom_point()+
  geom_errorbar(aes(ymin = range-sd, ymax = range+sd))+
  scale_color_manual(values = laurapalsite)+
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~season,1,4)

mean_sum %>%
  ggplot(aes(x = logger, y = mean, color = logger, group = logger, shape = logger))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd))+
  scale_color_manual(values = laurapalsite)+
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~season,2,2)

max_sum %>%
  ggplot(aes(x = logger, y = max, color = logger, group = logger, shape = logger))+
  geom_point()+
  geom_errorbar(aes(ymin = max-sd, ymax = max+sd))+
  scale_color_manual(values = laurapalsite)+
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~season,2,2)

min_sum %>%
  ggplot(aes(x = logger, y = min, color = logger, group = logger, shape = logger))+
  geom_point()+
  geom_errorbar(aes(ymin = min-sd, ymax = min+sd))+
  scale_color_manual(values = laurapalsite)+
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~season,2,2)



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


#### Compare External Logger Data with Our Data ####

# our dtv data:
daily.mmm <- read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/DailyTempRangeData.csv")
head(daily.mmm)

# external logger combined dtv data:
all.loggers <-read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Field_Hobo_Loggers/data_sheets/AllLoggers.csv")
head(all.loggers)

# re-order and re-name external logger data

col_order <- c("site", "date", "DayMax", "DayMin", "DayMean", "DayRange")
all.loggers <- all.loggers[, col_order]
head(all.loggers)

newnames <- c('logger', 'date', 'max', 'min', 'mean', 'range')
names(all.loggers) <- newnames
head(all.loggers)

#convert date to same format as our dtv dataframe
all.loggers$date <- as.Date(all.loggers$date,"%Y-%m-%d")

#convert dtv dataframe - some data classes seem to be messed up when reading in csv's of data vs. running entire code
daily.mmm$date <- as.Date(daily.mmm$date, "%Y-%m-%d")
daily.mmm$logger <- factor(daily.mmm$logger)

all.loggers$logger <- factor(all.loggers$logger)

str(all.loggers)
str(daily.mmm)

combined.loggers <- rbind(daily.mmm, all.loggers)
head(combined.loggers)

# plot boxplot of daily ranges
combined.loggers %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  ggplot(aes(x=logger, y=range, fill=logger)) +
  geom_boxplot()+
  scale_fill_brewer(palette="RdBu") +
  geom_jitter(shape=16, position=position_jitter(0.2), color="gray", alpha=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6))

# plot ggridges of daily ranges of all data
combined.loggers %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  dplyr::filter(range < 3) %>%
  ggplot(aes(x = range, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_brewer(palette="RdBu") +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  scale_x_continuous(breaks = seq(0,3,0.5))+
  theme_ridges(center = T)

# plot ggridges of daily means of all data
combined.loggers %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  ggplot(aes(x = mean, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_brewer(palette="RdBu") +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  #scale_x_continuous(breaks = seq(0,3,0.5))+
  theme_ridges(center = T)

# plot ggridges of daily max of all data
combined.loggers %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  ggplot(aes(x = max, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_brewer(palette="RdBu") +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  #scale_x_continuous(breaks = seq(0,3,0.5))+
  theme_ridges(center = T)

# plot ggridges of daily mins of all data
combined.loggers %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  ggplot(aes(x = min, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_brewer(palette="RdBu") +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  scale_x_continuous(breaks = seq(0,3,0.5))+
  theme_ridges(center = T)

#### Plot new proposed sites only ####

propsites <- combined.loggers %>%
  dplyr::filter(logger != "Cayo.OR3.arr3") %>%
  dplyr::filter(logger != "Drago.OR4") %>%
  dplyr::filter(logger != "Bella_Vista") %>%
  dplyr::filter(logger != "Cayo_Roldan") %>%
  dplyr::filter(logger != "Cayo_Wilson") %>%
  dplyr::filter(logger != "Punta_Caracol") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Spring") %>%
  dplyr::filter(logger != "Tranquilo_Bay_Summer") %>%
  mutate(logger = factor(logger, levels = c("Cristo.IR3.arr1","Punta.IR1","Hospital_Point","STRI.IR2")))

levels(propsites$logger) <- c("Cristobal Island", "Punta Donato", "Hospital Point", "STRI Point")

# plot daily range as boxplot
colors <- c("#CA0020", "#F4A582", "#0571B0","#92C5DE")

dtv_boxplot <- ggplot(propsites, aes(x=logger, y=range)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.6,
              aes(color = logger)) +
  scale_color_manual(values = colors) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.75,
               aes(fill = logger))+
  scale_fill_manual(values = colors) + # for boxplot
  ylab("Daily Temperature Range (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_classic() +
  theme(legend.position = "none")
dtv_boxplot

ggsave(dtv_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Proposed_Sites/DailyTempRange_colorswap.pdf", width=4.3, height=4, units=c("in"), useDingbats=FALSE)

# plot mean temp with boxplot
mean_boxplot <- ggplot(propsites, aes(x=logger, y=mean)) +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.6,
              aes(color = logger)) +
  scale_color_manual(values = colors) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.75,
               aes(fill = logger))+
  scale_fill_manual(values = colors) + # for boxplot
  ylab("Daily Mean Temperature (°C)") +
  xlab("Site") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7, hjust=.6)) +
  theme_classic() +
  theme(legend.position = "none")
mean_boxplot

ggsave(mean_boxplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Proposed_Sites/DailyMeanTemp_colorswap.pdf", width=4.3, height=4, units=c("in"), useDingbats=FALSE)

# plot daily range as ggridge

dtv_ridgeplot <- ggplot(propsites, aes(x = range, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1)+
  scale_fill_manual(values = colors) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  scale_x_continuous(breaks = seq(0,3,0.5))+
  xlab("Daily Temperature Range (°C)") +
  ylab("Site") +
  theme_ridges(center = T) +
  theme(legend.position = "none")
dtv_ridgeplot

ggsave(dtv_ridgeplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Proposed_Sites/DailyTempRange_ggridges_colorswap.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)

# plot daily mean as ggridge

mean_ridgeplot <- ggplot(propsites, aes(x = mean, y = logger, fill = logger))+
  geom_density_ridges(jittered_points=F, scale = 1.1, rel_min_height = .01, point_shape = "|", point_size = 2, size = 0.1, linetype = "dashed")+
  scale_fill_manual(values = colors) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, scale = 1.1)+
  #scale_x_continuous(breaks = seq(0,3,0.5))+
  xlab("Daily Mean Temperature (°C)") +
  ylab("Site") +
  theme_ridges(center = T) +
  theme(legend.position = "none")
mean_ridgeplot

ggsave(mean_ridgeplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/TemperatureData/Proposed_Sites/DailyMeanTemp_ggridges_colorswap.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)

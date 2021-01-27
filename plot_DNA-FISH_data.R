##################################################
# plot DNA_FISH data created by ImageJ script
##################################################

# define variable: PicturesPerWell
PicturesPerWell = 64

##################################################
# install packages
##################################################
# ipak function: install and load multiple R packages.
#check to see if packages are installed, then install them if needed.
.ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  lapply(pkg, library, character.only = TRUE)
}

# list of packages you want to install/use
packages <- c("ggplot2", "platetools", "dplyr","plotly")

# install and load packages
.ipak(packages)

##################################################
# set working directory and load files
##################################################

# this section requires you to change file/folder paths accordingly!

########################################
# win10 only:
########################################
# set working directory
setwd(choose.dir(default = "", caption = "Set working directory"))
WorkDir = getwd()

# define path to data file
filepathdata <- choose.files(default = "", caption = "Select data file", multi = FALSE)
# make sure that sample or target names do not contain spaces
# skip first two rows of file
data <-  read.table(filepathdata, 
                    quote = "", 
                    sep = ",", 
                    header=TRUE)

# define path to file containing information about plate layout/ information for treatment per well
filepathPlateLayout <- choose.files(default = "", caption = "Select file for PlateLayout", multi = FALSE)
# read in plate layout/ information for treatment per well
PlateLayout <- read.table(filepathPlateLayout,
                          sep = ",",
                          row.names = 1,
                          na.strings = c("", "NA"),
                          header = TRUE)
#rename columns
colnames(PlateLayout) = c(1:12)

########################################
# the usual way:
########################################
# set working directory
#setwd("C:/Users/Manuel/Documents/PhD/lab_book_supplemental/20200615_MT030")

# read in data
#data <-  read.table("ME014_PLA_results.txt", 
#                    quote = "", 
#                    sep = ",", 
#                    header=TRUE, 

# read in plate layout/ information for treatment per well
#PlateLayout <- read.table("Well_Order_ME014_PLA.txt",
#                          sep = ",",
#                          row.names = 1,
#                          na.strings = c("", "NA"),
#                          header = TRUE)

##################################################
# data preparation
##################################################

# remove datasets/rows that are below threshold e.g. c3mean (=IF channel)
#data = filter(data, data$c3mean >600)
#data = filter(data, data$c3mean <1000)

########################################
# write filenumber into added column
########################################
# extract number of original file
#select column containing the filename, then extract last part of the string, which is the number of the original file
filenameNumber <- gsub(".*?_F","",data$File)

# make sure number is numeric
filenameNumber <- as.numeric(filenameNumber)

# Note: be aware that some pictures contained more than one cell, hence they share same file number, but different ROI

# integer division by number of images per well defined at the beginning of the script 
#Modulo operation: find remainder after division of filename number by number of images per well, return and add 1. E.g.: 2/100 gives no remainder, therefore returns 0. Finally add 1.
wellID <- (filenameNumber %/% PicturesPerWell) + 1

# write results to dataScreen list as new column, accordingly. 
data$well_ID <- wellID

########################################
# create data frame: WellData
########################################
# rearange PlateLayout as vector
#listing treatment names on plate from left to right, up to bottom.
WellData = as.vector(t(PlateLayout))

# add well_ID
WellData = as.data.frame(WellData)
WellData$well_ID = 1:96

#rename columns
colnames(WellData) = c("treatment", "well_ID")

# remove rows with empty treatment
WellData = na.omit(WellData)

# add ROI/cell count per well
#count rows per well_ID in data; one row is equal to one ROI/cell; add to WellData
WellData$roi_count <- 0
for(i in unique(data$well_ID)){
  RoiCount <- nrow(data[(data$well_ID == i),])
  WellData[i,]$roi_count <- RoiCount
}

# convert well_ID into ID writen on plate. e.g.: A1,B1..etc
WellData$well_ID = num_to_well(WellData$well_ID)

########################################
# assign treatments, convert well_ID
########################################
# The microscope software is imaging the plate from left to right, row by row. E.g.: A1-A12...
# Hence treatment is assigned in the same direction, since first well imaged is assigned to ID: 1 and so on.

# assign treatment according WellData and write to data as new column
data$treatment = WellData[data$well_ID,"treatment"]

# assign actual well ID to rows in data, according to treatment
data$well_ID = WellData[data$well_ID,"well_ID"]

########################################
# optional: add sub-treatment
########################################
# might be empty or not useful sometimes
data$subtreatment = sub("_.*","", data$treatment)
# same for WellData
WellData$subtreatment = sub("_.*","", WellData$treatment)

##################################################
# construct figures
##################################################
# needed for keeping same x-order in al figures
xform <- list(categoryorder = "array",
              categoryarray = data$treatment)

# convert well_ID into xy coordinates split in two columns
WellData <- mutate(WellData,
                   Row=as.numeric(match(toupper(substr(well_ID, 1, 1)), LETTERS)),
                   Column=as.numeric(substr(well_ID, 2, 5)))

########################################
# heatmap: ROI/cell count per well
########################################
# draw 96 well plate as a heatmap indicating cell count per well
#threshold can be used to identify wells with low cell count more easily
threshold = 0
FigPlateRoiCount <-  plot_ly(data = WellData, x = ~Column, y = ~Row*-1, z = ~roi_count,
                              type = "heatmap",
                              zmin = threshold,
                              zmax = max(WellData$roi_count),
                              colors = c("white","orange"),
                              text = ~roi_count, hoverinfo = "text", hovertext= paste("well:",WellData$well_ID,"<br>treatment:",WellData$treatment,"<br>roi_count:",WellData$roi_count), hoverlabel = list(bgcolor="white")
                      ) %>%
                      layout(#title = paste("Wells containing", threshold,"or less cells appear white"),
                              xaxis = list(ticks="", showline = TRUE,side= "top", showgrid = FALSE, title = "", autotick = FALSE, mirror = "ticks", range = c(0.5,12.5)),
                              yaxis = list(ticks="", tickmode="array", ticktext=LETTERS[1:8],tickvals=c(-1:-8), showline = TRUE, showgrid = FALSE, title = "", autorange = "FALSE", mirror = "ticks", range = c(-8.5,-0.5))
                      )

FigPlateRoiCount

########################################
# barplot: ROI/cell count per well
########################################
FigRoiCountBar <-      plot_ly(data = WellData, x = ~well_ID, y = ~roi_count,
                            type = "bar", 
                            color= ~subtreatment,
                            showlegend = FALSE
) %>%
                          layout(#title = "Channel 2",
                            xaxis = list(tickangle = -45, showgrid = FALSE, title = ""),
                            yaxis = list(#title = "Average ROI intensity", 
                            showline = TRUE)
  )
FigRoiCountBar

########################################
# ROI count per frame per well
########################################

# construct new data frame to hold information
WellFrameInfo = data[!duplicated(data$File),]
WellFrameInfo = WellFrameInfo[,c("File","well_ID","treatment","subtreatment")]
colnames(WellFrameInfo) = c("File","well_ID","treatment","subtreatment")
# View(WellFrameInfo)

# Add ROI count per frame, which is equivalent to filename column. Since One file is one frame.
#count ROIs identified per frame
WellFrameInfo$ROI_count_per_frame <- 0
#set counter for loop
j=1
for(i in unique(data$File)){
  ROICount <- nrow(data[(data$File == i),])
  WellFrameInfo[j,]$ROI_count_per_frame <- ROICount
  j = j+1
}

# plot
#median ROI count per frame per well
FigFrameInfo <-    plot_ly(data = WellFrameInfo, x = ~treatment, y = ~ROI_count_per_frame,
                          type = "violin", jitter = 1, pointpos = 0,
                          marker = list(color = "gray", opacity = 1, size = 3),
                          color= ~subtreatment,
                          showlegend = FALSE, 
                          hoverinfo = "text", hovertext= paste("<br>filename:",WellFrameInfo$File, "<br>ROI count:",WellFrameInfo$ROI_count), hoverlabel = list(bgcolor="white")
) %>%
  layout(#title = "Channel 2",
    xaxis = list(tickangle = -45, showgrid = FALSE, title = ""),
    yaxis = list(#title = "Average ROI intensity", 
      showline = TRUE)
  )
FigFrameInfo

########################################
# scatter: 
########################################
# the bigger the ROI/cell the more foci and what about signal intensity?
FigC <-   plot_ly(data = data, x = ~c2mean, y = ~c2foci,
                  type = "scatter", mode = "markers",
                  color = ~treatment,
                  hoverinfo = "text", hovertext= paste("well:",data$well_ID,"<br>filename:",data$File,"<br>c2mean:",data$c2mean, "<br>ROI number:",data$ROI), hoverlabel = list(bgcolor="white")
          ) %>%
          layout(showlegend = TRUE,
                 yaxis = list(showline = TRUE)
          )

FigC

########################################
# boxplot: median intensity channel 2 per treatment
########################################

# median intensity per ROI per well
FigMeanInt3 <-    plot_ly(data = data, x = ~treatment, y = ~c3mean,
                          type = "box", boxpoints="outlier", jitter = 1, pointpos = 0, notched = TRUE,
                          marker = list(color = "gray", opacity = 1, size = 3),
                          color= ~subtreatment,
                          showlegend = FALSE, 
                          hoverinfo = "text", hovertext= paste("well:",data$well_ID,"<br>filename:",data$File,"<br>c3mean:",data$c3mean, "<br>ROI number:",data$ROI), hoverlabel = list(bgcolor="white")
                  ) %>%
                  layout(#title = "Channel 2",
                          xaxis = list(tickangle = -45, showgrid = FALSE, title = ""),
                          yaxis = list(#title = "Average ROI intensity", 
                          showline = TRUE)
                  )
FigMeanInt3 

########################################
# boxplot: Median Foci count per ROI per treatment
########################################

FigFociCount <-   plot_ly(data = data, x = ~treatment, y = ~c2mean,
                          type = "box", boxpoints="outliers", jitter = 1, pointpos = 0, notched = TRUE,
                          marker = list(color = "gray", opacity = 1, size = 3),
                          color= ~subtreatment,
                          showlegend = FALSE,
                          hoverinfo = "text", hovertext= paste("well:",data$well_ID,"<br>filename:",data$File,"<br>c2foci:",data$c2foci, "<br>ROI number:",data$ROI), hoverlabel = list(bgcolor="white")
                ) %>%
                layout(xaxis = list(tickangle = -45, showgrid = FALSE, title = ""),
                          yaxis = list(#title = "Foci per ROI", 
                          showline = TRUE)
                  )   
FigFociCount

########################################
# Construct final figure
########################################

s1 = subplot(FigMeanInt3, FigFociCount, shareY = FALSE, titleY = TRUE, nrows = 2, shareX = TRUE)
s2 = subplot(FigRoiCountBar,FigFrameInfo, shareY = FALSE, titleY = TRUE, nrows = 2)
s3 = subplot(FigPlateRoiCount, FigC, shareY = FALSE, titleY = TRUE, shareX = FALSE, nrows = 2, titleX = TRUE)


sp = subplot(s1, s2, s3, shareY = FALSE, titleY = TRUE, titleX = TRUE)
sp

# save interactive plot as html in workig directory
#htmlwidgets::saveWidget(sp, file = paste(WorkDir, "/Summary.html",sep=""))

# show PlateLayout
#View(PlateLayout)

##################################################
# save files
##################################################

# save data as csv
#write.table(data, "res.csv", sep=",", col.names = NA)



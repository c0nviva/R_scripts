# visualize foci per frame devided by roi count

# define variable: PicturesPerWell
PicturesPerWell = 57

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
# foci count per frame devided by roi count
########################################

# construct new data frame to hold information
WellFrameInfo = data[!duplicated(data$File),]
WellFrameInfo = WellFrameInfo[,c("File","well_ID","treatment","ndotsframe", "subtreatment")]
colnames(WellFrameInfo) = c("File","well_ID","treatment","ndotsframe", "subtreatment")
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

# divide ndotsframe by ROI_count_per_frame and write into new column
WellFrameInfo$FociDivROI = WellFrameInfo$ndotsframe/WellFrameInfo$ROI_count_per_frame

# plot
#median ROI count per frame per well
FigFrameInfo <-    plot_ly(data = WellFrameInfo, x = ~treatment, y = ~FociDivROI,
                           type = "box", boxpoints="all", jitter = 1, pointpos = 0, notched = TRUE,
                           marker = list(color = "gray", opacity = 1, size = 3),
                           color= ~subtreatment,
                           showlegend = TRUE, 
                           hoverinfo = "text", hovertext= paste("<br>filename:",WellFrameInfo$File, "<br>ROI count:",WellFrameInfo$ROI_count, "<br>Foci count:", WellFrameInfo$ndotsframe), hoverlabel = list(bgcolor="white")
) %>%
  layout(#title = "Channel 2",
    xaxis = list(tickangle = -45, showgrid = FALSE, title = ""),
    yaxis = list(#title = "Average ROI intensity", 
      showline = TRUE)
  )
FigFrameInfo


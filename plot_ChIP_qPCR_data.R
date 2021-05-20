##################################################
# plot ChIP-qPCR data
##################################################
# according to this analysis:
# https://www.thermofisher.com/de/de/home/life-science/epigenetics-noncoding-rna-research/chromatin-remodeling/chromatin-immunoprecipitation-chip/chip-analysis.html

# percentage of input
# e.g if you have taken 10 % as input: enter 10
INfactor = 10

##################################################
# install packages
##################################################
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
.ipak = function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  lapply(pkg, library, character.only = TRUE)
}

# list of packages you want to install/use
packages = c("ggplot2", "dplyr", "platetools", "plotly")

# install and load packages
.ipak(packages)

##################################################
# set working directory and file paths; read in files
##################################################

# this section requires you to change file/folder paths accordingly!

########################################
# win10 only:
########################################
# set working directory
setwd(choose.dir(default = "", caption = "Set working directory"))

# define path to data file
filepath = choose.files(default = "", caption = "Select file for data", multi = FALSE)

# define path to file containing information about plate layout
filepathPlateLayout = choose.files(default = "", caption = "Select file for PlateLayout", multi = FALSE)

##################################################
# read in file and data preparation
##################################################

########################################
# read in files
########################################
# make sure that sample or target names do not contain spaces
# skip first two rows of file
# also, set "Invalid" entries as "NA", which happens when detection limit is hit, e.g. in negative controls
data <-  read.table(filepath, quote = "", sep = "", skip = 2, fill=TRUE, na.strings = "Invalid", header = FALSE)

# extract file name
filename <- gsub(".*?\\\\","",filepath)
filename <- gsub("\\..*","",filename)

# read in PlateLayout file; extra row: first row contains primer pair information
PlateLayout <- read.table(filepathPlateLayout,
                          sep = ",",
                          row.names = 1,
                          na.strings = c("", "NA"),
                          header = TRUE)

########################################
# data preparation/ clean up
########################################
# connect sample name and wellID, write into data and get rid of not needed columns of original data file

# WellData
# rename columns of PlateLayout
colnames(PlateLayout) = c(1:ncol(PlateLayout))
# rearange PlateLayout as vector
#listing treatment names on plate from left to right, top to bottom.
sample = as.vector(t(PlateLayout))
# add well_ID
WellData = as.data.frame(sample)
WellData$well_ID = 1:384
# convert well_ID into ID written on plate. e.g.: A1,B1..etc
WellData$well_ID = num_to_well(WellData$well_ID, plate = 384)
# remove "NA" samples
WellData = na.omit(WellData)

# remove specific row, is necessary sometimes
#  i hate how the qpcr machine exports the results. it is different every time!
# additional outlier filtering
data = filter(data, data[1] == "True")
# data
#remove not needed colums; they usually contain error messages
data = data[,3:6]
# get rid of not cp column
data[3] = NULL
# rename columns of interest in data
colnames(data) = c("well_ID", "sample", "cp")
# write samples into data
data$sample = WellData$sample
# make sure that data$cp is numeric. sometimes this is not the case
data$cp = as.numeric(data$cp)

##################################################
# data analysis
##################################################

########################################
# filters
########################################
# save rawdata
rawdata = data

# exclude data that is to close to detection limit of 40, hence above threshold
data = filter(data, data$cp < 39)
# exclude data with cp = 0
data = filter(data, data$cp > 0)
# remove ddh2o samples
data = data[- grep ("ddh2o", data$sample, ignore.case = TRUE),]

########################################
# extract sample information
########################################
# get unique sample descriptions
# usually technical triplicates per sample description
# name: IP_primer_treatment

# extract IP and input samples
data$identity = gsub("_.*","",data$sample)
# extract targets/primer pairs and write into new column
# remove everything before first "_"
data$primer = sub(".*?_", "", data$sample)
# remove everything after first "_"
data$primer = sub("_.*", "", data$primer)
# extract treatment and write into a new column
data$treatment = gsub(".*_", "", data$sample)

# remove IP or IN and write into new column sample_name
data$sample_name = sub(".*?_","", data$sample)
#data$sample_name = sub(".*?_","", data$sample_name)
#data$sample_name = sub("_.*","", data$sample_name)

########################################
# calc % input
########################################

#additional filtering
# additional outlier filtering
remove = c("J17", "J5", "J1", "J13", "I3", "I18", "K15", "B15", "B17", "G20", "C14", "E4", "E19")
# new methopd easier. just plot cp per well_ID and color code for IN or IP
data = data[!data$well_ID %in% remove,]
# old method using rows.
#data = data[!row.names(data)%in%remove,]

# split data into IP and IN data frames
IP_data = filter(data, data$identity == "IP")
IN_data = filter(data, data$identity == "IN")

# adjusted input
IN_data$cp_adj = IN_data$cp - log(100/INfactor, 2)

# create results data frame
results = IP_data[,c("sample", "sample_name", "primer", "treatment")]

# calculate delta_cp only for sample that are present in IP_data
# sometimes IN data is filtered out
for( i in IP_data$sample_name) {
  IP = IP_data[IP_data$sample_name == i,]
  IN = IN_data[IN_data$sample_name == i,]
  
  # subtract each IP value from each IN value, all combinations
  delta_data = expand.grid(IN$cp_adj,IP$cp)
  # calculate delta_cp
  delta_data[3] = delta_data[1]-delta_data[2]
  colnames(delta_data) = c("IN","IP","delta_cp")
  
  for (j in 1: length(delta_data$delta_cp)){
    
    # extract delta_cp
    delta_cp = delta_data$delta_cp[j]
    
    # calculate percentage_input
    percentage_input = 100 * 2^delta_cp
    
    #write into new column
    delta_data$percentage_input[j] = percentage_input
  }
  # get mean and std
  mean = mean(delta_data$percentage_input)
  sd = sd(delta_data$percentage_input)
  
  # write into results
  results$percentage_input[results$sample_name == i] = mean
  results$sd[results$sample_name == i] = sd
}

##################################################
# plot qPCR data
##################################################

# plot raw cp values IN vs IP
ggplot(data, aes(fill = identity, x = sample_name, y = cp)) +
  geom_bar(position="dodge", stat="Identity") +
  theme_classic() +
  theme(
    axis.text.x= element_text(angle = 45,hjust = 1),
    axis.title.x=element_blank(),
    axis.ticks.x = element_blank()
  )

# plot plate layout
# convert well_ID into xy coordinates split in two columns
# for unfiltered data, set x= rawdata
x = rawdata
x <- mutate(x,
                   Row=as.numeric(match(toupper(substr(well_ID, 1, 1)), LETTERS)),
                   Column=as.numeric(substr(well_ID, 2, 5)))

# draw 96 well plate as a heatmap indicating cell count per well
#threshold can be used to identify wells with low cell count more easily
threshold = 0
FigPlateRoiCount <-  plot_ly(data = x, x = ~Column, y = ~Row*-1, z = ~cp,
                             type = "heatmap",
                             zmin = threshold,
                             zmax = 40,
                             colors = c("white","blue"),
                             text = ~cp, hoverinfo = "text", hovertext= paste("well:",x$well_ID,"<br>sample:",x$sample,"<br>cp:",x$cp), hoverlabel = list(bgcolor="white")
) %>%
  layout(#title = paste("Wells containing", threshold,"or less cells appear white"),
    xaxis = list(ticks="", showline = TRUE,side= "top", showgrid = FALSE, title = "", autotick = FALSE, mirror = "ticks", range = c(0.5,24.5)),
    yaxis = list(ticks="", tickmode="array", ticktext=LETTERS[1:17],tickvals=c(-1:-17), showline = TRUE, showgrid = FALSE, title = "", autorange = "FALSE", mirror = "ticks", range = c(-16.5,-0.5))
  )
FigPlateRoiCount

# plot percentage input
ggplot(results, aes(fill = treatment, x = primer, y = percentage_input)) +
  geom_bar(position="dodge", stat="Identity") +
  geom_errorbar(aes(ymin=percentage_input-sd, ymax = percentage_input+sd), width=.2, position=position_dodge(.9))+
  #  geom_text(aes(label=percentage_input), vjust= -0.5, size = 2) +
  ylim(-0.001,0.5) +
  #  scale_y_continuous(trans = 'log2') +
  #  labs(y = "IP/input") +
  theme_classic() +
  theme(
    axis.text.x= element_text(angle = 45,hjust = 1),
    axis.title.x=element_blank(),
    axis.ticks.x = element_blank()
  )

# save figure in working directory
#ggsave(paste(filename,"",".png",sep=""))
#write.table(results, "res.csv", sep=",", col.names = NA)






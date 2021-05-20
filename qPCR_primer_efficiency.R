##################################################
# primer efficiency curve ChIP-qPCR data
##################################################
# according to : https://toptipbio.com/calculate-primer-efficiencies/

# Dilution factor
DilFactor = 5

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
packages = c("ggplot2", "dplyr", "platetools", "plotly", "stringr")

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

# remove wells / remove outliers
# therefore look at plate layout plot to identify them
# write into list of well identifiers:
remove = c("D11", "A14", "E13", "B18")
# then remove them from the data 
data = data[!data$well_ID %in% remove,]

########################################
# calculate means
########################################
# calculate mean of technical replicates per primer and dilution
data = data %>%
      group_by(sample) %>%
      summarise(cpmean = mean(cp, na.rm = TRUE), cpsd = sd(cp)) 
  

########################################
# extract sample information
########################################
# get primer and dilution information
# usually technical triplicates per sample description
# name: primer_dilution

# split string of data$sample column and write into strv
strv = str_split_fixed(data$sample,"_",2)
# then write primer name into primer column
data$primer = strv[,1]
# and dilution identifier into dilution column
data$dilution = strv[,2]
# make sure dilution is numeric
data$dilution = as.numeric(data$dilution)

# convert dilution identifier into dilution factor, using DilFactor
data$dilution = 1/(DilFactor^data$dilution)

# log10 transformation of dilution column
data$dilloq = log(data$dilution,10)

# remove data above certain dilution step:
data = data[!data$dilloq <= -3,]

# calc slope for each primer pair
primerlist = unique(data$primer)
results = data.frame(primerlist, primerlist, primerlist)
colnames(results)= c("primer", "slope", "eff")

for (p in primerlist) {
  fit = lm(cpmean ~ dilloq, data = data, subset = data$primer == p)
  slope = as.numeric(fit$coefficients[[2]])
  results[results$primer == p,"slope"] = round(slope, 2)
  eff = as.numeric((10^(-1/slope)-1)*100,2)
  eff =  round(eff,2)
  results[results$primer == p,"eff"] = eff
}

results$slope = as.numeric(results$slope)
results$eff = as.numeric(results$eff)

# plot standard curves:
ggplot(data=data, aes(x = dilloq, y = cpmean, group = primer, color = primer))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(min=cpmean-cpsd, max=cpmean+cpsd), width = 0.1)+
  #xlim(-4.5,NA)+
  theme_classic()


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


# plot efficiency value per primer
ggplot(results, aes(primer, eff))+
  geom_bar(stat="identity")+
  geom_text(aes(label = eff),
            vjust = -0.5,
            size = 5) +
  theme_classic()

# write results
#write.table(results, "res.csv", sep=",", col.names = NA)



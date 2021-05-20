##################################################
# co-localization screen analysis
##################################################
# Screen was done by Magdalena Engl
# Data was pre-sorted by Andreas Ettinger
# Channel order: PolII, DNA, EdU
# Data set contains 3x3 plates (triplicates); DNA damage library; EdU positive only; EdU data not present

##################################################
# install packages
##################################################
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
.ipak = function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  lapply(pkg, library, character.only = TRUE)
}

# list of packages you want to install/use
packages = c("ggplot2", "dplyr", "platetools", "plotly", "gridExtra", "scales")

# install and load packages
.ipak(packages)

# furthermore make sure that bioconductor package is installed
# https://www.bioconductor.org/install/
# load bioconductor stuff
# genome wide anotations for human
library('org.Hs.eg.db')

##################################################
# settings
##################################################
# set global theme for ggplot2
theme_set(theme_classic())

##################################################
# setwd and file paths; read in files
##################################################
# set working directory
WorkDir = "C:/Users/Manuel/Documents/PhD/lab_book_supplemental/20210312_ME_coloc_screen"
setwd(WorkDir)

# set path to result file

# magdalena's data:
filepathdataME = "results-with-treatment-and-replicate-nr.csv"

dataME = read.table(filepathdataME,
                    quote = "",
                    sep = ",",
                    header = TRUE)

##################################################
# data preparation
##################################################
# replace Gene ID with Gene symbol
# non ENTREZID will result in NA value
dataME$temp = mapIds(org.Hs.eg.db, keys = dataME$treatment, column = 'SYMBOL', keytype = 'ENTREZID')
# replace NA value with old non ENTREZID entry
dataME$temp[is.na(dataME$temp)] = dataME$treatment[is.na(dataME$temp)]
# replace treatment column with temp column
dataME$treatment = dataME$temp
# finally remove temp column
dataME$temp = NULL

# reorganize data in list, use plateID as sub-category
plateIDnames = unique(dataME$plateID)
# initialize datalist
platelist = list()
# Fill list
for (item in plateIDnames) {
  platedata = filter(dataME, dataME$plateID == item)
  platelist[[item]] = platedata
}

# for compatibility reasons write dataME into data
data = dataME

##################################################
# data exploration
##################################################
# ROI count per plate
# count dataset rows per plate
platedata = data %>%
  group_by(plateID) %>%
  summarise(ROI_count = length(ROI))
# plot ROI count per plate
ggplot(data = platedata, aes(x = plateID, y = ROI_count))+
  geom_bar(stat = "Identity")+
  geom_text(aes(label = ROI_count),
            vjust = -0.5)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank())


# PLot number of ROIs identified per well for each plate
# create plotlist that will hold plots
plotlist = list()
# set counter
i = 1
for (plate in platelist){
  # load data for one plate
  platedata = plate
  # count number of ROIs per well
  platedata = platedata %>%
    group_by(wellID) %>%
    summarise(ROI_count = length(ROI))
  # convert wellID into coordinates and add to platedata
  platedata <- mutate(platedata,
                      Row=as.numeric(match(toupper(substr(wellID, 1, 1)), LETTERS)),
                      Column=as.numeric(substr(wellID, 2, 5)))
  # plot plate as heatmap
plotlist[[i]] =   ggplot(data = platedata, aes(x = Column, y = Row, fill = ROI_count))+
  geom_tile()+
  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8])+
  scale_x_continuous(position = "top", breaks=seq(1, 12))+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  ggtitle(plate$plateID[1])
# counter +1
i = i+1
}
# arrange plots
do.call("grid.arrange", c(plotlist, nrow = 3, ncol=3))

# Intensity histograms per plate
# normalize intensities using z-score
data$zc1mean = scale(data$c1mean, center = TRUE, scale = TRUE)
data$zc2mean = scale(data$c2mean, center = TRUE, scale = TRUE)
# plot histograms for c1mean
hist1 = ggplot(data, aes(x = c1mean))+ 
  geom_histogram(binwidth = 100,
                 color = "darkgreen", 
                 fill = "darkgreen")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0,5000), breaks = seq(0,4000,2000))+
  labs(x="c1mean (Pol II)")+
  theme_bw()
# plot histograms for c2mean
hist2 = ggplot(data, aes(x = c2mean))+ 
  geom_histogram(binwidth = 100,
                 color = "darkred", 
                 fill = "darkred")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0,5000), breaks = seq(0,4000,2000))+
  labs(x="c2mean (SiR-DNA)")+
  theme_bw()
# plot histograms for zc1mean
hist3 = ggplot(data, aes(x = zc1mean))+ 
  geom_histogram(binwidth = 1,
                 color = "darkgreen", 
                 fill = "darkgreen")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(limits = c(NA,5))+
  labs(x="zc1mean (Pol II normalized)")+
  theme_bw()
# plot histograms for zc2mean
hist4 = ggplot(data, aes(x = zc2mean))+ 
  geom_histogram(binwidth = 1,
                 color = "darkred", 
                 fill = "darkred")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(limits = c(NA,5))+
  labs(x="zc2mean (SiR-DNA normalized)")+
  theme_bw()
# arrange plots
grid.arrange(hist1, hist2, hist3, hist4, nrow = 2, ncol = 2)

##################################################
# analysis
##################################################

########################################
# quality control
########################################
# Zhang et al. 1999 (https://doi.org/10.1177/108705719900400206)

# calculate Z'-factor for each plate

# calculate SSMD for each plate

# Plot

########################################
# hit identification
########################################








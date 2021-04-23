##################################################
# co-localization screen analysis
##################################################
# Screen was done by Magdalena Engl
# Data was pre-sorted by Andreas Ettinger
# Channel order: PolII, DNA, EdU

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
packages = c("ggplot2", "dplyr", "platetools", "plotly")

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
setwd("C:/Users/Manuel/Documents/PhD/lab_book_supplemental/20210312_ME_coloc_screen")
# get Wworking directory
WorkDir = getwd()

# set path to results file
# my analysis:
filepathdata = "results.txt"
# magdalena's data:
filepathdataME = "results-with-treatment-and-replicate-nr.csv"

data = read.table(filepathdata,
                  quote = "",
                  sep = ",",
                  header = TRUE)

dataME = read.table(filepathdataME,
                    quote = "",
                    sep = ",",
                    header = TRUE)
##################################################
# data preparation
##################################################

########################################
# optional: filter data
########################################
# remove datasets/rows that are below threshold e.g. c3mean (=IF channel)
# This is needed because results-file is created by python script which does not include this filter, yet
#data = filter(data, data$c2mean >1000)
#data = filter(data, data$c3mean >560)
#data = filter(data, data$c3mean >700)
########################################


# Extract experimentID, plateID, wellID and fieldID from filename
strParts = strsplit(data$File, "-")
# Access first item of each element in list strParts, then truncate string if needed
# then append to data as new column
data$experimentID = gsub("_.*", "", lapply(strParts, `[[`, 1))
data$plateID = gsub(".*_", "", lapply(strParts, `[[`, 1))
data$wellID = lapply(strParts, `[[`, 2)
data$fieldID = lapply(strParts, `[[`, 3)

# reorganize data in list, use plateID as sub-category
plateIDnames = unique(dataME$plateID)
# initialize datalist
datalist = list()
# Fill list
for (item in plateIDnames) {
  dataplate = filter(dataME, dataME$plateID == item)
  datalist[[item]] = dataplate
}

# replace Gene ID with Gene symbol
#dataME$treatment = mapIds(org.Hs.eg.db, dataME$treatment, 'SYMBOL', 'ENTREZID')

##################################################
# visual data exploration
##################################################

# histogram: SiR-DNA channel
ggplot(data = data, aes(x = c2mean)) +
  geom_histogram(
    binwidth = 350,
    aes(y = ..density..),
    fill = "grey",
    color = "black"
  ) +
  #xlim(NA,5000) +
  geom_density()

# interactive scatter: select co-loc coefficients
plot_ly(
  data = data,
  x = ~ pearson_c1vsc3,
  y = ~ spearman_c1vsc3,
  type = "scatter",
  mode = "markers",
  color = ~ M1_manders_c1vsc3,
  hoverinfo = "text",
  hovertext = paste(
    "filename:",
    data$File,
    "<br>plateID:",
    data$plateID,
    "<br>wellID:",
    data$wellID,
    "<br>fieldID:",
    data$fieldID,
    "<br>ROI:",
    data$ROI
  ),
  hoverlabel = list(bgcolor = "white")
) %>%
  layout(showlegend = TRUE,
         yaxis = list(showline = FALSE))

# plot data for one plate
plot_ly(
  data = dataME,
  x = ~ treatment,
  y = ~ pearson_c1vsc3,
  type = "box",
  boxpoints = "outlier",
  jitter = 1,
  pointpos = 0,
  notched = TRUE,
  marker = list(
    color = "gray",
    opacity = 1,
    size = 3
  ),
  color = ~ treatment,
  showlegend = FALSE
)

# plot means
# create dataframe holding means per group
# for data of one plate only use: datalist[[plateIDnames[1]]]
dataMeans = dataME %>%
  group_by(treatment) %>%
  summarise(M1_manders_c1vsc3 = mean(M1_manders_c1vsc3, na.rm = TRUE))
# add column that holds means expressed as z-score
dataMeans$zScore = scale(dataMeans[2])
# reorder data based on z-score
dataMeans$treatment = factor(dataMeans$treatment, levels = dataMeans$treatment[order(dataMeans$zScore)])
# set color identifier
dataMeans$color = cut(
  dataMeans$zScore,
  breaks = c(-Inf, 1, Inf),
  labels = c("black", "red")
)
# plot raw data + means
ggplot(dataME, aes(x = treatment, y = M1_manders_c1vsc3)) +
  geom_jitter(width = 0.25, size = 0.01) +
  geom_bar(data = dataMeans, stat = "identity", alpha = 0.3) +
  theme(axis.text.x = element_text(angle = 90,
                                   
                                   hjust = 1))
# plot means as z-score
ggplot(dataMeans, aes(x = treatment, y = zScore)) +
  geom_point(color = dataMeans$color) +
  geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 5
  ))

# hiytogramm of all coefficients for one sample
untr = dataME[(dataME$treatment == "untr"),]
hist(untr$M1_manders_c1vsc3)

plot_ly(data= dataME,
        x = ~M1_manders_c1vsc3,
        type = "histogram",
        color = ~treatment)

nountr = dataME[(!dataME$treatment == "untr"),]
hist(nountr$M1_manders_c1vsc3)


# Export data
#write.table(dataMeans, "res.csv", sep=",", col.names = NA)

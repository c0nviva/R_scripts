##################################################
# co-localization screen analysis
##################################################
# Screen was done by Magdalena Engl
# Data was pre-sorted by Andreas Ettinger
# Channel order: PolII, EdU, DNA
# Data set contains 3x3 plates (triplicates); DNA damage library; EdU positive only; SiR-DNA data not present
# background substracted dataset

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
# gsea analysis using clusterProfiler
library('clusterProfiler')
# load topGO from bioconductor
# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
library('topGO')

##################################################
# settings
##################################################
# set global theme for ggplot2
theme_set(theme_classic())

# define negative control list
lneg = c("untr","si008","si009","si001","si010")
# define positive control list
lpos = c("ATRi","si002","si003","si004","si011")

##################################################
# functions
##################################################
# function to calculate z-score or robust z-score
# stole it from: https://github.com/HenrikBengtsson/R.basic/blob/master/R/zscore.R
zscore <- function(x, robust = FALSE){
  if (length(x) < 2)
    return(x);
  
  if (robust == TRUE) {
    xavg <- median(x, na.rm = TRUE);
    xdev <- mad(x, center=xavg, na.rm = TRUE);
  } else {
    xavg <- mean(x, na.rm = TRUE);
    xdev <- sd(x, na.rm = TRUE);
  }
  score = (x-xavg)/xdev;
  
  return (score)
}

# function to calculate z'-factor according to Zhang et al. 1999 (https://doi.org/10.1177/108705719900400206)
# for z-factor instead of z'-factor put positive control into x
zfactor <- function (x, neg){
  spc = sd(x, na.rm = TRUE)
  snc = sd(neg,na.rm = TRUE)
  upc = mean(x, na.rm = TRUE)
  unc = mean(neg,na.rm = TRUE)
  factor = 1- (3*spc + 3*snc)/abs(upc-unc)
  
  return (factor)
}

# function to calculate strictly standardized mean difference (SSDM)
# stole it from: https://github.com/Swarchal/phenoDist/blob/master/R/ssmd.R
# a = values
# b = negative control values
calc_ssmd <- function (a, b, verbose = TRUE) 
{
  if (length(a) < 2 | length(b) < 2) {
    stop(call. = FALSE, "Inputs need to be greater at least 2 elements long")
  }
  if (is.numeric(a) == FALSE | is.numeric(b) == FALSE) {
    stop(call. = FALSE, "Input needs to be numeric.")
  }
  
  mu_a <- mean(a, na.rm=TRUE)
  mu_b <- mean(b, na.rm=TRUE)
  var_a <- var(a, na.rm=TRUE)
  var_b <- var(b, na.rm=TRUE)
  
  # if lengths are equal assume correlation and calculate covariance
  if (length(a) == length(b)) {
    cov_ab <- cov(a, b)
    beta <- (mu_a - mu_b)/sqrt(var_a + var_b - 2 * cov_ab)
  } else{ # unequal lengths => not paired , cannot calc covariance
    beta <- (mu_a - mu_b) / sqrt(var_a + var_b)
    if(verbose == TRUE)
    {
      warning("a and b have different lengths. Calculations assumed no correlation.",
              call. = FALSE)
    }
  }
  
  return(beta)
}



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

########################################
# filter
########################################
# remove treatments with less then 100 datapoints available
# currently only works by re-running the script after first run
remove = results$treatment[results$datapoints < 100]
dataME = dataME[!dataME$treatment %in% remove,]

# for compatibility reasons write dataME into data
data = dataME

##################################################
# analysis
##################################################
# normalize intensities per plate, using z-score
for (ID in unique(data$plateID)){
  data$zc1mean[data$plateID == ID] = scale(data$c1mean[data$plateID == ID], center = TRUE, scale = TRUE)
  data$zc2mean[data$plateID == ID] = scale(data$c2mean[data$plateID == ID], center = TRUE, scale = TRUE)
}

# add column to data holding information about control status
data$control = "sample"
data$control[data$treatment %in% lneg] = "neg"
data$control[data$treatment %in% lpos] = "pos"


# reorganize data in list, use plateID as sub-category
plateIDnames = unique(data$plateID)
# initialize datalist
platelist = list()
# Fill list
for (item in plateIDnames) {
  platedata = filter(data, data$plateID == item)
  platelist[[item]] = platedata
}

# calculate means of co-loc coefficients per treatment combining the entire dataset
# M2 does weird things!
results = data %>%
  group_by(treatment) %>%
  summarise(mean_pearson_c1vsc3 = mean(pearson_c1vsc3, na.rm = TRUE),
            mean_spearman_c1vsc3 = mean(spearman_c1vsc3, na.rm = TRUE),
            mean_icq_c1vsc3 = mean(icq_c1vsc3, na.rm = TRUE),
            mean_M1_manders_c1vsc3 = mean(M1_manders_c1vsc3, na.rm = TRUE),
            mean_M2_manders_c3vsc1 = mean(M2_manders_c3vsc1, na.rm = TRUE),
            mean_zc1mean = mean(zc1mean, na.rm = TRUE),
            mean_zc2mean = mean(zc2mean, na.rm = TRUE),
            # count how many data points have been available for all the above calculations
            datapoints = length(treatment))

# calculate z-score for all co-loc coefficients
results$zscore_pearson_c1vsc3 = scale(results$mean_pearson_c1vsc3, center = TRUE, scale = TRUE)
results$zscore_spearman_c1vsc3 = scale(results$mean_spearman_c1vsc3, center = TRUE, scale = TRUE)
results$zscore_icq_c1vsc3 = scale(results$mean_icq_c1vsc3, center = TRUE, scale = TRUE)
results$zscore_M1_manders_c1vsc3 = scale(results$mean_M1_manders_c1vsc3, center = TRUE, scale = TRUE)
results$zscore_M2_manders_c3vsc1 = scale(results$mean_M2_manders_c3vsc1, center = TRUE, scale = TRUE)
results$zscore_zc1mean = scale(results$mean_zc1mean, center = TRUE, scale = TRUE)
results$zscore_zc2mean = scale(results$mean_zc2mean, center = TRUE, scale = TRUE)

# calculate robust z-score for all co-loc coefficients       
results$rz_pearson_c1vsc3 = zscore(results$zscore_pearson_c1vsc3, robust = TRUE)
results$rz_spearman_c1vsc3 = zscore(results$zscore_spearman_c1vsc3, robust = TRUE)
results$rz_icq_c1vsc3 = zscore(results$zscore_icq_c1vsc3, robust = TRUE)
results$rz_M1_manders_c1vsc3 = zscore(results$zscore_M1_manders_c1vsc3, robust = TRUE)
results$rz_M2_manders_c3vsc1 = zscore(results$zscore_M2_manders_c3vsc1, robust = TRUE)
results$rz_zc1mean = zscore(results$mean_zc1mean, robust = TRUE)
results$rz_zc2mean = zscore(results$mean_zc2mean, robust = TRUE)

# add information if control
results$control = "sample"
results$control[results$treatment %in% lneg] = "neg"
results$control[results$treatment %in% lpos] = "pos"

########################################
# data exploration
########################################
# ROI count per plate
# count dataset rows per plate
platedata = data %>%
  group_by(plateID) %>%
  summarise(ROI_count = length(ROI))
# plot ROI count per plate
ggplot(data = platedata, aes(x = plateID, y = ROI_count))+
  geom_bar(stat = "Identity", color = "black")+
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
    scale_fill_gradient(low = "white", high = "darkblue")+
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
                 color = "darkorange", 
                 fill = "darkorange")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0,5000), breaks = seq(0,4000,2000))+
  labs(x="c2mean (EdU)")+
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
                 color = "darkorange", 
                 fill = "darkorange")+
  facet_wrap(~plateID)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3))+
  scale_x_continuous(limits = c(NA,5))+
  labs(x="zc2mean (EdU normalized)")+
  theme_bw()
# arrange plots
grid.arrange(hist1, hist2, hist3, hist4, nrow = 2, ncol = 2)

########################################
# quality control
########################################
# Zhang et al. 1999 (https://doi.org/10.1177/108705719900400206)

listmeasurements = c("pearson_c1vsc3", 
                     "spearman_c1vsc3", 
                     "icq_c1vsc3", 
                     "M1_manders_c1vsc3", 
                     "M2_manders_c3vsc1", 
                     "c1mean", 
                     "c2mean")
# create quality df
quality = data.frame(matrix(nrow = length(platelist), ncol = length(listmeasurements)*3+1))
# counter row
i = 1
for (plate in platelist) {
  platedata = plate
  quality[i,1] = plate$plateID[1]
  # counter column
  j = 2
  for (measurement in listmeasurements){
    dat = platedata[c("treatment",measurement)]
    datneg = dat[dat$treatment %in% lneg,]
    datpos = dat[dat$treatment %in% lpos,]
    # calculate Z'-factor for each plate
    # takes into account pos and neg control
    z_f = zfactor(datpos[,2], datneg[,2])
    # calculate Z-factor for each plate
    # sample vs neg control subdataset
    zf = zfactor(dat[,2], datneg[,2])
    # calculate SSMD for each plate using negative control data
    ssmd = calc_ssmd(dat[,2], datneg[,2])
    # write into quality df
    quality[i,j:(j+2)] = c(z_f,zf,ssmd)
    j = j+3
  }
  i = i+1
}

# Plot
# loop was not working for some reason that i could not figure out :(
titlelist = c("z'-factor", "z-factor", "ssmd")

# set t
t = 1
# get title
title = titlelist[t]
# plot for intensities
plt1 = ggplot(data = quality)+
  # some hardcoded positions that i need to change
  geom_point(aes(x = quality[,1], y = quality[,2+3*5], color = "c1mean (Pol II)"))+
  geom_point(aes(x = quality[,1], y = quality[,2+3*6], color = "c2mean (EdU)"))+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  labs(y=title)+
  ggtitle(title)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x= element_text(angle = 90,
                                  hjust = 1,
                                  size = 7))

  # plot for coefficients
  plt2 = ggplot(data = quality)+
    # some hardcoded positions that i need to change
    geom_point(aes(x = quality[,1], y = quality[,2], color = "pearson"))+
    geom_point(aes(x = quality[,1], y = quality[,2+3], color = "spearman"))+
    geom_point(aes(x = quality[,1], y = quality[,2+3*2], color = "ICQ"))+
    geom_point(aes(x = quality[,1], y = quality[,2+3*3], color = "M1"))+
    geom_point(aes(x = quality[,1], y = quality[,2+3*4], color = "M2"))+
    scale_color_manual(values = c("grey", "red", "blue", "green", "black"))+
    labs(y=title)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x= element_text(angle = 90,
                                    hjust = 1,
                                    size = 7))
  
  # set t
  t = 2
  # get title
  title = titlelist[t]
  # plot for intensities
  plt3 = ggplot(data = quality)+
    # some hardcoded positions that i need to change
    geom_point(aes(x = quality[,1], y = quality[,3+3*5], color = "c1mean (Pol II)"))+
    geom_point(aes(x = quality[,1], y = quality[,3+3*6], color = "c2mean (EdU)"))+
    scale_color_manual(values = c("darkgreen", "darkorange"))+
    labs(y=title)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x= element_text(angle = 90,
                                    hjust = 1,
                                    size = 7))
  
  # plot for coefficients
  plt4 = ggplot(data = quality)+
    # some hardcoded positions that i need to change
    geom_point(aes(x = quality[,1], y = quality[,3], color = "pearson"))+
    geom_point(aes(x = quality[,1], y = quality[,3+3], color = "spearman"))+
    geom_point(aes(x = quality[,1], y = quality[,3+3*2], color = "ICQ"))+
    geom_point(aes(x = quality[,1], y = quality[,3+3*3], color = "M1"))+
    geom_point(aes(x = quality[,1], y = quality[,3+3*4], color = "M2"))+
    scale_color_manual(values = c("grey", "red", "blue", "green", "black"))+
    labs(y=title)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x= element_text(angle = 90,
                                    hjust = 1,
                                    size = 7))
  
  # set t
  t = 3
  # get title
  title = titlelist[t]
  # plot for intensities
  plt5 = ggplot(data = quality)+
    # some hardcoded positions that i need to change
    geom_point(aes(x = quality[,1], y = quality[,4+3*5], color = "c1mean (Pol II)"))+
    geom_point(aes(x = quality[,1], y = quality[,4+3*6], color = "c2mean (EdU)"))+
    scale_color_manual(values = c("darkgreen", "darkorange"))+
    labs(y=title)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x= element_text(angle = 90,
                                    hjust = 1,
                                    size = 7))
  
  # plot for coefficients
  plt6 = ggplot(data = quality)+
    # some hardcoded positions that i need to change
    geom_point(aes(x = quality[,1], y = quality[,4], color = "pearson"))+
    geom_point(aes(x = quality[,1], y = quality[,4+3], color = "spearman"))+
    geom_point(aes(x = quality[,1], y = quality[,4+3*2], color = "ICQ"))+
    geom_point(aes(x = quality[,1], y = quality[,4+3*3], color = "M1"))+
    geom_point(aes(x = quality[,1], y = quality[,4+3*4], color = "M2"))+
    scale_color_manual(values = c("grey", "red", "blue", "green", "black"))+
    labs(y=title)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x= element_text(angle = 90,
                                    hjust = 1,
                                    size = 7))

# loop was not working hence this weird solution emerged
grid.arrange(plt1, plt3, plt5, plt2, plt4, plt6, nrow = 2, ncol = 3)

########################################
# hit identification
########################################

# plot means per treatment (not sure if that is a good representation)
# reorder data based on mean_pearson_c1vsc3
results$treatment = factor(results$treatment, levels = results$treatment[order(results$mean_pearson_c1vsc3)])
ggplot(data = results)+
  geom_point(aes(x = treatment, y = mean_pearson_c1vsc3, color = "pearson"))+
  geom_point(aes(x = treatment, y = mean_spearman_c1vsc3, color = "spearman"))+
  geom_point(aes(x = treatment, y = mean_icq_c1vsc3, color = "ICQ"))+
  geom_point(aes(x = treatment, y = mean_M1_manders_c1vsc3, color = "M1"))+
  geom_point(aes(x = treatment, y = mean_M2_manders_c3vsc1, color = "M2"))+
  scale_color_manual(values = c("grey", "red", "blue", "green", "black"))+
  labs(y = "mean Co-loc coefficient")+
  theme(legend.title = element_blank(),
    axis.text.x= element_text(angle = 90,
                              hjust = 1,
                              size = 7))

# plot normalized PolII int vs normalized EdU int, entire dataset
ggplot(data = results, aes(x = mean_zc1mean, y = mean_zc2mean, color = control))+
  geom_point()+
  scale_color_manual(values = c("black", "orange", "grey"))


# get data
dat = data[,c("treatment", "zc1mean", "control")]
# remove NA values
dat = na.omit(dat)
# reorder data based on mean intensity per treatment
dat$treatment = reorder(dat$treatment, dat$zc1mean, median)
# plot normalized intensity means, ranked
box1 =ggplot(data = dat, aes(x = treatment, y = zc1mean, color = control))+
  geom_boxplot(outlier.shape = NA)+
  ylim(NA,5)+
  scale_color_manual(values = c("black", "orange", "grey"))+
  theme(axis.text.x= element_blank(),
        axis.ticks = element_blank())
  
# get data
dat = data[,c("treatment", "zc2mean", "control")]
# remove NA values
dat = na.omit(dat)
# reorder data based on mean intensity per treatment
dat$treatment = reorder(dat$treatment, dat$zc2mean, median)
# plot normalized intensity means, ranked
box2 =ggplot(data = dat, aes(x = treatment, y = zc2mean, color = control))+
  geom_boxplot(outlier.shape = NA)+
  ylim(NA,5)+
  scale_color_manual(values = c("black", "orange", "grey"))+
  theme(axis.text.x= element_blank(),
        axis.ticks = element_blank())

# plot the two boxplots
grid.arrange(box1, box2, nrow = 1, ncol = 2)


# z-score per co-loc coefficient
# create plotlist that will hold plots
plotlist = list()
# define which coefficients to plot
coefflist = c("zscore_pearson_c1vsc3", "zscore_spearman_c1vsc3", "zscore_icq_c1vsc3", "zscore_M1_manders_c1vsc3", "zscore_M2_manders_c3vsc1", "zscore_zc1mean", "zscore_zc2mean")
for (coeff in coefflist){
  # select dataset
  dataset = results[c("treatment",coeff)]
  # rename columns
  colnames(dataset) = c("treatment", "score")
  
  # reorder data based on z-score
  dataset$treatment = factor(dataset$treatment, levels = dataset$treatment[order(dataset$score)])
  # set color for hits
  dataset$hitcolor = cut(dataset$score,
                         breaks = c(-Inf, 1, Inf),
                         labels = c("black", "red"))
  # write plot into plotlist
  plotlist[[coeff]] = ggplot(data = dataset, aes(x = treatment, y = score))+
    geom_point(color = dataset$hitcolor)+
    geom_hline(yintercept = 1)+
    ggtitle(coeff)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}
# arrange plots
do.call("grid.arrange", c(plotlist, nrow = 2, ncol=4))

# z*-score (robust z-score) 
# create plotlist that will hold plots
plotlist = list()
# define which coefficients to plot
coefflist = c("rz_pearson_c1vsc3", "rz_spearman_c1vsc3", "rz_icq_c1vsc3", "rz_M1_manders_c1vsc3", "rz_M2_manders_c3vsc1", "rz_zc1mean", "rz_zc2mean")
for (coeff in coefflist){
  # select dataset
  dataset = results[c("treatment",coeff)]
  # rename columns
  colnames(dataset) = c("treatment", "score")
  
  # reorder data based on z-score
  dataset$treatment = factor(dataset$treatment, levels = dataset$treatment[order(dataset$score)])
  # set color for hits
  dataset$hitcolor = cut(dataset$score,
                         breaks = c(-Inf, 1, Inf),
                         labels = c("black", "red"))
  # write plot into plotlist
  plotlist[[coeff]] = ggplot(data = dataset, aes(x = treatment, y = score))+
    geom_point(color = dataset$hitcolor)+
    geom_hline(yintercept = 1)+
    ggtitle(coeff)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}
# arrange plots
do.call("grid.arrange", c(plotlist, nrow = 2, ncol=4))


# reorder data based on mean_pearson_c1vsc3
results$treatment = factor(results$treatment, levels = results$treatment[order(results$rz_pearson_c1vsc3)])
ggplot(data = results)+
  geom_point(aes(x = treatment, y = rz_pearson_c1vsc3, color = "pearson"))+
  geom_point(aes(x = treatment, y = rz_spearman_c1vsc3, color = "spearman"))+
  geom_point(aes(x = treatment, y = rz_icq_c1vsc3, color = "ICQ"))+
  geom_point(aes(x = treatment, y = rz_M1_manders_c1vsc3, color = "M1"))+
  geom_point(aes(x = treatment, y = rz_M2_manders_c3vsc1, color = "M2"))+
  geom_hline(yintercept = 1)+
  scale_color_manual(values = c("grey", "red", "blue", "green", "black"))+
  labs(y = "z*-score")+
  #ylim(1,NA)+
  theme(legend.title = element_blank(),
        axis.text.x= element_text(angle = 90,
                                  hjust = 1,
                                  size = 7))


# gene set enrichment anaylsis (GSEA)
# following this tutorial: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# prepare vector using robust z-score for:
original_gene_list = results$rz_spearman_c1vsc3
# name vector
names(original_gene_list) = results$treatment
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 100, # higher number makes it more accurate (default was 10000)
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


# go-term analysis

# create topGOdata object
# ontology:
# BP biological process
# MF molecular function
# CC celular component

# quick function for filtering data
filterthresh <- function(x, thresh =1){
  out = x[x >= thresh]
  return(out)
}

GOdata <- new("topGOdata",
              ontology = "MF",
              allGenes = gene_list,
              geneSel = filterthresh,
              nodeSize = 10,
              mapping = 'org.Hs.eg.db',
              ID = "symbol",
              annot = annFUN.org)











# Export data
#write.table(results, "analysis.csv", sep=",", col.names = NA)


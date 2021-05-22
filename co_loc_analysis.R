##################################################
# co-localization screen analysis
##################################################
# Screen was done by Magdalena Engl
# Data was pre-sorted by Andreas Ettinger
# Channel order: PolII, EdU, DNA
# Data set contains 3x3 plates (triplicates); DNA damage library; EdU positive only; SiR-DNA data not present

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
z_factor <- function (x, neg){
  spc = sd(x)
  snc = sd(neg)
  upc = mean(x)
  unc = mean(neg)
  factor = 1- (3*spc + 3*snc)/abs(upc-unc)
  
  return (factor)
}

# function to calculate strictly standardized mean difference (SSDM)
# stole it from: https://github.com/Swarchal/phenoDist/blob/master/R/ssmd.R
# a = values
# b = negative control values
ssmd <- function (a, b, verbose = TRUE, ...) 
{
  if (length(a) < 2 | length(b) < 2) {
    stop(call. = FALSE, "Inputs need to be greater at least 2 elements long")
  }
  if (is.numeric(a) == FALSE | is.numeric(b) == FALSE) {
    stop(call. = FALSE, "Input needs to be numeric.")
  }
  
  mu_a <- mean(a, ...)
  mu_b <- mean(b, ...)
  var_a <- var(a, ...)
  var_b <- var(b, ...)
  
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
  
  if (verbose == TRUE) {
    ssmd_effect_message(beta)
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
# analysis
##################################################
# normalize intensities per plate, using z-score
for (ID in unique(data$plateID)){
  data$zc1mean[data$plateID == ID] = scale(data$c1mean[data$plateID == ID], center = TRUE, scale = TRUE)
  data$zc2mean[data$plateID == ID] = scale(data$c2mean[data$plateID == ID], center = TRUE, scale = TRUE)
}

# calculate means of co-loc coefficients per treatment combining the entire dataset
results = data %>%
  group_by(treatment) %>%
  summarise(mean_pearson_c1vsc3 = mean(pearson_c1vsc3, na.rm = TRUE),
            mean_spearman_c1vsc3 = mean(spearman_c1vsc3, na.rm = TRUE),
            mean_icq_c1vsc3 = mean(icq_c1vsc3, na.rm = TRUE),
            mean_M1_manders_c1vsc3 = mean(M1_manders_c1vsc3, na.rm = TRUE),
            mean_M2_manders_c3vsc1 = mean(M2_manders_c3vsc1, na.rm = TRUE),
            mean_zc1mean = mean(zc1mean, na.rm = TRUE),
            mean_zc2mean = mean(zc2mean, na.rm = TRUE),
            # calcuate z-score for all co-loc coefficients
            zscore_pearson_c1vsc3 = scale(pearson_c1vsc3, center = TRUE, scale = TRUE),
            zscore_spearman_c1vsc3 = scale(spearman_c1vsc3, center = TRUE, scale = TRUE),
            zscore_icq_c1vsc3 = scale(icq_c1vsc3, center = TRUE, scale = TRUE),
            zscore_M1_manders_c1vsc3 = scale(M1_manders_c1vsc3, center = TRUE, scale = TRUE),
            zscore_M2_manders_c3vsc1 = scale(M2_manders_c3vsc1, center = TRUE, scale = TRUE),
            # calculate robust z-score for all co-loc coefficients
            rz_pearson_c1vsc3 = zscore(pearson_c1vsc3, robust = TRUE),
            rz_spearman_c1vsc3 = zscore(spearman_c1vsc3, robust = TRUE),
            rz_icq_c1vsc3 = zscore(icq_c1vsc3, robust = TRUE),
            rz_M1_manders_c1vsc3 = zscore(M1_manders_c1vsc3, robust = TRUE),
            rz_M2_manders_c3vsc1 = zscore(M2_manders_c3vsc1, robust = TRUE),
            # count how many data points have been available for all the above calculations
            datapoints = length(treatment))

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
                 color = "darkred", 
                 fill = "darkred")+
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
                 color = "darkred", 
                 fill = "darkred")+
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

# calculate Z'-factor for each plate

# calculate SSMD for each plate

# Plot

########################################
# hit identification
########################################

# plot means per treatment (not sure if that is a good representation)
ggplot(data = results)+
  geom_point(aes(x = treatment, y = mean_pearson_c1vsc3), color = "grey")+
  geom_point(aes(x = treatment, y = mean_spearman_c1vsc3), color = "red")+
  geom_point(aes(x = treatment, y = mean_icq_c1vsc3), color = "blue")+
  geom_point(aes(x = treatment, y = mean_M1_manders_c1vsc3), color = "green")+
  geom_point(aes(x = treatment, y = mean_M2_manders_c3vsc1), color = "darkgreen")+
  labs(y = "mean Co-loc coefficient")+
  theme(
    axis.text.x= element_text(angle = 90,
                              hjust = 1,
                              size = 7))

# plot normalized PolII int vs normalized EdU int, entire dataset
ggplot(data = results, aes(x = mean_zc1mean, y = mean_zc2mean))+
  geom_point()


# z-score

ggplot(data = results, aes(x = treatment, y = rz_pearson_c1vsc3))+
  geom_point()+
  geom_hline(yintercept = 1)+
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1,
                               size = 7))

# z*-score (robust z-score) 








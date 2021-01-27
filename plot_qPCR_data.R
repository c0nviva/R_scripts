# plot qPCR data

##################################################
# install packages
##################################################
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
.ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  lapply(pkg, library, character.only = TRUE)
}

# list of packages you want to install/use
packages <- c("ggplot2","scales", "platetools", "dplyr","plotly")

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

# define path to data file
filepath <- choose.files(default = "", caption = "Select file for results_screen", multi = FALSE)

########################################
# the usual way:
########################################
# set working directory
#setwd("C:/Users/Manuel/Documents/PhD/lab_book_supplemental/20200615_MT030")

# read in data
#data <- read.table("ME014_PLA_results.txt", header = TRUE, sep = ",")

# extract file name
filename <- gsub(".*?\\\\","",filepath)
filename <- gsub("\\..*","",filename)

##################################################
# read in file and data preparation
##################################################

########################################
# read in file
########################################
# make sure that sample or target names do not contain spaces
# skip first two rows of file
# also, set "Invalid" entries as "NA"
data <-  read.table(filepath, quote = "", sep = "", skip = 2, fill=TRUE, na.strings = "Invalid")

########################################
# data preparation/ clean up
########################################

# replace "NA" or empty entries and set them to "0"
data[is.na(data)] <- 0

# combine probe name and target name, overwrite column 1
#data[1] <- paste(data$V3, data$V4, sep = "_")

# split data sets:
data$V1 = grepl("1_",data$V3)



# replace oriP with M13, since ther was a mixup in the experiment
# not working
#data <- data.frame(lapply(data, function(x) {sub("oriP", "M13", x)}))

##################################################
# plot qPCR data
##################################################

########################################
# plot data with plotly
########################################
p <- plot_ly(data, x = ~V1, y = ~V10,
        type = "bar", 
        marker = list(color='rgb(288, 0, 58)'),
        error_y = list(array = ~V11, color = '#000000'),
        hoverlabel = list(bgcolor="white")
) %>%
  layout(xaxis = list(tickangle = -45, showgrid = FALSE, title = "", showline=TRUE),
         yaxis = list(title = "", showline = TRUE, type = "log")
  )
# show plot
#p

# save as png in working directory
#orca(p, file = paste(filename,".png",sep=""))

########################################
# plot data with qqplot2
########################################
# set my own custom theme
theme_mt <- function (base_size=7, base_family="arial"){
  theme_classic() %+replace%
    theme(
      text = element_text(face="bold"),
      strip.background = element_rect(fill="#e9e9e9", color="black", size=1),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 15),
      line = element_line(size = 1)
    )
}
# construct plot

ggplot(data, aes(V3, V10)) +
#  geom_hline(yintercept = 1, color="black", size=1)+
  geom_bar(stat="Identity", fill="gray", color="black", size=1) +
  geom_errorbar(aes(y=V10, ymin=V10,ymax=V10+V11), width=0.5, size=1) +
#  scale_y_continuous(trans = "pseudo_log",breaks=10^(0:log10(max(data$V10))), expand = c(0,0)) +
#  coord_trans(y="pseudo_log")+
#  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  labs(y = "Target/Ref") +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0, 0.1))) +
#  annotation_logticks(base=10, sides="l", scaled = TRUE, size = 1)+
  facet_wrap( ~ V4,scales = "free") + 
  theme_mt()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="NONE",
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank()
  )+
  scale_fill_brewer(palette="Spectral")
  
# View data
View(data)

# save figure in working directory
#ggsave(paste(filename,"",".png",sep=""))


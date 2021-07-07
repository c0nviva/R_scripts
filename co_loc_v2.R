##################################################
# Co-localization screen analysis v2
##################################################
# Screen was done by Magdalena Engl
# Images have been pre-processed into tiff files by Andreas Ettinger
# Channel order: PolII, EdU, DNA
# Files have been re-analyzed by my python script:
# co-loc has been done using channel 1-2: PolII vs EdU

##################################################
# Install packages
##################################################
# ipak function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
.ipak = function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  lapply(pkg, library, character.only = TRUE)
}

# List of packages you want to install/use
packages = c("ggplot2", "dplyr", "platetools", "plotly", "gridExtra", "scales")

# Install and load packages
.ipak(packages)

##################################################
# Settings
##################################################
# Set global theme for ggplot2
theme_set(theme_classic())

##################################################
# Functions
##################################################

##################################################
# setwd and file paths; read in files
##################################################

WorkDir = "C:/Users/Manuel/Documents/PhD/lab_book_supplemental/20210312_ME_coloc_screen"
setwd(WorkDir)

# Folder containing csv-files holding plate layouts
LayoutDir = "platelayouts"

# Path to data created by python script
filepathdataMT = "20210607_results.txt"
# Path to file holding indexed files by python script
filepathIndex = "20210607_filepathlist.txt"
# Path to pre-sorted data by Andreas
filepathdataME = "results-with-treatment-and-replicate-nr.csv"

# Read in data created by python script
dataMT = read.table(filepathdataMT,
                    quote = "",
                    sep = ",",
                    header = TRUE)
# Read in indexed files
indexedPaths = scan(file = filepathIndex, what = character(1), sep = ",")
# Read in pre-sorted data
dataME = read.table(filepathdataME,
                    quote = "",
                    sep = ",",
                    header = TRUE)

# Read in treatment information if available
# Get a list of filepaths of available layout files
layoutfilepathlist = list.files(path = LayoutDir, pattern = ".csv", full.names = TRUE)
# Initialize empty layoutList
layoutList = vector(mode = "list", length = 0)
# Load layoutfiles into a list
for (filepath in layoutfilepathlist) {
  # Extract plateID from filename: "..._plateID.csv"
  plateID = sub(".*_(\\S+)\\..*","\\1",filepath)
  # Read in layout file 
  plateLayout = read.table(filepath,
                               sep = ",",
                               row.names = 1,
                               na.strings = c("", "NA"),
                               header = TRUE)
  # Re-arange plateLayout as vector
  # Listing treatment names on plate from left to right, up to bottom.
  plateLayout = as.vector(t(plateLayout))
  # Convert back to dataframe and add well_ID column
  plateLayout = as.data.frame(plateLayout)
  plateLayout$well_ID = 1:96
  # Rename columns
  colnames(plateLayout) = c("treatment", "well_ID")
  # Remove rows containing empty treatment
  plateLayout = na.omit(plateLayout)
  # Convert well_ID into ID as written on plate. e.g.: A1,B1..etc
  plateLayout$well_ID = num_to_well(plateLayout$well_ID)
  # use well_IDs as rownames, need this for treatment asignment, later
  rownames(plateLayout) = plateLayout$well_ID
  # Add to layout list
  layoutList[[plateID]] = plateLayout
}

##################################################
# Data preparation
##################################################
# Extract plate_ID and well_ID from filename and write into new column, if possible
# Hence select only rows that contain filenames containing plateIDs in this regex format: "plate[A-N][1-3]-"
# Then extract plateID using sub
dataMT[grepl("plate[A-N][1-3]-", dataMT$File),"plate_ID"] = sub("-.*","",dataMT[grepl("plate[A-N][1-3]-", dataMT$File),"File"])
# Do the same to extract well_ID
dataMT[grepl("plate[A-N][1-3]-", dataMT$File),"well_ID"] = sub(".*plate[A-N][0-9]-(.\\d+)-.*","\\1",dataMT[grepl("plate[A-N][1-3]-", dataMT$File),"File"])
# Convert A1 into A01 and so on, this is needed for correct treatment asignment 
dataMT$well_ID = gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", dataMT$well_ID, perl = TRUE)

# Add treatment information, if available
dataMT$treatment[grepl("plateA\\d", dataMT$plate_ID)] = 1

for (plateID in names(layoutList)) {
  # Pull plateLayout for plateID
  playout = layoutList[[plateID]]
  # Construct string for data subset selection based on plateID, will include replicates e.g. plateA1, plateA2...
  platestring = paste(plateID, "\\d", sep="")
  # Select datasubset based on plateID, then match well_ID and write treatment into new column
  # Treatment info is pulled from layoutList item with same name as plateID
  # Next line is working, because it matches well_IDs and rownames of playout (had to set rownames for this way earlier!)
  dataMT$treatment[grepl(platestring, dataMT$plate_ID)] = playout[dataMT$well_ID[grepl(platestring, dataMT$plate_ID)], "treatment"]
  
}

# Multiple datasets have been analyzed and written into one file.
# Need to separate them first: DNA damage response (DDR), chromatin regulatory (CR) and ME017
# DDR: Plate A-C; plates in triplicates e.g. A1, A2, A3
# CR: Plate D-N; plates in triplicates, see DDR
# ME017: Filenames start with ME017_

# Get all rows containing data from ME017 experiment based on filename column
#dataME017 = dataMT[grepl("ME017_", dataMT$File),]

# Get CR data
#dataCR = dataMT[grepl("plate[D-N][1-3]-", dataMT$File),]

# Get DDR data
dataDDR = dataMT[grepl("plate[A-C][1-3]-", dataMT$File),]





# Calculate percentage of indexed files that have been processed by the python script
# Get files processed for DDR
fprocessedDDR = unique(dataDDR$File)
# Get files indexed for DDR
findexedDDR = indexedPaths[grepl("plate[A-C][1-3]-", indexedPaths)]
percentageDDR = length(fprocessedDDR)/length(findexedDDR)
print(percentageDDR)





## general setup script to load all packages and install them if the are not already on the system
  #install.packages('easypackages')
  rm(list = ls()) # clear the workspace
  library(easypackages) # then we can do the rest in one go
  # get a list of useful packages
  list.of.packages <- c("Hmisc","ggplot2","caret","gplots","Rmisc","dplyr",
            "MatchIt","optmatch","data.table","plotrix","ggthemes",
            "viridis","coin","plyr","psytabs","RColorBrewer",
            "msir","lmtest", "ggpubr","stats", "reshape2","xtable",
            "ez","apa","parallel", "jmuOutlier","Rtsne","fpc", "cluster",
            "RCurl")

  # check if they are already installed and otherwise install them (not that
  # this doesn't work for biocLite tools!)
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) { install.packages(new.packages)}

  # then load them all
  libraries(list.of.packages)
  rm(list.of.packages, new.packages)

  # set the working directory
  setwd("~")

# note that that the merging step (1a) is run manually as it is somewhat dependent
# on the folder structure and thus might require some editing while running it
# but for the other steps:
# run forest run!

source("./Scripts/variancePart.R")
source("./Scripts/2_localRegression.R")
source("./Scripts/3_Stats.R")
source("./Scripts/4_SymptomCorrelations.R")

# get some basic group comparisons done
Job1 = mcparallel(basicStats(measure = "CT", parcellation = "500aparc"))
Job2 = mcparallel(basicStats(measure = "CT", parcellation = "HCP"))
mccollect(list(Job1, Job2))

# can't do this in parallel cause the underlying function also assigns multi-core processing
varianceStats(measure = "CT", parcellation = "HCP")
varianceStats(measure = "CT", parcellation = "500aparc")

# compute w-scores
Job1 = mcparallel(localRegression(measure = "CT", parcellation = "HCP"))
Job2 = mcparallel(localRegression(measure = "CT", parcellation = "500aparc"))
mccollect(list(Job1, Job2))

# getting Symptom*W-Score correlation
Job1 = mcparallel(symptomCorrelation(measure = "CT", parcellation = "HCP"))
Job2 = mcparallel(symptomCorrelation(measure = "CT", parcellation = "500aparc"))
mccollect(list(Job1, Job2))

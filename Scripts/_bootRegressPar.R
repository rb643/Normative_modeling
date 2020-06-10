bootRegression <- function(measure, parcellation, threshold,...) {
  
  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  
  if (threshold == TRUE){
    load(paste(basedir,"Matched_Age_IQ_CommonPheno_Thresholded.RData",sep=""))  
    bootdir <- paste(basedir,"boot_W_Threshold/",sep="")
  } else if (threshold == FALSE) {
    load(paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))
    bootdir <- paste(basedir,"boot_W/",sep="")
  } 
  
  if (!dir.exists(bootdir))(
    dir.create(bootdir)
  )
  
  # load the two functions to get the SD and LOESS coeffient
  source("./Scripts/_calcSD.R")
  source("./Scripts/_calcLOESS.R")
  source("./Scripts/_normscores.R")
  
  # select the columns that we want to compute z-scores on
  columnnames <- networkDataNames$V1
  
  # order the data by age
  combinedData <- combinedData[order(combinedData$AGE_AT_SCAN),]
  # split the data into equal size age bins so that we can allign the data later and do some subsetting
  agerange <- seq(4.5, 69.5, 1)
  combinedData$agebins <- cut(combinedData$AGE_AT_SCAN, agerange)
  combinedData$numbins <- as.numeric(combinedData$agebins)
  
  combinedData.M <- subset(combinedData, SEX == "Male")
  
  combinedData.M.CTL <- subset(combinedData.M, DX_GROUP == "Control")
  combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
  
  nboot = 1000

library(doParallel)
library(abind)
availableCPU <- detectCores()
cl <- makeCluster(availableCPU-2)
registerDoParallel(cl)
acomb <- function(...) abind(..., along=3)

bt <- foreach(r = 1:nboot, .combine='acomb', .multicombine=TRUE, .packages=c("magrittr","dplyr","msir","Rmisc")) %dopar% {
  #get a random sample
  bootsample = combinedData.M.CTL[sample(nrow(combinedData.M.CTL),dim(combinedData.M.CTL)[1]*1, replace  = TRUE),]
  # compute the scores
  tlt <- normscores(columnnames = columnnames, input_norm = bootsample, input_case = combinedData.M.ASD)
  tlt
}
stopCluster(cl)

# compute the sd
test <- apply(bt, c(1,2), sd, na.rm = TRUE)

save(test, bt, combinedData.M.ASD, file = paste(bootdir,"bootParOutput.RData",sep=""))


}

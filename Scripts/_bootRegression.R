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

  nboot = 100
  
  # create a matrix for the bootstrapped w-scores
  bootmat = array(NA,dim=c(nboot,length(unique(columnnames)),dim(combinedData.M.ASD)[1]))
  
  for (b in 1:nboot){
    #get a random sample
  bootsample = combinedData.M.CTL[sample(nrow(combinedData.M.CTL),dim(combinedData.M.CTL)[1]*0.80, replace  = TRUE),]
  # compute the scores for all males
  j = 1
  print(paste("working on bootsample:",b))
  for (i in unique(columnnames)) {
    #print(paste("working on:",i))
    # get the normative scores for that column
    normativeScores <- calcLOESS(bootsample[,as.character(i)],bootsample$AGE_AT_SCAN)
    
    ctltest <- data.frame(agebins = bootsample$agebins, normativeScores)
    
    # need to check for bins that only have one value as they won't provide a standard deviation
    test <- summarySE(ctltest, measurevar = "Mu", groupvars = "agebins")
    test <- subset(test, N > 2)
    
    normative <- ctltest %>% group_by(agebins) %>% summarise_all(mean)
    colnames(normative) <- c("agebins",paste(i,"_mu",sep = ""),paste(i,"_sd",sep = ""))
    
    keepBins <- match(test$agebins, normative$agebins)
    normative <- normative[keepBins,]
    # in the ASD data compute a z-score for FIQ based on the normative data.
    # first add the normative data to the ASD data
    temp <- merge(combinedData.M.ASD,normative, by = "agebins", all = TRUE)
    duplicates <- match(temp$SUB_ID,combinedData.M.ASD$SUB_ID)
    combinedData.M.ASD <- temp[!is.na(duplicates),]
    
    # then compute the z-score
    name1 <- paste(i,"_mu",sep = "")
    name2 <- paste(i,"_sd",sep = "")
    name3 <- paste(i,"_z",sep = "")
    
    combinedData.M.ASD[,name3] <- ((combinedData.M.ASD[,as.character(i)]-combinedData.M.ASD[,name1])/combinedData.M.ASD[,name2])
    combinedData.M.ASD[,name1] <- NULL
    combinedData.M.ASD[,name2] <- NULL
    
    bootmat[b,j,] <- combinedData.M.ASD[,name3]
    j <- j+1
    
  }
  rm(i, name1, name2, name3, normativeScores, normative, ctltest)
  }  
    
}

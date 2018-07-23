normscores <- function(columnnames, input_norm, input_case, ...) {
  
  j = 1
  outBoot <- array(NA,c(length(unique(columnnames)),dim(input_case)[1]))
  colnames(outBoot) <- input_case$SUB_ID
  rownames(outBoot) <- columnnames
  
  for (i in unique(columnnames)) {
    #print(paste("working on:",i))
    # get the normative scores for that column
    normativeScores <- calcLOESS(input_norm[,as.character(i)],input_norm$AGE_AT_SCAN)
    
    ctltest <- data.frame(agebins = input_norm$agebins, normativeScores)
    
    # need to check for bins that only have one value as they won't provide a standard deviation
    test <- summarySE(ctltest, measurevar = "Mu", groupvars = "agebins")
    test <- subset(test, N > 2)
    
    normative <- ctltest %>% group_by(agebins) %>% summarise_all(mean)
    colnames(normative) <- c("agebins",paste(i,"_mu",sep = ""),paste(i,"_sd",sep = ""))
    
    keepBins <- match(test$agebins, normative$agebins)
    normative <- normative[keepBins,]
    # in the ASD data compute a z-score for FIQ based on the normative data.
    # first add the normative data to the ASD data
    temp <- merge(input_case,normative, by = "agebins", all = TRUE)
    duplicates <- match(temp$SUB_ID,input_case$SUB_ID)
    input_case <- temp[!is.na(duplicates),]
    
    # then compute the z-score
    name1 <- paste(i,"_mu",sep = "")
    name2 <- paste(i,"_sd",sep = "")
    name3 <- paste(i,"_z",sep = "")
    
    input_case[,name3] <- ((input_case[,as.character(i)]-input_case[,name1])/input_case[,name2])
    input_case[,name1] <- NULL
    input_case[,name2] <- NULL
    
    #bootmat[b,j,] <- combinedData.M.ASD[,name3]
    outBoot[j,] <- input_case[,name3]
    j <- j+1
  }
  return(outBoot)
}
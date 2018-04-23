
varianceStats <- function(measure, parcellation, ...) {
#measure <- "MeanCurv"
#parcellation <- "HCP"
  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  vdir <- paste(basedir,"Variance/",sep="")
  load(paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))

  if (!dir.exists(vdir))(
    dir.create(vdir)
  )

library(variancePartition)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

nmeasure <- length(networkDataNames$V1)
vars <- ncol(combinedData)

df <- t(as.matrix(combinedData[,((vars-nmeasure)+1):vars-2]))
meta <- combinedData[,1:((vars-nmeasure)-2)]
form <- ~ AGE_AT_SCAN + (1|DX_GROUP) + FIQ + (1|SEX) + SRS_TOTAL_T + VIQ + HANDEDNESS_SCORES + (1|SET) + (1|SITE_ID)


varPart <- fitExtractVarPartModel(df, form, meta)

vp <- sortCols(varPart)
bars <- plotPercentBars(vp[1:20,])
violins <- plotVarPart(vp)


  pdf(file = paste(vdir,"VarianceContribution.pdf",sep=""), width = 24, height = 12)
    multiplot(bars,violins,cols = 2)
  dev.off()
  
  pdf(file = paste(vdir,"VarianceContribution_s.pdf",sep=""), width = 12, height = 7)
    multiplot(bars,violins,cols = 2)
  dev.off()

}

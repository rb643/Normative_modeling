
basicStats <- function(measure, parcellation, threshold, ...) {
  
basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")

if (threshold == TRUE){
  zdir <- paste(basedir,"W_Threshold/",sep="")
  anovadir <- paste(basedir,"ANOVA_Threshold/",sep="")
  load(paste(basedir,"Matched_Age_IQ_CommonPheno_Thresholded.RData",sep=""))
} else if (threshold == FALSE){
  zdir <- paste(basedir,"W/",sep="")
  anovadir <- paste(basedir,"ANOVA/",sep="")
  load(paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))
}

if (!dir.exists(anovadir))(
  dir.create(anovadir)
)

combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")

# select the columns that we want to compute z-scores on
columnnames <- networkDataNames$V1

# might be interesting to check broad group differences as well
# check for significant differences
df <- melt(combinedData, id.vars=c("DX_GROUP","SEX","SUB_ID","SITE_ID","AGE_AT_SCAN","Euler_Left","Euler_Right"), measure.vars = columnnames)

## Linear Mixed-Effect models
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/mvlombardo/utils/master/cohens_d.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

Fv <- data.frame(Intercept=double(),
                 Dx=double(),
                 Sx=double(),
                 Age=double(),
                 stringsAsFactors=FALSE)
Pv <- Fv
DFv <- Fv
CohensD <- data.frame(CohensD=double())
for (i in unique(df$variable)){
  
  df2 <- subset(df, variable == i)
  m <- lme(fixed=value ~ DX_GROUP+SEX+AGE_AT_SCAN+Euler_Left+Euler_Right, 
           random = ~ 1|SITE_ID, data = df2, 
           control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  a <- anova(m)
  Fv <- rbind(Fv,a$`F-value`)
  Pv <- rbind(Pv,a$`p-value`)
  DFv <- rbind(DFv,a$denDF)
  
  tmp_m = lme(fixed=value ~ SEX+AGE_AT_SCAN+Euler_Left+Euler_Right, 
              random = ~ 1|SITE_ID, data = df2, 
              control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  df2$resid = residuals(tmp_m)
  d = cohens_d(df2$resid[df2$DX_GROUP=="Autism"],df2$resid[df2$DX_GROUP=="Control"])
  CohensD <- rbind(CohensD,d)
}

adjustedP <- p.adjust(Pv[,2], method = "fdr")

MixedModel_Dx <- cbind(Fv[,2], Pv[,2],adjustedP,CohensD)
rownames(MixedModel_Dx) <- unique(df$variable)
colnames(MixedModel_Dx) <- c("F","P","CorrectedP","CohensD")
write.csv(MixedModel_Dx, file = paste(anovadir,"/MixedModelDx.csv", sep=""), row.names = FALSE, col.names = FALSE)

## Linear Mixed-Effect models removing outliers
load(paste(zdir,"CommonPheno_Wscores.RData",sep=""))
columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
regressionData <- combinedData.M.ASD[,as.character(columnnames2)]

df <- melt(combinedData.M, id.vars=c("DX_GROUP","SEX","SUB_ID","SITE_ID","AGE_AT_SCAN","Euler_Left","Euler_Right"), measure.vars = columnnames)
Fv <- data.frame(Intercept=double(),
                 Dx=double(),
                 Age=double(),
                 stringsAsFactors=FALSE)
Pv <- Fv
DFv <- Fv
CohensD <- data.frame(CohensD=double())
for (i in unique(df$variable)){
  
  df2 <- subset(df, variable == i)
  zname <- paste(i,"_z",sep="")
  removeS <- combinedData.M.ASD[abs(combinedData.M.ASD[,zname]) > 2,]$SUB_ID
  df2 <- df2[ ! df2$SUB_ID %in% removeS, ]
  
  m <- lme(fixed=value ~ DX_GROUP+AGE_AT_SCAN+Euler_Left+Euler_Right, 
           random = ~ 1|SITE_ID, data = df2, 
           control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  a <- anova(m)
  Fv <- rbind(Fv,a$`F-value`)
  Pv <- rbind(Pv,a$`p-value`)
  DFv <- rbind(DFv,a$denDF)
  
  tmp_m = lme(fixed=value ~ AGE_AT_SCAN+Euler_Left+Euler_Right, 
              random = ~ 1|SITE_ID, data = df2, 
              control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  df2$resid = residuals(tmp_m)
  d = cohens_d(df2$resid[df2$DX_GROUP=="Autism"],df2$resid[df2$DX_GROUP=="Control"])
  CohensD <- rbind(CohensD,d)
}

adjustedP2 <- p.adjust(Pv[,2], method = "fdr")
sum(adjustedP2 < 0.05)

MixedModel_Dx <- cbind(Fv[,2], Pv[,2],adjustedP2,CohensD)
rownames(MixedModel_Dx) <- unique(df$variable)
colnames(MixedModel_Dx) <- c("F","P","CorrectedP","CohensD")
write.csv(MixedModel_Dx, file = paste(anovadir,"/MixedModelDx_Outliers.csv",sep=""), row.names = FALSE, col.names = FALSE)

## one sample test on w-scores
load(paste(zdir,"CommonPheno_Wscores.RData",sep=""))
columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
regressionData <- combinedData.M.ASD[,as.character(columnnames2)]
regressionData[is.infinite(as.matrix(regressionData))] <- 0 # this is needed cause there is an occasional rounding error of very low numbers...

r <- lm(as.matrix(regressionData)~as.factor(combinedData.M.ASD$SITE_ID))$residuals+colMeans(regressionData)
ps <- matrix(NA,nrow = 308, ncol = 1)
ts <- matrix(NA,nrow = 308, ncol = 1)
e <- matrix(NA,nrow = 308, ncol = 1)
for (i in 1:ncol(r)) {
  
  data <- r[,i]
  data <- data[abs(data)<2]
  result <- t.test(data)
  ps[i,] <- result$p.value
  ts[i,] <- result$statistic
  e[i,] <- mean(data)/sd(data)
}

adjustedP <- p.adjust(ps, method = "fdr")
WModel <- cbind(ps, adjustedP,e)
colnames(WModel) <- c("p","adjusted_p","d")
write.csv(WModel, file = paste(anovadir,"/WModel_Outliers_TTest.csv",sep=""), row.names = FALSE, col.names = FALSE)


## LME test on w-scores
load(paste(zdir,"CommonPheno_Wscores.RData",sep=""))
columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")

df <- melt(combinedData.M.ASD, id.vars=c("DX_GROUP","SEX","SUB_ID","SITE_ID","AGE_AT_SCAN","Euler_Left","Euler_Right"), measure.vars = columnnames2)
Fv <- data.frame(Intercept=double(),
                 Dx=double(),
                 Age=double(),
                 stringsAsFactors=FALSE)
Pv <- Fv
DFv <- Fv
CohensD <- data.frame(CohensD=double())

for (i in unique(df$variable)){
  
  df2 <- subset(df, variable == i)
  df2 <- df2[abs(df2$value) <= 2,]
  
  m <- lme(fixed=value ~ 1+Euler_Left+Euler_Right,
            random = ~ 1|SITE_ID, data = df2, 
            control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  
  a <- anova(m)
  Fv <- rbind(Fv,a$`F-value`)
  Pv <- rbind(Pv,a$`p-value`)
  DFv <- rbind(DFv,a$denDF)
  CohensD <- rbind(CohensD,(mean(df2$value)/sd(df2$value)))
  
}

adjustedP_LME <- p.adjust(Pv[,1], method = "fdr")
WModel_LME <- cbind(Pv[,1], adjustedP_LME,CohensD)
colnames(WModel_LME) <- c("p","adjusted_p","d")
write.csv(WModel_LME, file = paste(anovadir,"/WModel_Outliers_LME.csv",sep=""), row.names = FALSE, col.names = FALSE)

}

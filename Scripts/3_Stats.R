
basicStats <- function(measure, parcellation, ...) {
  
basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
zdir <- paste(basedir,"W/",sep="")

anovadir <- paste(basedir,"ANOVA/",sep="")
if (!dir.exists(anovadir))(
  dir.create(anovadir)
)

load(paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))

# select the columns that we want to compute z-scores on
columnnames <- networkDataNames$V1

# might be interesting to check broad group differences as well
# check for significant differences
df <- melt(combinedData, id.vars=c("DX_GROUP","SEX","SUB_ID","SITE_ID","AGE_AT_SCAN"), measure.vars = columnnames)

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
  m <- lme(fixed=value ~ DX_GROUP+SEX+AGE_AT_SCAN, 
           random = ~ 1|SITE_ID, data = df2, 
           control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  a <- anova(m)
  Fv <- rbind(Fv,a$`F-value`)
  Pv <- rbind(Pv,a$`p-value`)
  DFv <- rbind(DFv,a$denDF)
  
  tmp_m = lme(fixed=value ~ SEX+AGE_AT_SCAN, 
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
load("./Output_500aparc/CT_Age_IQ_Match/W/CommonPheno_Wscores.RData")
columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
regressionData <- combinedData.M.ASD[,as.character(columnnames2)]
df <- melt(combinedData, id.vars=c("DX_GROUP","SEX","SUB_ID","SITE_ID","AGE_AT_SCAN"), measure.vars = columnnames)
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
  zname <- paste(i,"_z",sep="")
  removeS <- combinedData.M.ASD[abs(combinedData.M.ASD[,zname]) > 2,]$SUB_ID
  df2 <- df2[ ! df2$SUB_ID %in% removeS, ]
  
  m <- lme(fixed=value ~ DX_GROUP+SEX+AGE_AT_SCAN, 
           random = ~ 1|SITE_ID, data = df2, 
           control=lmeControl(singular.ok=TRUE,returnObject=TRUE))
  a <- anova(m)
  Fv <- rbind(Fv,a$`F-value`)
  Pv <- rbind(Pv,a$`p-value`)
  DFv <- rbind(DFv,a$denDF)
  
  tmp_m = lme(fixed=value ~ SEX+AGE_AT_SCAN, 
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
regressionData <- combinedData.M.ASD[,as.character(columnnames2)]
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
write.csv(MixedModel_Dx, file = paste(anovadir,"/WModel_Outliers.csv",sep=""), row.names = FALSE, col.names = FALSE)
}

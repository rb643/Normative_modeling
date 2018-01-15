mergeData <- function(measure, parcellation,...) {

basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
if (!dir.exists(basedir))(
  dir.create(basedir)
)

# load both phenotype files
ABIDE1 <- read.csv("./Data/Phenotypic_V1_0b_preprocessed1.csv") # in these files the column names between the two datasets have been manually matched
ABIDE2 <- read.csv("./Data/Phenotypic_V2_0b_preprocessed1.csv")
ABIDE1$SET <- "ABIDE1" # add a site identifier
ABIDE2$SET <- "ABIDE2"

# check common variables
commonColumns <- pmatch(colnames(ABIDE1), colnames(ABIDE2))
commonColumns <- commonColumns[!is.na(commonColumns)]
keepColumns <- colnames(ABIDE2)[commonColumns]
rm(commonColumns)

# only keep the common ones
ABIDE1 <- ABIDE1[,keepColumns]
ABIDE2 <- ABIDE2[,keepColumns]
rm(keepColumns)

# combine them and clear some important datatypes
combinedData <- rbind(ABIDE1,ABIDE2)
combinedData$DX_GROUP <- as.factor(combinedData$DX_GROUP)
combinedData$DX_GROUP <- revalue(combinedData$DX_GROUP, c("1"="Autism", "2"="Control"))
combinedData$SEX <- revalue(as.factor(combinedData$SEX), c("1"="Male", "2"="Female"))
combinedData[combinedData == -9999] <- NA # annoyingly some missing data points were coded as -9999 so recode them
rm(ABIDE1,ABIDE2)

# now load the network data (there is no network data yet so just volume for now!)
#subjects
networkDataSubs1 <- read.csv("./Data/ABIDE-I/Phenotypic.csv", header = TRUE)# this loads the Phenotypic data for only the subjects that were pre-processed
networkDataSubs2 <- read.csv("./Data/ABIDE_II/Phenotype.csv", header = TRUE)
networkDataSubs1 <- as.data.frame(networkDataSubs1$SUB_ID)
networkDataSubs2 <- as.data.frame(networkDataSubs2$SUB_ID)
colnames(networkDataSubs1) <- "SUB_ID"
colnames(networkDataSubs2) <- "SUB_ID"

#data:
if (parcellation == "500aparc") {
  networkData1 <- read.csv(paste("./Data/ABIDE-I/500.aparc/raw",measure,".csv",sep=""), header = FALSE)
  networkData2 <- read.csv(paste("./Data/ABIDE_II/500.aparc/raw",measure,".csv",sep=""), header = FALSE)
  networkDataNames <- read.csv("./Data/columns_id_308.txt", header = FALSE)
} else if (parcellation == "HCP"){
  networkData1 <- read.csv(paste("./Data/ABIDE-I/HCP.fsaverage.aparc/raw",measure,".csv",sep=""), header = FALSE)
  networkData2 <- read.csv(paste("./Data/ABIDE_II/HCP.fsaverage.aparc/raw",measure,".csv",sep=""), header = FALSE)
  networkDataNames <- read.csv("./Data/columns_id_360.txt", header = FALSE)
}

colnames(networkData1) <- networkDataNames$V1
colnames(networkData2) <- networkDataNames$V1

#combine
networkData1 <- cbind(networkDataSubs1,networkData1)
networkData2 <- cbind(networkDataSubs2,networkData2)
rm(networkDataSubs1,networkDataSubs2)
networkData <- rbind(networkData1,networkData2)
rm(networkData1,networkData2)

combinedData <- merge(combinedData, networkData, by = "SUB_ID")
rm(networkData)

# check the variance of all numeric variables as a quick sanity check
nums <- sapply(combinedData, is.numeric)
nums <- combinedData[,nums]
checkVar <- nearZeroVar(nums, freqCut = 75/25, uniqueCut = 50, saveMetrics = FALSE,
                        names = TRUE, foreach = FALSE, allowParallel = TRUE)
rm(nums, checkVar)
# note that this will probably return some variables that should be coded as categorical not numeric

# test for missing data in the potentially important matching variables
sapply(combinedData,function(x) sum(is.na(x)))
# most likely FIQ will have some missing values, so remove all rows with missing values...
combinedData <- combinedData[complete.cases(combinedData[ ,"FIQ"]),]
combinedData <- combinedData[complete.cases(combinedData[ ,300]),] #extra check in case any of the raw values are missing

# check the site distribution
site <- demographicTable(combinedData$DX_GROUP, combinedData$SITE_ID, count = TRUE, percent = TRUE)
site <- site[-c(1),]
write.csv(site,paste(basedir,"site_distribution_PreMatch.csv",sep=""))
rm(site)

# visualize the age distribition
pdf(paste(basedir,"Age_Distribution_PreMatch.pdf",sep=""))
ggplot(combinedData, aes(AGE_AT_SCAN, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("Age") +
  ylab("Density") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme(legend.position = "top")
dev.off()

# visualize the IQ distribition
pdf(paste(basedir,"IQ_Distribution_PreMatch.pdf",sep=""))
ggplot(combinedData, aes(FIQ, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("FIQ") +
  ylab("Density") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme(legend.position = "top")
dev.off()

# # test age differences (using permutation tests)
# independence_test(AGE_AT_SCAN ~ DX_GROUP, data = combinedData)
# # test IQ differences (using permutation tests)
# independence_test(FIQ ~ DX_GROUP, data = combinedData)

save.image(file = paste(basedir,"Raw_CommonPheno.RData",sep=""))

## matching
combinedData$Group <- (combinedData$DX_GROUP == "Autism")# matchit needs a logical
combinedData$Site <- as.numeric(combinedData$SITE_ID)
tempData <- combinedData[,c("SUB_ID","AGE_AT_SCAN","Group","FIQ","Site")] # matchit can't deal with missing data in non-covariates so create a temporary subset
#match.it <- matchit(Group ~ AGE_AT_SCAN + FIQ + Site, data = tempData, method = "genetic") # genetic seems to be the most conservative option
match.it <- matchit(Group ~ AGE_AT_SCAN + FIQ, data = tempData, method = "nearest")
#summary(match.it) # view the difference

# write the matching stats to a csv file
write.csv(summary(match.it)$nn,paste(basedir,"SampleMatchingN.csv",sep=""))
write.csv(summary(match.it)$sum.all,paste(basedir,"SampleStatsPreMatch.csv",sep=""))
write.csv(summary(match.it)$sum.matched,paste(basedir,"SampleStatsPostMatch.csv",sep=""))

save.image(file = paste(basedir,"Raw_CommonPheno.RData",sep=""))

newData <- match.data(match.it) # get new data frame for all matched samples
keepSubs <- match(newData$SUB_ID, combinedData$SUB_ID) # only keep those from the original data
combinedData <- combinedData[keepSubs,]

rm(keepSubs,match.it,newData,tempData)

save.image(file = paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))

### do all checks again after matching ###
# visualize the age distribition
cbPalette2 <- viridis(2)
pdf(paste(basedir,"Age_Distribution_PostMatch.pdf",sep=""))
ggplot(combinedData, aes(AGE_AT_SCAN, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("Age") +
  ylab("Density") +
  scale_fill_manual(values = cbPalette2) +
  theme(
    plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray30"),
    panel.background = element_rect(fill="gray30"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color="white", size=rel(1)),
    axis.title.y  = element_text(color="white",hjust=0.5),
    axis.title.x  = element_text(color="white",hjust=0.5),
    axis.text.y  = element_text(color="white",hjust=1),
    axis.text.x = element_text(hjust=0.95,vjust=0.2),
    legend.text = element_text(color="white", size=rel(1)),
    legend.background = element_rect(fill="gray30"),
    legend.position = "bottom",
    legend.title=element_blank()
  )
dev.off()

# visualize the IQ distribition
pdf(paste(basedir,"IQ_Distribution_PostMatch.pdf",sep=""))
ggplot(combinedData, aes(FIQ, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("FIQ") +
  ylab("Density") +
  scale_fill_manual(values = cbPalette2) +
  theme(
    plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray30"),
    panel.background = element_rect(fill="gray30"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color="white", size=rel(1)),
    axis.title.y  = element_text(color="white",hjust=0.5),
    axis.title.x  = element_text(color="white",hjust=0.5),
    axis.text.y  = element_text(color="white",hjust=1),
    axis.text.x = element_text(hjust=0.95,vjust=0.2),
    legend.text = element_text(color="white", size=rel(1)),
    legend.background = element_rect(fill="gray30"),
    legend.position = "bottom",
    legend.title=element_blank()
  )
dev.off()
rm(cbPalette2)

# check the site distribution
site <- demographicTable(combinedData$DX_GROUP, combinedData$SITE_ID, count = TRUE, percent = TRUE)
site <- site[-c(1),]
write.csv(site,paste(basedir,"/site_distribution_PostMatch.csv",sep=""))
rm(site)
}
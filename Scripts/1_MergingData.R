mergeData <- function(measure, parcellation,...) {

base <- paste("./Output_",parcellation,sep="")
if (!dir.exists(base))(
    dir.create(base)
  )

basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
if (!dir.exists(basedir))(
  dir.create(basedir)
)

# load both phenotype files
ABIDE1 <- read.csv("./Data/Phenotypic_V1_0b_preprocessed1.csv") # in these files the column names between the two datasets have been manually matched
ABIDE2 <- read.csv("./Data/Phenotypic_V2_0b_preprocessed1.csv")
# add them
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

# add the Euler QC metrics
Euler1 <- read.csv("./Data/abide_1_holes.csv", header = FALSE)
Euler2 <- read.csv("./Data/abide_2_holes.csv", header = FALSE)
Euler1 <- cbind(networkDataSubs1,Euler1)
Euler2 <- cbind(networkDataSubs2,Euler2)
colnames(Euler1) <- colnames(Euler2) <- c("SUB_ID","Euler_Left","Euler_Right")
Euler <- rbind(Euler1,Euler2)
rm(Euler1, Euler2)

#combine
networkData1 <- cbind(networkDataSubs1,networkData1)
networkData2 <- cbind(networkDataSubs2,networkData2)
rm(networkDataSubs1,networkDataSubs2)
networkData <- rbind(networkData1,networkData2)
rm(networkData1,networkData2)

combinedData <- merge(combinedData, networkData, by = "SUB_ID")
combinedData <- merge(combinedData, Euler, by = "SUB_ID")
rm(networkData,Euler)

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
p <- ggplot(combinedData, aes(AGE_AT_SCAN, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("Age") +
  ylab("Density") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme(legend.position = "top")

pdf(paste(basedir,"Age_Distribution_PreMatch.pdf",sep=""))
  print(p)
dev.off()

# visualize the IQ distribition
p <- ggplot(combinedData, aes(FIQ, fill=DX_GROUP)) +
  geom_density(alpha=.5) +
  xlab("FIQ") +
  ylab("Density") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme(legend.position = "top")
pdf(paste(basedir,"IQ_Distribution_PreMatch.pdf",sep=""))
  print(p)
dev.off()

# # test age differences (using permutation tests)
# independence_test(AGE_AT_SCAN ~ DX_GROUP, data = combinedData)
# # test IQ differences (using permutation tests)
# independence_test(FIQ ~ DX_GROUP, data = combinedData)

save(networkDataNames,combinedData,file = paste(basedir,"Raw_CommonPheno.RData",sep=""))

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

save(networkDataNames,combinedData,file = paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))

### do all checks again after matching ###
# visualize the age distribition
cbPalette2 <- viridis(2)
p <- ggplot(combinedData, aes(AGE_AT_SCAN, fill=DX_GROUP)) +
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
pdf(paste(basedir,"Age_Distribution_PostMatch.pdf",sep=""))
  print(p)
dev.off()

# visualize the IQ distribition
p <- ggplot(combinedData, aes(FIQ, fill=DX_GROUP)) +
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
pdf(paste(basedir,"IQ_Distribution_PostMatch.pdf",sep=""))
  print(p)
dev.off()

rm(cbPalette2)

# check the site distribution
site <- demographicTable(combinedData$DX_GROUP, combinedData$SITE_ID, count = TRUE, percent = TRUE)
site <- site[-c(1),]
write.csv(site,paste(basedir,"/site_distribution_PostMatch.csv",sep=""))
rm(site)


# check QC
script <- getURL("https://raw.githubusercontent.com/mvlombardo/utils/master/cohens_d.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
script <- getURL("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

x <- subset(combinedData, DX_GROUP == "Autism")$Euler_Left
y <- subset(combinedData, DX_GROUP == "Control")$Euler_Left
t1 <- perm.t.test(x,y,statistic = "mean", B = 20000)
d1 <- cohens_d(x,y, DIM = 1)

P_left <- ggplot(data = combinedData, aes(y = Euler_Left, x = DX_GROUP, fill = DX_GROUP)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6) +
  geom_point(aes(y = Euler_Left, color = DX_GROUP), position = position_jitter(width = .15),
             size = 1, alpha = 0.6) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  ggtitle(paste("MeanDiff = ",round(t1$statistic,3),
          "\n P = ", round(t1$p.value,6),
          "\n Cohens d =", round(d1,3))) +
  ylab("Euler - LH") +
  theme_bw() +
  theme(axis.title.x=element_blank())

x <- subset(combinedData, DX_GROUP == "Autism")$Euler_Right
y <- subset(combinedData, DX_GROUP == "Control")$Euler_Right
t2 <- perm.t.test(x,y,statistic = "t", B = 20000)
d2 <- cohens_d(x,y, DIM = 1)

P_right <- ggplot(data = combinedData, aes(y = Euler_Right, x = DX_GROUP, fill = DX_GROUP)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6) +
  geom_point(aes(y = Euler_Right, color = DX_GROUP), position = position_jitter(width = .15),
             size = 1, alpha = 0.6) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  ggtitle(paste("MeanDiff = ",round(t2$statistic,3),
                "\n P = ", round(t2$p.value,6),
                "\n Cohens d =", round(d2,3))) +
  ylab("Euler - RH") +
  theme_bw() +
  theme(axis.title.x=element_blank())

E_Site_Left <- ggplot(data = combinedData, aes(y = Euler_Left, x = SITE_ID, fill = SITE_ID)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6) +
  geom_point(aes(y = Euler_Left, color = SITE_ID), position = position_jitter(width = .15),
             size = 1, alpha = 0.6) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  ggtitle("Left Hemisphere") + ylab("Euler") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())

E_Site_Right <- ggplot(data = combinedData, aes(y = Euler_Right, x = SITE_ID, fill = SITE_ID)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6) +
  geom_point(aes(y = Euler_Right, color = SITE_ID), position = position_jitter(width = .15),
             size = 1, alpha = 0.6) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  ggtitle("Right Hemisphere") + ylab("Euler") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())

P_Euler <- ggarrange(ggarrange(P_left, P_right, ncol = 2,labels = c("A", "B")),
          E_Site_Left, E_Site_Right, ncol = 1, nrow = 3,labels = c("","C", "D"))

pdf(paste("./Plots/",measure,"Euler_Number.pdf",sep=""))
  print(P_Euler)
dev.off()

## also exclude subjects with high (>300) Euler number
combinedData_Thresholded <- subset(combinedData, Euler_Left <= 300)
combinedData_Thresholded <- subset(combinedData_Thresholded, Euler_Right <= 300)

## matching
combinedData_Thresholded$Group <- (combinedData_Thresholded$DX_GROUP == "Autism")# matchit needs a logical
combinedData_Thresholded$Site <- as.numeric(combinedData_Thresholded$SITE_ID)
tempData <- combinedData_Thresholded[,c("SUB_ID","AGE_AT_SCAN","Group","FIQ","Site")] # matchit can't deal with missing data in non-covariates so create a temporary subset
match.it <- matchit(Group ~ AGE_AT_SCAN + FIQ, data = tempData, method = "nearest")

## check matching on Euler
df <- combinedData[,c("SUB_ID","Group","Euler_Left","Euler_Right")]
zz <- matchit(Group ~ Euler_Left + Euler_Right, data=df, method="genetic",
                         distance="mahalanobis", replace=F)

                         write.csv(summary(zz)$nn,paste(basedir,"GeneMatchN_PostThreshold.csv",sep=""))
                         write.csv(summary(zz)$sum.all,paste(basedir,"GeneMatchStatsPreMatch_PostThreshold.csv",sep=""))
                         write.csv(summary(zz)$sum.matched,paste(basedir,"SGeneMatchStatsPostMatch_PostThreshold.csv",sep=""))

# write the matching stats to a csv file
write.csv(summary(match.it)$nn,paste(basedir,"SampleMatchingN_PostThreshold.csv",sep=""))
write.csv(summary(match.it)$sum.all,paste(basedir,"SampleStatsPreMatch_PostThreshold.csv",sep=""))
write.csv(summary(match.it)$sum.matched,paste(basedir,"SampleStatsPostMatch_PostThreshold.csv",sep=""))

newData <- match.data(match.it) # get new data frame for all matched samples
keepSubs <- match(newData$SUB_ID, combinedData_Thresholded$SUB_ID) # only keep those from the original data
combinedData_Thresholded <- combinedData_Thresholded[keepSubs,]


rm(keepSubs,match.it,newData,tempData)

combinedData <- combinedData_Thresholded # just keep the same name for future steps
save(networkDataNames,combinedData,file = paste(basedir,"Matched_Age_IQ_CommonPheno_Thresholded.RData",sep=""))

# check the site distribution
site_id <- ddply(combinedData,~DX_GROUP+SITE_ID,summarise,count=n())

site_distr <- ggplot(site_id, aes(x = SITE_ID, y = count, fill = DX_GROUP)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  guides(fill=guide_legend(title="Dx"))

pdf(paste(basedir,"Site_Distribution_PostMatch_PostThreshold.pdf",sep=""))
  print(site_distr)
dev.off()

site <- demographicTable(combinedData$DX_GROUP, combinedData$SITE_ID, count = TRUE, percent = TRUE)
site <- site[-c(1),]
write.csv(site,paste(basedir,"/site_distribution_PostMatch_PostThreshold.csv",sep=""))
rm(site, site_distr, site_id)

# create a descriptive table for age and one for IQ
age_table <- ddply(combinedData,~DX_GROUP+SEX,summarise,mean=mean(AGE_AT_SCAN),sd=sd(AGE_AT_SCAN),count=n(),median=median(AGE_AT_SCAN))
write.csv(age_table,paste(basedir,"/Descriptive_Age_PostMatch_PostThreshold.csv",sep=""))
iq_table <- ddply(combinedData,~DX_GROUP+SEX,summarise,mean=mean(FIQ),sd=sd(FIQ),count=n(),median=median(FIQ))
write.csv(iq_table,paste(basedir,"/Descriptive_IQ_PostMatch_PostThreshold.csv",sep=""))

# create a descriptive table for SRS and one for ADOS
srs_table <- ddply(combinedData[(complete.cases(combinedData$SRS_TOTAL_T)),],~DX_GROUP+SEX,summarise,mean=mean(SRS_TOTAL_T),sd=sd(SRS_TOTAL_T),count=n(),median=median(SRS_TOTAL_T))
write.csv(srs_table,paste(basedir,"/Descriptive_SRS_PostMatch_PostThreshold.csv",sep=""))
ados_table <- ddply(combinedData[(complete.cases(combinedData$ADOS_G_TOTAL)),],~DX_GROUP+SEX,summarise,mean=mean(ADOS_G_TOTAL),sd=sd(ADOS_G_TOTAL),count=n(),median=median(ADOS_G_TOTAL))
write.csv(ados_table,paste(basedir,"/Descriptive_ADOS_PostMatch_PostThreshold.csv",sep=""))
}

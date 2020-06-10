localRegression <- function(measure, parcellation, threshold,...) {

  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  
  if (threshold == TRUE){
    load(paste(basedir,"Matched_Age_IQ_CommonPheno_Thresholded.RData",sep=""))  
    zdir <- paste(basedir,"W_Threshold/",sep="")
   } else if (threshold == FALSE) {
    load(paste(basedir,"Matched_Age_IQ_CommonPheno.RData",sep=""))
     zdir <- paste(basedir,"W/",sep="")
   } 

  if (!dir.exists(zdir))(
    dir.create(zdir)
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
combinedData.F <- subset(combinedData, SEX == "Female")

combinedData.M.CTL <- subset(combinedData.M, DX_GROUP == "Control")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")

print("Calculating W-Scores for Male group")
# now in a loop compute the scores for all males
for (i in unique(columnnames)) {
  print(paste("working on:",i))
  # get the normative scores for that columns
  normativeScores <- calcLOESS(combinedData.M.CTL[,i],combinedData.M.CTL$AGE_AT_SCAN)

  ctltest <- data.frame(agebins = combinedData.M.CTL$agebins, normativeScores)

  # need to check for bins that only have one value as they won't provide a standard deviation
  test <- summarySE(ctltest, measurevar = "Mu", groupvars = "agebins")
  test <- subset(test, N > 2)

  normative <- ctltest %>% group_by(agebins) %>% summarise_all(mean)
  colnames(normative) <- c("agebins",paste(i,"_mu",sep = ""),paste(i,"_sd",sep = ""))

  keepBins <- match(test$agebins, normative$agebins)
  normative <- normative[keepBins,]
  # in the ASD data compute a z-score for FIQ based on the normative data.
  # first add the normative data to the ASD data
  combinedData.M.ASD <- merge(combinedData.M.ASD,normative, by = "agebins")
  # then compute the z-score
  name1 <- paste(i,"_mu",sep = "")
  name2 <- paste(i,"_sd",sep = "")
  name3 <- paste(i,"_z",sep = "")

  combinedData.M.ASD[,name3] <- ((combinedData.M.ASD[,i]-combinedData.M.ASD[,name1])/combinedData.M.ASD[,name2])
  combinedData.M.ASD[,name1] <- NULL
  combinedData.M.ASD[,name2] <- NULL

  # add the normative scores to the original as well
  colnames(normativeScores) <- c(name1, name2)
  combinedData.M.CTL[,name1] <- normativeScores[,name1]
  combinedData.M.CTL[,name2] <- normativeScores[,name2]
  combinedData.M.CTL[,name3] <- ((combinedData.M.CTL[,i]-combinedData.M.CTL[,name1])/combinedData.M.CTL[,name2])
}
rm(i, name1, name2, name3, normativeScores, normative, ctltest)

combinedData.F.CTL <- subset(combinedData.F, DX_GROUP == "Control")
combinedData.F.ASD <- subset(combinedData.F, DX_GROUP == "Autism")
print("Calculating W-Scores for Female group")
# now in a loop compute the scores
for (i in unique(columnnames)) {
  print(paste("working on:",i))
  # get the normative scores for that columns
  normativeScores <- calcLOESS(combinedData.F.CTL[,i],combinedData.F.CTL$AGE_AT_SCAN)

  ctltest <- data.frame(agebins = combinedData.F.CTL$agebins, normativeScores)

  # need to check for bins that only have one value as they won't provide a standard deviation
  test <- summarySE(ctltest, measurevar = "Mu", groupvars = "agebins")
  test <- subset(test, N > 2)

  normative <- ctltest %>% group_by(agebins) %>% summarise_all(mean)
  colnames(normative) <- c("agebins",paste(i,"_mu",sep = ""),paste(i,"_sd",sep = ""))

  keepBins <- match(test$agebins, normative$agebins)
  normative <- normative[keepBins,]

  # in the ASD data compute a z-score for FIQ based on the normative data.
  # first add the normative data to the ASD data
  combinedData.F.ASD <- merge(combinedData.F.ASD,normative, by = "agebins")
  # then compute the z-score
  name1 <- paste(i,"_mu",sep = "")
  name2 <- paste(i,"_sd",sep = "")
  name3 <- paste(i,"_z",sep = "")

  combinedData.F.ASD[,name3] <- ((combinedData.F.ASD[,i]-combinedData.F.ASD[,name1])/combinedData.F.ASD[,name2])
  combinedData.F.ASD[,name1] <- NULL
  combinedData.F.ASD[,name2] <- NULL

  # add the normative scores to the original as well
  colnames(normativeScores) <- c(name1, name2)
  combinedData.F.CTL[,name1] <- normativeScores[,name1]
  combinedData.F.CTL[,name2] <- normativeScores[,name2]
  combinedData.F.CTL[,name3] <- ((combinedData.F.CTL[,i]-combinedData.F.CTL[,name1])/combinedData.F.CTL[,name2])
}
rm(i, name1, name2, name3, normativeScores, normative, ctltest)

combinedData <- rbind.fill(combinedData.F.ASD, combinedData.F.CTL, combinedData.M.CTL, combinedData.M.ASD)
rm(combinedData.F.ASD, combinedData.F.CTL, combinedData.M.CTL, combinedData.M.ASD, combinedData.M, combinedData.F)

save(combinedData,networkDataNames,file = paste(zdir,"CommonPheno_Wscores.RData",sep=""))

## now create some plots of the above
regions <- columnnames
# plot all curves per region per group
pcol <- viridis(2)
Vol <- melt(combinedData, id.vars = c("DX_GROUP","AGE_AT_SCAN","SEX"), measure.vars = regions)
colnames(Vol) <- c("Group","Age","Sex","Region","Value")

plots <- list()
for(i in unique(Vol$Region)){
  Voltemp <- subset(Vol, Region == i)
  print(paste("working on:",i))

plots[[i]] <- ggplot(Voltemp, aes(x=Age, y=Value, fill = Group)) + # make basic plot
  geom_point(aes(colour = Group), size = 3, alpha = 0.7) + # add scatterplot
  geom_smooth(method=loess,span = 0.4,   # Add local regression line
              se=TRUE) +  # Add shaded confidence region
              scale_color_manual(values=pcol) +
              scale_fill_manual(values=pcol) +
  theme_dark() +
  ggtitle(paste("Region:",i)) +
  theme(text = element_text(size = 18),
        plot.background = element_rect(fill="gray30"),
        panel.background = element_rect(fill="gray30"),
        axis.text.y = element_text(color="white", size=rel(1)),
        axis.text.x = element_text(color="white", size=rel(1)),
        axis.title.x = element_text(color="white", size=rel(1)),
        axis.title.y = element_text(color="white", size=rel(1)),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(color="white", size=rel(0.7)),
        legend.background = element_rect(fill="gray30"),
        legend.position = "bottom",
        plot.title = element_text(color="white",hjust = 0.5)) +
  facet_wrap(~Sex)
}

## optional printing of all curves for each brain region
## this takes really long so optional!
# print("Printing local regression plots")
# pdf(paste(zdir,"SmoothedCurves.pdf",sep=""), onefile = TRUE, paper = "a4r")
# for (i in seq(length(plots))) {
#   print(plots[[i]])
# }
# dev.off()
# rm(pcol,plots,i,Vol,Voltemp)

# plot z-scores for autism group
zregions <- sapply(columnnames, function(x) paste(x, "_z", sep=""))
ASD <- subset(combinedData, DX_GROUP == "Autism")
zVol <- melt(ASD, id.vars = c("SEX","AGE_AT_SCAN"), measure.vars = zregions)
colnames(zVol) <- c("Sex","Age","Region","ZScore")
pcol <- viridis(1)

plots <- list()
for(i in unique(zVol$Region)){
  zVoltemp <- subset(zVol, Region == i)
  print(paste("working on:",i))

plots[[i]] <- ggplot(zVoltemp, aes(x=Age, y=ZScore)) + # make basic plot
  geom_point(colour = "white", fill = pcol,
             stroke = 1.5,pch = 21, size = 3,shape = 16, alpha = 0.6) + # add scatterplot
  geom_smooth(method=loess,span = 0.4,   # Add local regression line
              se=TRUE) +  # Add shaded confidence region
              theme_dark() +
              ggtitle(paste("Region:",i)) +
              theme(text = element_text(size = 18),
                    plot.background = element_rect(fill="gray30"),
                    panel.background = element_rect(fill="gray30"),
                    axis.text.y = element_text(color="white", size=rel(1)),
                    axis.text.x = element_text(color="white", size=rel(1)),
                    axis.title.x = element_text(color="white", size=rel(1)),
                    axis.title.y = element_text(color="white", size=rel(1)),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.x = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_text(color="white", size=rel(0.7)),
                    legend.background = element_rect(fill="gray30"),
                    legend.position = "bottom",
                    plot.title = element_text(color="white",hjust = 0.5)) +
  facet_wrap(~Sex, scales = "free")
}

## optional printing of all curves for each brain region
## this takes really long so optional!
# print("Printing z-score plots")
# pdf(paste(zdir,"ASD_Wscores.pdf",sep=""), onefile = TRUE, paper = "a4r")
# for (i in seq(length(plots))) {
#   print(plots[[i]])
# }
# dev.off()
# rm(zregions, ASD,zVol,plots,i,zVoltemp, pcol)
}

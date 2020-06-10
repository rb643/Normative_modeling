computeRatios <- function(measure, parcellation, threshold, ...) {
  #measure <- "MeanCurv"
  #parcellation <- "HCP"
  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  
  if (threshold == TRUE){
    rdir <- paste(basedir,"Ratio_Threshold/",sep="")
    load(paste(basedir,"W_Threshold/CommonPheno_Wscores.RData",sep=""))
  } else if (threshold == FALSE){
    rdir <- paste(basedir,"Ratio/",sep="")
    load(paste(basedir,"W/CommonPheno_Wscores.RData",sep=""))
  }
  
  if (!dir.exists(rdir))(
    dir.create(rdir)
  )

columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
regressionData <- combinedData.M.ASD[,as.character(columnnames2)]

## Ratios
ratio <- matrix(NA,nrow = nrow(regressionData),ncol = 1)
for (i in 1:nrow(regressionData)){
  data <- abs(regressionData[i,])
  x <- length(data[data>2])
  y <- length(data[data<2])
  ratio[i,] <- x/y
}
colnames(ratio) <- "ABS_Ratio"

ratioPos <- matrix(NA,nrow = nrow(regressionData),ncol = 1)
for (i in 1:nrow(regressionData)){
  data <- regressionData[i,]
  data <- data[data > 0]
  x <- length(data[data>2])
  y <- length(data[data<2])
  ratioPos[i,] <- x/y
}
colnames(ratioPos) <- "ABS_Ratio_pos"

ratioNeg <- matrix(NA,nrow = nrow(regressionData),ncol = 1)
for (i in 1:nrow(regressionData)){
  data <- regressionData[i,]
  data <- data[data < 0]
  x <- length(data[data< -2])
  y <- length(data[data> -2])
  ratioNeg[i,] <- x/y
}
colnames(ratioNeg) <- "ABS_Ratio_neg"

combinedData.M.ASD <- cbind(combinedData.M.ASD,ratio, ratioPos, ratioNeg)

pcol <- viridis(1)
ratioPos <- ggplot(combinedData.M.ASD, aes(x=agebins, y=ABS_Ratio_pos)) + # make basic plot
  geom_boxplot(color = pcol,alpha = 0.5) +
  geom_point(color = pcol, size = 3, alpha = 0.7) + # add scatterplot
  geom_smooth(method=loess,span = 0.4,   # Add local regression line
              se=TRUE) +  # Add shaded confidence region
  scale_color_manual(values=pcol) +
  scale_fill_manual(values=pcol) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  ylab("Positive Ratio") + xlab("Age") +
  theme(
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)
  ) +
  coord_flip()

ratioNeg <- ggplot(combinedData.M.ASD, aes(x=agebins, y=ABS_Ratio_neg)) + # make basic plot
  geom_boxplot(color = pcol,alpha = 0.5) +
  geom_point(color = pcol, size = 3, alpha = 0.7) + # add scatterplot
  geom_smooth(method=loess,span = 0.4,   # Add local regression line
              se=TRUE) +  # Add shaded confidence region
  scale_color_manual(values=pcol) +
  scale_fill_manual(values=pcol) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  ylab("Negative Ratio") + xlab("Age") +
  theme(
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)
  ) +
  coord_flip()

ratioABS <- ggplot(combinedData.M.ASD, aes(x=agebins, y=ABS_Ratio)) + # make basic plot
  geom_boxplot(color = pcol,alpha = 0.5) +
  geom_point(color = pcol, size = 3, alpha = 0.7) + # add scatterplot
  geom_smooth(method=loess,span = 0.4,   # Add local regression line
              se=TRUE) +  # Add shaded confidence region
  scale_color_manual(values=pcol) +
  scale_fill_manual(values=pcol) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  ylab("Ratio") + xlab("Age") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)
  )

Ratio <- ggarrange(ratioABS,
                  ggarrange(ratioNeg, ratioPos, ncol = 2,labels = c("B", "C")),
                     ncol = 1, nrow = 2,labels = c("A",""))

pdf(paste(rdir,"Ratios.pdf",sep=""))
  print(Ratio)
dev.off()

## stratify for each region
ratioSpatial <- matrix(NA,ncol = ncol(regressionData),nrow = 1)
for (i in 1:ncol(regressionData)){
  data <- abs(regressionData[,i])
  x <- length(data[data>2])
  y <- length(data[data<2])
  ratioSpatial[,i] <- x/y
}

dataf <- as.data.frame(t(ratioSpatial))
prev <- ggplot(dataf,aes(V1, fill = pcol)) + 
  geom_density(alpha = 0.8) + theme_minimal() +
  xlab("Percentage of subjects with atypical w-score") + ylab("number of \n brain regions") +
  ggtitle(paste0("Median prevalence: ", round(median(dataf$V1),3))) +
  theme(legend.position = 'none') 

pdf(paste(rdir,"Prevalence.pdf",sep=""), width = 10, height = 2.5)
  print(prev)
dev.off()

rs <- cbind((ratio), combinedData.M.ASD$SUB_ID)
write.csv(rs,paste(rdir,"Prevalence.csv"),sep = "", col.names = FALSE)


rs2 <- cbind((ratioSpatial), combinedData.M.ASD$SUB_ID)
write.csv(rs2,paste(rdir,"Region_Prevalence.csv"),sep = "", col.names = FALSE)

}


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
df2 <- df
colnames(df) <- c("Group","Sex","ID","Site","Age","Region","value")
#main anova
print("Running ANOVA")
aov1 <- ezANOVA(data = df, dv = value, wid = ID, between = .(Group,Sex), between_covariates = .(Site,Age), within = Region, detailed = TRUE, type = 3)
write.csv(xtable(aov1$ANOVA),paste(anovadir,"GroupMeanANOVA_",measure,"_table.csv",sep=""))
write.csv(xtable(aov1$`Sphericity Corrections`),paste(anovadir,"GroupMeanANOVA_SperCorr_",measure,"_table.csv",sep=""))
aov1.out <- xtable(anova_apa(aov1, sph_corr = "gg"))
write.csv(aov1.out,paste(anovadir,"GroupMeanANOVA_",measure,"_APA.csv",sep=""))
print(aov1)
rm(aov1, aov1.out)

print("Running Pairwise Tests")
#region-wise wilcoxon rank test with FDR correction
result <- compare_means(value~DX_GROUP, data = df2, method = "wilcox.test", paired = FALSE,
              group.by = c("variable","SEX"), p.adjust.method = "BH")
result <- subset(result, p.adj < 0.05)

# plot all significant differences as boxplots and write out tables of p-values
for (i in unique(result$SEX)) {
tempDF <- subset(result, SEX == i)
write.csv(tempDF,paste(anovadir,"GroupMeanTests_",i,"_FDR.csv",sep=""))

tempDF2 <- subset(combinedData, SEX == i)
keepcols <- match(colnames(tempDF2),tempDF$variable)
keepcols <- !is.na(keepcols)
keepcols <- c("DX_GROUP","SEX",colnames(tempDF2)[keepcols])
tempDF2 <- tempDF2[,keepcols]

tempDF2 <- melt(tempDF2,id.vars=c("DX_GROUP","SEX"), measure.vars = tempDF$variable)
p <- ggboxplot(tempDF2, x="DX_GROUP", y="value",notch = TRUE,
               add = "jitter", color = "DX_GROUP") +
               facet_wrap(~variable, scales = "free")

pdf(file = paste(anovadir,"GroupMeanTests_",i,"_FDR.pdf",sep=""), width = 30, height = 30)
print(p)
dev.off()
}
rm(tempDF2,tempDF,keepcols,result,p,df,df2,i)


# but also visualize the overall scores in a more accescible way
SData <- as.data.frame(list())
for (i in unique(columnnames)) {
  print(paste("working on:",i))
  temp <- summarySE(combinedData, measurevar = i, na.rm = TRUE, groupvars = c("DX_GROUP","SEX"))
  temp$point <- i
  colnames(temp) <- c("DX_GROUP","SEX","N","Mean","sd","se","ci","region")
  SData <- rbind(SData,temp)
}
dfl1 <- melt(SData, id.vars=c("DX_GROUP","SEX","region","ci"), measure.vars = c("Mean"))
pdf(file = paste(anovadir,"GroupMeans.pdf",sep=""), height = 50, width = 15)
ggplot(dfl1, aes(x=region, y=value, ymin=(value-ci), ymax=(value + ci), color = DX_GROUP)) +
  geom_pointrange(alpha = 0.6, shape = 21, size = 1.5) +
  # geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Mean (CI)") + ggtitle(measure) +
  theme(axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        legend.background = element_rect(fill='white'),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~SEX)
dev.off()
rm(dfl1, SData, temp,i)

}

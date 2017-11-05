symptomCorrelation <- function(measure, parcellation,...) {

  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  outdir <- paste(basedir,"SymptomCorrelation/",sep="")
  load(paste(basedir,"W/","CommonPheno_Wscores.RData",sep=""))

  if (!dir.exists(outdir))(
    dir.create(outdir)
  )

# select main severity variables
#columnnames <- networkDataNames$V1
columnnames2 <- as.factor(paste(networkDataNames$V1,"_z",sep=""))
#variables <- c("AGE_AT_SCAN","FIQ","AQ_TOTAL", "SCQ_TOTAL", "SRS_TOTAL_T","ADOS_G_TOTAL", "ADOS_2_TOTAL",as.character(columnnames))
variables2<- c("AGE_AT_SCAN","FIQ","AQ_TOTAL", "SCQ_TOTAL", "SRS_TOTAL_T","ADOS_G_TOTAL", "ADOS_2_TOTAL",as.character(columnnames2))

combinedData.M <- subset(combinedData, SEX == "Male")
combinedData.F <- subset(combinedData, SEX == "Female")
combinedData.M.CTL <- subset(combinedData.M, DX_GROUP == "Control")
combinedData.M.ASD <- subset(combinedData.M, DX_GROUP == "Autism")
combinedData.F.CTL <- subset(combinedData.F, DX_GROUP == "Control")
combinedData.F.ASD <- subset(combinedData.F, DX_GROUP == "Autism")

# compute pairwise correlations
correlationData <- combinedData.M.ASD[,variables2]

# check for remaining NA's
check <- !is.na(sapply(correlationData,function(x) standarderror((x))))
correlationData <- correlationData[,check]

correlationMatrix <- rcorr(as.matrix(correlationData), type = "spearman")

r <- as.data.frame(correlationMatrix$r)
p <- as.data.frame(correlationMatrix$P)
n <- as.data.frame(correlationMatrix$n)

nn <- ncol(p)
nvars <- 7

p <- p[1:nvars,(nvars+1):nn]
r <- r[1:nvars,(nvars+1):nn]
n <- n[1:nvars,(nvars+1):nn]

# and... presto!
pnorm <- melt(as.matrix(p))
rnorm <- melt(as.matrix(r))
nnorm <- melt(as.matrix(n))

pplot <- ggplot(pnorm, aes(Var1, Var2)) +
geom_tile(aes(fill = value)) +
scale_fill_viridis(option='B') +
ggtitle("Correlation P-Value") +
theme(
  strip.text = element_text(colour = "white",face="bold", size=9,lineheight=5.0),
  strip.background = element_rect(fill="black", colour="white",size=1),
  plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
  plot.background = element_rect(fill="gray30"),
  panel.background = element_rect(fill="gray30"),
  panel.border = element_rect(fill=NA,color="gray30", size=0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.2),
  panel.grid.minor = element_line(size = 0.2),
  axis.text = element_text(color="white", size=rel(1)),
  axis.title.y  = element_blank(),
  axis.title.x  = element_blank(),
  axis.text.y  = element_text(color="white",hjust=1),
  axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
  legend.text = element_text(color="white", size=rel(0.7)),
  legend.background = element_rect(fill="gray30"),
  legend.position = "bottom",
  legend.title=element_blank()
) +
coord_flip()

rplot <- ggplot(rnorm, aes(Var1, Var2)) +
geom_tile(aes(fill = value)) +
scale_fill_viridis(option='B') +
ggtitle("Correlation R-Value") +
theme(
  strip.text = element_text(colour = "white",face="bold", size=9,lineheight=5.0),
  strip.background = element_rect(fill="black", colour="white",size=1),
  plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
  plot.background = element_rect(fill="gray30"),
  panel.background = element_rect(fill="gray30"),
  panel.border = element_rect(fill=NA,color="gray30", size=0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.2),
  panel.grid.minor = element_line(size = 0.2),
  axis.text = element_text(color="white", size=rel(1)),
  axis.title.y  = element_blank(),
  axis.title.x  = element_blank(),
  axis.text.y  = element_text(color="white",hjust=1),
  axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
  legend.text = element_text(color="white", size=rel(0.7)),
  legend.background = element_rect(fill="gray30"),
  legend.position = "bottom",
  legend.title=element_blank()
) +
coord_flip()

nplot <- ggplot(nnorm, aes(Var1, Var2)) +
geom_tile(aes(fill = value)) +
scale_fill_viridis(option='B') +
ggtitle("Correlation N") +
theme(
  strip.text = element_text(colour = "white",face="bold", size=9,lineheight=5.0),
  strip.background = element_rect(fill="black", colour="white",size=1),
  plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
  plot.background = element_rect(fill="gray30"),
  panel.background = element_rect(fill="gray30"),
  panel.border = element_rect(fill=NA,color="gray30", size=0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.2),
  panel.grid.minor = element_line(size = 0.2),
  axis.text = element_text(color="white", size=rel(1)),
  axis.title.y  = element_blank(),
  axis.title.x  = element_blank(),
  axis.text.y  = element_text(color="white",hjust=1),
  axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
  legend.text = element_text(color="white", size=rel(0.7)),
  legend.background = element_rect(fill="gray30"),
  legend.position = "bottom",
  legend.title=element_blank()
) +
coord_flip()

pdf(file = paste(outdir,"SymptomCorrelations.pdf",sep=""), height = 20, width = 40)
multiplot(pplot, rplot, nplot, ncol = 1)
dev.off()


p.adjusted <- matrix(p.adjust(as.numeric(data.matrix(p)),method="fdr"),nrow=nrow(p),ncol=ncol(p))
colnames(p.adjusted) <- colnames(p)
rownames(p.adjusted) <- rownames(p)

p_fdr <- melt(p.adjusted)
p_fdr <- merge(p_fdr,rnorm, by = c("Var1","Var2"))
p_fdr <- subset(p_fdr, value.x < 0.05)

adjustedR <- merge(p_fdr,rnorm, by = c("Var1","Var2"), all = TRUE)
adjustedR[is.na(adjustedR[,"value.x"]),"value"] <- 0

fdrplot1 <- ggplot(p_fdr, aes(Var1, Var2)) +
geom_tile(aes(fill = value.x)) +
scale_fill_viridis(option='B') +
ggtitle("Correlation FDR Corrected P-Value") +
theme(
  strip.text = element_text(colour = "white",face="bold", size=9,lineheight=5.0),
  strip.background = element_rect(fill="black", colour="white",size=1),
  plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
  plot.background = element_rect(fill="gray30"),
  panel.background = element_rect(fill="gray30"),
  panel.border = element_rect(fill=NA,color="gray30", size=0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.2),
  panel.grid.minor = element_line(size = 0.2),
  axis.text = element_text(color="white", size=rel(1)),
  axis.title.y  = element_blank(),
  axis.title.x  = element_blank(),
  axis.text.y  = element_text(color="white",hjust=1),
  axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
  legend.text = element_text(color="white", size=rel(0.7)),
  legend.background = element_rect(fill="gray30"),
  legend.position = "bottom",
  legend.title=element_blank()
) +
coord_flip()

fdrplot2 <- ggplot(p_fdr, aes(Var1, Var2)) +
geom_tile(aes(fill = value.y)) +
scale_fill_viridis(option='B') +
ggtitle("Correlation R-Value") +
theme(
  strip.text = element_text(colour = "white",face="bold", size=9,lineheight=5.0),
  strip.background = element_rect(fill="black", colour="white",size=1),
  plot.title = element_text(color="white",hjust=0.5,vjust=1, size=rel(2)),
  plot.background = element_rect(fill="gray30"),
  panel.background = element_rect(fill="gray30"),
  panel.border = element_rect(fill=NA,color="gray30", size=0.5, linetype="solid"),
  panel.grid.major = element_line(size = 0.2),
  panel.grid.minor = element_line(size = 0.2),
  axis.text = element_text(color="white", size=rel(1)),
  axis.title.y  = element_blank(),
  axis.title.x  = element_blank(),
  axis.text.y  = element_text(color="white",hjust=1),
  axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
  legend.text = element_text(color="white", size=rel(0.7)),
  legend.background = element_rect(fill="gray30"),
  legend.position = "bottom",
  legend.title=element_blank()
) +
coord_flip()

pdf(file = paste(outdir,"SymptomCorrelations_FDR.pdf",sep=""), height = 20, width = 40)
multiplot(fdrplot1, fdrplot2, ncol = 1)
dev.off()

save(p_fdr,pnorm,rnorm,nnorm,adjustedR,file = paste(outdir,"CorrelationScores.RData",sep=""))
adjustedR1 <- subset(adjustedR, Var1 == "SRS_TOTAL_T")
write.table(adjustedR1$Var2, file = paste(outdir,"FDR_RScores_Label.txt",sep=""), row.names = FALSE, col.names = FALSE)
write.table(adjustedR1$value, file = paste(outdir,"FDR_RScores_SRS.txt",sep=""), row.names = FALSE, col.names = FALSE)

adjustedR2 <- subset(adjustedR, Var1 == "ADOS_G_TOTAL")
write.table(adjustedR2$Var2, file = paste(outdir,"FDR_RScores_Label_ADOSG.txt",sep=""), row.names = FALSE, col.names = FALSE)
write.table(adjustedR2$value, file = paste(outdir,"FDR_RScores_ADOSG.txt",sep=""), row.names = FALSE, col.names = FALSE)


}

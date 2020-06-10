tsneParcellation <- function(measure, parcellation, threshold, ...) {

  basedir <- paste("./Output_",parcellation,"/",measure,"_Age_IQ_Match/",sep="")
  
  if (threshold == TRUE){
    zdir <- paste(basedir,"W_Threshold/",sep="")
    outdir <- paste(basedir,"Subgrouping_Threshold/",sep="")
    load(paste(zdir,"CommonPheno_Wscores.RData",sep=""))
  } else if (threshold == FALSE){
    zdir <- paste(basedir,"W/",sep="")
    outdir <- paste(basedir,"Subgrouping/",sep="")
    load(paste(zdir,"CommonPheno_Wscores.RData",sep=""))
  }
  
  if (!dir.exists(outdir))(
    dir.create(outdir)
  )
  

# grab the data
  x <- combinedData
  measures <- networkDataNames$V1
  z <- x[,as.character(measures)]

# run tsne to reduce to 2 dimensions
  sm.tsne <- Rtsne(as.matrix(z), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
# extract the distance matrix from the tsne
  t.dist <- as.matrix(dist(sm.tsne$Y))
# run partitioning around medoids with silhouette estimation to get the number of optimal clusters
  pamk.best <- pamk(t.dist)
# run PAM with that number of clusters
  pam.res <- pam(t.dist, pamk.best$nc)
# put the cluster in a separate variable
  groups <- as.data.frame(pamk.best$pamobject$clustering)
  groups <- as.data.frame(pam.res$clustering)

# optional plotting
  # pam.colors <- rainbow(pamk.best$nc)
  # #point.colors <- sapply(pam.res$clustering, function(i) { pam.colors[i] })
  # point.colors <- sapply(x$DX_GROUP, function(i) { pam.colors[i] })

# create an empty plot for grid arrange
  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") +
    theme(plot.background = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks = element_blank())

# create 2D density plots
test <- cbind(as.data.frame(sm.tsne$Y),(groups))
colnames(test) <- c("V1","V2","Cl")
test$Cl <- as.factor(test$Cl)
dens1 <- ggplot(test,aes(x=V1,y=V2))+
        stat_density2d(aes(fill = Cl, alpha = ..level..), geom="polygon") +
  theme_void() +
  theme(legend.position = "none") +
    scale_colour_continuous(guide = FALSE)

test <- cbind(as.data.frame(sm.tsne$Y),(x$DX_GROUP))
colnames(test) <- c("V1","V2","Dx")
test$Dx <- as.factor(test$Dx)
dens2 <- ggplot(test,aes(x=V1,y=V2))+
    stat_density2d(aes(fill = Dx, alpha = ..level..), geom="polygon") +
  theme_void() +
  theme(legend.position = "none") +
scale_colour_continuous(guide = FALSE)

# create the plots for the left panel with colours based on the PAM output
  test <- cbind(as.data.frame(sm.tsne$Y),(groups))
  colnames(test) <- c("V1","V2","Cl")
  test$Cl <- as.factor(test$Cl)
  # scatterplot
  l1 <- ggplot(test,aes(V1,V2)) + geom_point(aes(color = Cl)) +
    theme_void() +
    theme(legend.position = "none")

  # density plot for dim 1
  t1 <- ggplot(test, aes(V1, fill = Cl)) + geom_density(alpha = 0.5) +
    theme_void() +
    theme(legend.position = "bottom",legend.justification = 'center',
                          legend.background = element_rect(fill="white",
                                                           size=0.5, linetype="solid", 
                                                           colour ="black"))
  # density plot for dim 2
  l2 <- ggplot(test, aes(V2, fill = Cl)) + geom_density(alpha = 0.5) +
    theme_void() +
    theme(legend.position = "none") + 
  coord_flip()

# create the plots for the right panel with colours based on diagnosis
  test <- cbind(as.data.frame(sm.tsne$Y),x$DX_GROUP)
  colnames(test) <- c("V1","V2","Dx")
  # scatterplot
  l3 <- ggplot(test,aes(V1,V2)) + geom_point(aes(color = Dx)) +
    theme_void() +
    theme(legend.position = "none")
  # density plot for dim 1
  t2 <- ggplot(test, aes(V1, fill = Dx)) + geom_density(alpha = 0.5) +
    theme_void() +
    theme(legend.position = "bottom", legend.justification = 'center',
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black"))
  # density plot for dim 2
  l4 <- ggplot(test, aes(V2, fill = Dx)) + geom_density(alpha = 0.5) +
    coord_flip()  +
    theme_void() +
    theme(legend.position = "none")
  

tsnePlot <- ggarrange(ggarrange(t1, dens1, l1, l2, ncol = 2, nrow = 2, labels = c(""), widths = c(3,0.9), heights = c(1.5,3)),
                      ggarrange(t2, dens2, l3, l4, ncol = 2, nrow = 2,  labels = c(""), widths = c(3,0.9), heights = c(1.5,3)),
                      ncol = 2, nrow = 1,labels = c("A","B"))  
  
# put it all into one plot
savename <- paste(outdir,"PamClusters.pdf", sep = "")
pdf(savename, width = 12, height = 4)
  print(tsnePlot)
dev.off()


  }

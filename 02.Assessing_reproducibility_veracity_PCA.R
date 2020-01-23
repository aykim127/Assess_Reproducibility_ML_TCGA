# load required packages
library(TCGA2STAT)
library(stats) # to use prcomp function
library(factoextra) # to produce 2D PCA plot

##########################################################################################
# 0. Load data RNASeq RPKM data 
##########################################################################################

# blca0 <- getTCGA(disease="BLCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
# brca0 <- getTCGA(disease="BRCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
# coad0 <- getTCGA(disease="COAD", data.type="RNASeq", type="RPKM", clinical=TRUE)    # no recurrent or normal cases in the dataset
# coadread0 <- getTCGA(disease="COADREAD", data.type="RNASeq", type="RPKM", clinical=TRUE)  # no recurrent or normal cases in the dataset
# esca0 <- getTCGA(disease="ESCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
# hnsc0 <- getTCGA(disease="HNSC", data.type="RNASeq", type="RPKM", clinical=TRUE)
# kipan0 <- getTCGA(disease="KIPAN", data.type="RNASeq", type="RPKM", clinical=TRUE)
# kirc0 <- getTCGA(disease="KIRC", data.type="RNASeq", type="RPKM", clinical=TRUE)
# kirp0 <- getTCGA(disease="KIRP", data.type="RNASeq", type="RPKM", clinical=TRUE)  # no recurrent or normal cases in the dataset
# laml0 <- getTCGA(disease="LAML", data.type="RNASeq", type="RPKM", clinical=TRUE)  # empty dataset 
# lihc0 <- getTCGA(disease="LIHC", data.type="RNASeq", type="RPKM", clinical=TRUE)
# luad0 <- getTCGA(disease="LUAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
# lusc0 <- getTCGA(disease="LUSC", data.type="RNASeq", type="RPKM", clinical=TRUE)
# ov0 <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE)      # no normal cases in the dataset
# read0 <- getTCGA(disease="READ", data.type="RNASeq", type="RPKM", clinical=TRUE)  # no normal cases in the dataset
# stad0 <- getTCGA(disease="STAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
# thca0 <- getTCGA(disease="THCA", data.type="RNASeq", type="RPKM", clinical=TRUE)  # no normal cases in the dataset
# ucec0 <- getTCGA(disease="UCEC", data.type="RNASeq", type="RPKM", clinical=TRUE) 
# 
# save.image("/extra/akim127/Project1/TCGA_all.RData")

## Load the data files obtained from the TCGA2STAT library
load("C:/Users/akim127/Google Drive/Personal/TCGA_all.RData")


data_all <- list(blca0, brca0, esca0, hnsc0, kipan0, kirc0, lihc0, luad0, lusc0, stad0, ucec0)
cancer_name <- c("BLCA", "BRCA", "ESCA", "HNSC", "KIPAN", "KIRC", "LIHC", "LUAD", "LUSC", "STAD", "UCEC")

## Create a data frame that store mean and standard deviations of RPKM reads for the datasets listed above
stat_all <- data.frame(cancer_name, N_genes=rep(0, length(data_all)), 
                          N_tumor=rep(0, length(data_all)), 
                          tumor_mean=rep(0, length(data_all)), tumor_sd=rep(0, length(data_all)),
                          N_normal=rep(0, length(data_all)), 
                          normal_mean=rep(0, length(data_all)), normal_sd=rep(0, length(data_all)))
## Store the PCA results
pca_fit_all <- list()


for (i in 1:length(data_all)){
  d0 <- data_all[[i]]
  d1 <- SampleSplit(d0$dat)
  stat_all[i, 3] <- dim(d1$primary.tumor)[2]
  stat_all[i, 6] <- dim(d1$normal)[2]
  
  tumor <- t(d1$primary.tumor)
  normal <- t(d1$normal)
  
  stat_all[i, 4] <- mean(tumor)
  stat_all[i, 5] <- sd(tumor)
  stat_all[i, 7] <- mean(normal)
  stat_all[i, 8] <- sd(normal)
  
  ## Create a label column
  tumor_group <- rep("tumor", nrow(tumor))
  normal_group <- rep("normal", nrow(normal))
  
  tumor <- data.frame(tumor, tumor_group)
  normal <- data.frame(normal, normal_group)
  
  colnames(tumor)[ncol(tumor)] <- "group"
  colnames(normal)[ncol(tumor)] <- "group"
  
  ## Combine the tumor and normal datasets 
  d <- rbind(tumor, normal)
  
  ## Remove columns with all zero values
  d <- d[, colSums(d != 0) > 10]
  stat_all[i, 2] <- dim(d)[2]
  
  pca.fit <- prcomp(d[,-ncol(d)], center = TRUE, scale = TRUE)
  pca_fit_all[[i]] <- pca.fit
  
  ## Need to change the file path
  mypath <- file.path("C:","Users","akim127", "Google Drive", "Personal", "For_revision",
                      paste("scree_", cancer_name[i], ".png", sep = ""))
  png(file=mypath)
  screeplot(pca.fit, type = "l", npcs = min(length(pca.fit$sdev), 500), 
            main = paste("Screeplot of the first", min(length(pca.fit$sdev), 500), "for", cancer_name[i]), pch=19)
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=1)
  dev.off()

  prop.var <- pca.fit$sdev^2 / sum(pca.fit$sdev^2)
  cumpro <- cumsum(pca.fit$sdev^2 / sum(pca.fit$sdev^2))
  cumpro_index <- min(which(cumpro > 0.9))
  
  ## Need to change the file path
  mypath2 <- file.path("C:","Users","akim127", "Google Drive", "Personal", "For_revision",
                      paste("cumpro_", cancer_name[i], ".png", sep = ""))
  png(file=mypath2)
  plot(cumpro[0:min(length(pca.fit$sdev), 500)], xlab = "PC #", ylab = "Amount of explained variance", 
       main = paste("Cumulative variance plot for", cancer_name[i]), pch=19)
  
  abline(v = cumpro_index, col="blue", lty=5)
  abline(h = cumpro[cumpro_index], col="blue", lty=5)
  legend("bottomright", legend=paste("Cut-off @", "PC", cumpro_index),
         col=c("blue"), lty=5, cex=1)
  dev.off()
  
  ## 2D PCA plot (1st and 2nd PC) from all the features (genes) in the dataset
  ## Need to change the file path
  mypath3 <- file.path("C:","Users","akim127", "Google Drive", "Personal", "For_revision",
                       paste("2d_pca_", cancer_name[i], ".png", sep = ""))
  
  mygraph <- fviz_pca_ind(pca.fit, geom.ind = "point", pointshape = 21, 
               pointsize = 3, 
               fill.ind = d$group, 
               col.ind = "black", 
               palette = "jco", 
               addEllipses = TRUE,
               label = "var",               
               col.var = "black",
               repel = TRUE,
               legend.title = "Diagnosis") +
    ggtitle(paste("2D PCA-plot from", ncol(d)-1, "feature", cancer_name[i], "data")) +
    theme(plot.title = element_text(hjust = 0.5, size=20), 
          legend.title=element_text(size = 20),
          legend.text=element_text(size = 20),
          text = element_text(size=20))
  
  ggsave(mypath3, mygraph)
}

save.image("C:/Users/akim127/Google Drive/Personal/Assessing_reproducibility_veracity_PCA_results.RData")

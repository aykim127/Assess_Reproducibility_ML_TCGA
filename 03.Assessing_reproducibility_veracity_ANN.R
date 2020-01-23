# load required packages
require(TCGA2STAT)

##########################################################################################
# 0. Load data RNASeq RPKM data 
##########################################################################################
load("/extra/akim127/Project1/TCGA_all.RData")

## Choose one dataset to be d0 and run the rest of the code to get the results for 
## the chosen dataset. After that, repeat the process for other datasets
## Make sure to change the very last 4 lines of the code to save the RData file with a correct name
d0 <- blca0
# d0 <- brca0
# d0 <- lihc0
# d0 <- luad0

##########################################################################################
# 1. Split the dataset into tumor and normal groups
##########################################################################################
# Split the TCGA data into different sample types (tumor, recurrent, normal)
d1 <- SampleSplit(d0$dat)
dim(d1$primary.tumor)
dim(d1$recurrent.tumor)
dim(d1$normal)

tumor <- t(d1$primary.tumor)
normal <- t(d1$normal)

mean(tumor)
sd(tumor)
mean(normal)
sd(normal)

# Create a label column
tumor_group <- rep("tumor", nrow(tumor))
normal_group <- rep("normal", nrow(normal))

tumor <- data.frame(tumor, tumor_group)
normal <- data.frame(normal, normal_group)

colnames(tumor)[ncol(tumor)] <- "group"
colnames(normal)[ncol(tumor)] <- "group"

# Combine the tumor and normal datasets 
d <- rbind(tumor, normal)
dim(d)

# Remove columns with all zero values: 198 columns have all zero values
d <- d[, colSums(d != 0) > 10]
dim(d)

##########################################################################################
# 2. Set the number of simulation
##########################################################################################
nsim <- 1

##########################################################################################
# 6. Neural network
##########################################################################################
library(RCurl)
library(bitops)
library(rjson)
library(statmod)
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '100g')

# create null matrices to store values from the loop
ann_mod1 <- ann_results1 <- list()
ann_pred_err1 <- matrix(0, ncol = 4, nrow = nsim)
colnames(ann_pred_err1) <- c("normal_train_err", "tumor_train_err", "train_err", "test_err")

# For BLCA data
ann_train_pred1 <- ann_train_prob1 <- ann_train_y1 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
ann_test_pred1 <- ann_test_prob1 <- ann_test_y1 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# # For BRCA, LIHC, and LUAD data
# ann_train_pred1 <- ann_train_prob1 <- ann_train_y1 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# ann_test_pred1 <- ann_test_prob1 <- ann_test_y1 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

ann_mod2 <- ann_results2 <- list()
ann_pred_err2 <- matrix(0, ncol = 4, nrow = nsim)
colnames(ann_pred_err2) <- c("normal_train_err", "tumor_train_err", "train_err", "test_err")

ann_train_pred2 <- ann_train_prob2 <- ann_train_y2 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
ann_test_pred2 <- ann_test_prob2 <- ann_test_y2 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))


for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future
  set.seed(101)
  
  # selecting 60% of data as sample from total 'n' rows of the data to construct the train set
  sample <- c(sample(x = which(d$group=="tumor"), size = floor(.60*sum(d$group=="tumor")), replace = F), 
              sample(x = which(d$group=="normal"), size = floor(.60*sum(d$group=="normal")), replace = F))
  sample <- sort(sample)
  
  train <- d[sample,]
  test <- d[-sample,]
  
  # get the # of columns and rows of the train set
  p <- ncol(train)
  n <- nrow(train)
  
  x_train <- train[,-p]
  y_train <- train[,p]
  
  x_train <- as.matrix(x_train)
  y_train <- as.matrix(y_train)
  
  x_test <- test[,-p]
  y_test <- test[,p]
  
  x_test <- as.matrix(x_test)
  y_test <- as.matrix(y_test)
  
  train$group <- as.factor(train$group)
  test$group <- as.factor(test$group)
  
  # fit ANN1 with 200 nodes in a hidden layer
  ann_fit <- h2o.deeplearning(x = 1:p-1, y = p, 
                              training_frame = as.h2o(train), epochs = 10,
                              hidden = 200, adaptive_rate = T, rho = 0.8, epsilon = 1e-10, 
                              export_weights_and_biases = T, nfolds = 3, stopping_tolerance = 0.001)
  
  ann_mod1[[i]] <- ann_fit
  
  confusion_mat <- h2o.confusionMatrix(ann_fit)
  ann_results1[[i]] <- confusion_mat
  
  # make predictions
  train_pred_ann <- predict(ann_fit, as.h2o(train))$predict
  test_pred_ann <- predict(ann_fit, as.h2o(test))$predict
  
  train_prob_ann <- predict(ann_fit, as.h2o(train))$tumor
  test_prob_ann <- predict(ann_fit, as.h2o(test))$tumor
  
  ann_train_pred1[,i] <- as.vector(train_pred_ann)
  ann_test_pred1[,i] <- as.vector(test_pred_ann)
  
  ann_train_prob1[,i] <- as.vector(train_prob_ann)
  ann_test_prob1[,i] <- as.vector(test_prob_ann)
  
  # train error
  train_err <- mean(y_train!=as.vector(train_pred_ann))
  # test error 
  test_err <- mean(y_test!=as.vector(test_pred_ann))
  
  # store the train error and the test error
  ann_pred_err1[i,] <- c(train_err, test_err)
  
  # fit ANN2 with 500 nodes in a hidden layer
  ann_fit <- h2o.deeplearning(x = 1:p-1, y = p, 
                              training_frame = as.h2o(train), epochs = 10,
                              hidden = 500, adaptive_rate = T, rho = 0.8, epsilon = 1e-10, 
                              export_weights_and_biases = T, nfolds = 3, stopping_tolerance = 0.001)
  
  ann_mod2[[i]] <- ann_fit
  
  confusion_mat <- h2o.confusionMatrix(ann_fit)
  ann_results2[[i]] <- confusion_mat
  
  # make predictions
  train_pred_ann <- predict(ann_fit, as.h2o(train))$predict
  test_pred_ann <- predict(ann_fit, as.h2o(test))$predict
  
  train_prob_ann <- predict(ann_fit, as.h2o(train))$tumor
  test_prob_ann <- predict(ann_fit, as.h2o(test))$tumor
  
  ann_train_pred2[,i] <- as.vector(train_pred_ann)
  ann_test_pred2[,i] <- as.vector(test_pred_ann)
  
  ann_train_prob2[,i] <- as.vector(train_prob_ann)
  ann_test_prob2[,i] <- as.vector(test_prob_ann)
  
  # train error
  train_err <- mean(y_train!=as.vector(train_pred_ann))
  # test error 
  test_err <- mean(y_test!=as.vector(test_pred_ann))
  
  # store the train error and the test error
  ann_pred_err2[i,] <- c(train_err, test_err)
}

save.image("/extra/akim127/Project1/blca_ann.RData")
# save.image("/extra/akim127/Project1/brca_ann.RData")
# save.image("/extra/akim127/Project1/luad_ann.RData")
# save.image("/extra/akim127/Project1/lihc_ann.RData")




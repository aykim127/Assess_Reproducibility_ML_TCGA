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
# 7. Gradient boosting
##########################################################################################
library(RCurl)
library(bitops)
library(rjson)
library(statmod)
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '100g')

# create null matrices to store values from the loop
gbm_mod1 <- gbm_results1 <- list()
gbm_pred_err1 <- matrix(0, ncol = 4, nsim)
colnames(gbm_pred_err1) <- c("normal_train_err", "tumor_train_err", "train_err", "test_err")

# For BLCA data
gbm_train_pred1 <- gbm_train_prob1 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
gbm_test_pred1 <- gbm_test_prob1 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# gbm_train_pred1 <- gbm_train_prob1 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# gbm_test_pred1 <- gbm_test_prob1 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

gbm_mod2 <- gbm_results2 <- list()
gbm_pred_err2 <- matrix(0, ncol = 4, nsim)
colnames(gbm_pred_err2) <- c("normal_train_err", "tumor_train_err", "train_err", "test_err")

# For BLCA data
gbm_train_pred2 <- gbm_train_prob2 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
gbm_test_pred2 <- gbm_test_prob2 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# gbm_train_pred2 <- gbm_train_prob2 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# gbm_test_pred2 <- gbm_test_prob2 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

gbm_mod3 <- gbm_results3 <- list()
gbm_pred_err3 <- matrix(0, ncol = 4, nsim)
colnames(gbm_pred_err3) <- c("normal_train_err", "tumor_train_err", "train_err", "test_err")

# For BLCA data
gbm_train_pred3 <- gbm_train_prob3 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
gbm_test_pred3 <- gbm_test_prob3 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# gbm_train_pred3 <- gbm_train_prob3 <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# gbm_test_pred3 <- gbm_test_prob3 <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future also
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
  
  # fit GBM using the training set
  ## model #1
  gbm_fit <- h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                     learn_rate=0.1, min_rows = 3, sample_rate=0.5,
                     stopping_tolerance = 0.001, stopping_metric = "misclassification")
  
  gbm_mod1[[i]] <- gbm_fit
  
  confusion_mat <- h2o.confusionMatrix(gbm_fit)
  gbm_results1[[i]] <- confusion_mat
  
  # make predictions
  train_pred_gbm <- predict(gbm_fit, as.h2o(train))$predict
  test_pred_gbm <- predict(gbm_fit, as.h2o(test))$predict
  
  train_prob_gbm <- predict(gbm_fit, as.h2o(train))$tumor
  test_prob_gbm <- predict(gbm_fit, as.h2o(test))$tumor
  
  gbm_train_pred1[,i] <- as.vector(train_pred_gbm)
  gbm_test_pred1[,i] <- as.vector(test_pred_gbm)
  
  gbm_train_prob1[,i] <- as.vector(train_prob_gbm)
  gbm_test_prob1[,i] <- as.vector(test_prob_gbm)
  
  # train error
  train_err <- mean(y_train!=as.vector(train_pred_gbm))
  # test error
  test_err <- mean(y_test!=as.vector(test_pred_gbm))
  
  # store the train error and the test error
  gbm_pred_err1[i,] <- c(train_err, test_err)
  
  ## model #2
  gbm_fit <- h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                     learn_rate=0.01, min_rows = 3, sample_rate=0.5,
                     stopping_tolerance = 0.001, stopping_metric = "misclassification")
  
  gbm_mod2[[i]] <- gbm_fit
  
  confusion_mat <- h2o.confusionMatrix(gbm_fit)
  gbm_results2[[i]] <- confusion_mat
  
  # make predictions
  train_pred_gbm <- predict(gbm_fit, as.h2o(train))$predict
  test_pred_gbm <- predict(gbm_fit, as.h2o(test))$predict
  
  train_prob_gbm <- predict(gbm_fit, as.h2o(train))$tumor
  test_prob_gbm <- predict(gbm_fit, as.h2o(test))$tumor
  
  gbm_train_pred2[,i] <- as.vector(train_pred_gbm)
  gbm_test_pred2[,i] <- as.vector(test_pred_gbm)
  
  gbm_train_prob2[,i] <- as.vector(train_prob_gbm)
  gbm_test_prob2[,i] <- as.vector(test_prob_gbm)
  
  # train error
  train_err <- mean(y_train!=as.vector(train_pred_gbm))
  # test error
  test_err <- mean(y_test!=as.vector(test_pred_gbm))
  
  # store the train error and the test error
  gbm_pred_err2[i,] <- c(train_err, test_err)
  
  ## model #3
  gbm_fit <- h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                     learn_rate=0.001, min_rows = 3, sample_rate=0.5,
                     stopping_tolerance = 0.001, stopping_metric = "misclassification")
  
  gbm_mod3[[i]] <- gbm_fit
  
  confusion_mat <- h2o.confusionMatrix(gbm_fit)
  gbm_results3[[i]] <- confusion_mat
  
  # make predictions
  train_pred_gbm <- predict(gbm_fit, as.h2o(train))$predict
  test_pred_gbm <- predict(gbm_fit, as.h2o(test))$predict
  
  train_prob_gbm <- predict(gbm_fit, as.h2o(train))$tumor
  test_prob_gbm <- predict(gbm_fit, as.h2o(test))$tumor
  
  gbm_train_pred3[,i] <- as.vector(train_pred_gbm)
  gbm_test_pred3[,i] <- as.vector(test_pred_gbm)
  
  gbm_train_prob3[,i] <- as.vector(train_prob_gbm)
  gbm_test_prob3[,i] <- as.vector(test_prob_gbm)
  
  # train error
  train_err <- mean(y_train!=as.vector(train_pred_gbm))
  # test error
  test_err <- mean(y_test!=as.vector(test_pred_gbm))
  
  # store the train error and the test error
  gbm_pred_err3[i,] <- c(train_err, test_err)
}


save.image("/extra/akim127/Project1/blca_gbm.RData")
# save.image("/extra/akim127/Project1/brca_gbm.RData")
# save.image("/extra/akim127/Project1/luad_gbm.RData")
# save.image("/extra/akim127/Project1/lihc_gbm.RData")


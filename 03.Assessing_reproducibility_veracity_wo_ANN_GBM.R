# load required packages
require(TCGA2STAT)

library(MASS)
library(glmnet)
library(e1071)
library(randomForest)
library(rpart)

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
nsim <- 30

##########################################################################################
# 3. Fit lasso 
##########################################################################################
# create null matrices to store values and results from each loop
lasso_coef <- matrix(0, nrow=nsim, ncol=ncol(d)-1)
colnames(lasso_coef) <- colnames(d)[-ncol(d)]
lasso_pred_err <- matrix(0, ncol=3, nrow=nsim)
colnames(lasso_pred_err) <- c("lambda_min", "train_err", "test_err")

# For BLCA
lasso_train_pred <- lasso_train_prob <- sample_train_y <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
lasso_test_pred <- lasso_test_prob <- sample_test_y <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# lasso_train_pred <- lasso_train_prob <- sample_train_y <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# lasso_test_pred <- lasso_test_prob <- sample_test_y <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

sample_train_x <- sample_test_x <- list()
lasso_mod <- list()

for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future also
  set.seed(100+i)
  
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
  
  # run lines 114 to 123 for LIHC data only
  # the sample size of the LIHC is so small that the model cannot perform 
  # the 3-fold CV properfly without manually assigning fold values by running these lines 
  
  # nfold <- 3
  # 
  # # assign folds evenly using the modulus operator
  # set.seed(i)
  # fold0 <- sample.int(sum(y_train=="normal")) %% nfold
  # fold1 <- sample.int(sum(y_train=="tumor")) %% nfold
  # foldid <- numeric(length(y_train))
  # foldid[y_train=="normal"] <- fold0
  # foldid[y_train=="tumor"] <- fold1
  # foldid <- foldid + 1
  
  # perform cross-validation
  set.seed(i)
  lasso_cv <- cv.glmnet(x_train, y_train, foldid = foldid, family = "binomial", standardize = FALSE)

  # set seed for the 3-fold cross-validation
  # we may use 10-fold or 5-fold CV for a bigger sample size 
  # set.seed(i)
  # lasso_cv <- cv.glmnet(x_train, y_train, alpha=1, nfolds = 3, family="binomial", standardize = FALSE)
  lasso_mod[[i]] <- lasso_cv
  
  # choose lambda that minimizes the cv
  lamb_min <- lasso_cv$lambda.min
  
  # store coefficients
  lasso_coef[i,] <- coef(lasso_cv, s = lamb_min)[-1] # exclude the intercept
  
  # make predictions
  train_pred_lasso <- predict(lasso_cv, newx=x_train, type="class", s = lamb_min)
  test_pred_lasso <- predict(lasso_cv, newx=x_test, type="class", s = lamb_min)
  
  lasso_train_pred[,i] <- train_pred_lasso
  lasso_test_pred[,i] <- test_pred_lasso
  
  lasso_train_prob[,i] <- predict(lasso_cv, newx=x_train, type="response", s = lamb_min)
  lasso_test_prob[,i] <- predict(lasso_cv, newx=x_test, type="response", s = lamb_min)
  
  sample_train_x[[i]] <- x_train
  sample_test_x[[i]] <- x_test
  
  sample_train_y[,i] <- y_train
  sample_test_y[,i] <- y_test
  
  # train error
  train_err <- mean(y_train!=train_pred_lasso)
  # test error 
  test_err <- mean(y_test!=test_pred_lasso)
  
  # store the minimum lambda, the train error, and the test error
  lasso_pred_err[i,] <- c(lamb_min, train_err, test_err)
}

##########################################################################################
# 4. Fit elastic net 
##########################################################################################
# create null matrices to store values from the loop
enet_coef <- matrix(0, nrow=nsim, ncol=ncol(d)-1)
colnames(enet_coef) <- colnames(d)[-ncol(d)]
enet_pred_err <- matrix(0, ncol=4, nrow=nsim)
colnames(enet_pred_err) <- c("alpha", "lambda_min", "train_err", "test_err")

# For BLCA data
enet_train_pred <- enet_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
enet_test_pred <- enet_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# enet_train_pred <- enet_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# enet_test_pred <- enet_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

enet_mod <- list()

for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future also
  set.seed(100+i)
  
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
  
  # run lines 217 to 226 for LIHC data only
  # the sample size of the LIHC is so small that the model cannot perform 
  # the 3-fold CV properfly without manually assigning fold values by running these lines 
  
  # nfold <- 3
  # 
  # # assign folds evenly using the modulus operator
  # set.seed(i)
  # fold0 <- sample.int(sum(y_train=="normal")) %% nfold
  # fold1 <- sample.int(sum(y_train=="tumor")) %% nfold
  # foldid <- numeric(length(y_train))
  # foldid[y_train=="normal"] <- fold0
  # foldid[y_train=="tumor"] <- fold1
  # foldid <- foldid + 1
  
  # set the sequence of alpha values that we are going to scan through
  alpha_seq <- seq(0.9, 0.1, -0.1)
  cvm_min <- lamb <- rep(0, 9)
  
  for (k in 1:length(alpha_seq)){
    set.seed(i)
    enet_cv <- cv.glmnet(x_train, y_train, foldid = foldid, family="binomial", alpha=alpha_seq[k], standardize=FALSE)
    lamb[k] <- enet_cv$lambda.min
    cvm_min[k] <- enet_cv$cvm[which.min(enet_cv$lambda)]
  }
  
  # select alpha/lambda that minimizes mean cv
  par_index <- which.min(cvm_min)
  cvm_enet <- min(cvm_min)
  alpha_enet <- alpha_seq[par_index]
  lambda_enet <- lamb[par_index]
  
  # now, fit elastic net using the selected parameters alpha and lambda
  set.seed(i)
  enet_cv <- cv.glmnet(x_train, y_train, alpha=alpha_enet, foldid = foldid, family="binomial", standardize=FALSE)
  enet_mod[[i]] <- enet_cv
  
  # minimum cv mean and lambda that minimizes the cv mean
  enet_mcv <- enet_cv$cvm[which.min(enet_cv$cvm)]
  lamb_min <- lambda_enet
  
  # store coefficients
  enet_coef[i,] <- coef(enet_cv, s = lamb_min)[-1] # exclude the intercept
  
  # make predictions
  train_pred_enet <- predict(enet_cv, newx=x_train, type="class", s = lamb_min)
  test_pred_enet <- predict(enet_cv, newx=x_test, type="class", s = lamb_min)
  
  enet_train_pred[,i] <- train_pred_enet
  enet_test_pred[,i] <- test_pred_enet
  
  enet_train_prob[,i] <- predict(enet_cv, newx=x_train, type="response", s = lamb_min)
  enet_test_prob[,i] <- predict(enet_cv, newx=x_test, type="response", s = lamb_min)
  
  # train error
  train_err <- mean(y_train!=train_pred_enet)
  # test error 
  test_err <- mean(y_test!=test_pred_enet)
  
  # store the minimum lambda, the train error, and the test error
  enet_pred_err[i,] <- c(alpha_enet, lamb_min, train_err, test_err)
}

##########################################################################################
# 5. Fit SVM 
##########################################################################################
# create null matrices to store values from the loop
svm_pred_err <- matrix(0, ncol=4, nrow=nsim)
colnames(svm_pred_err) <- c("gamma", "cost", "train_err", "test_err")

# For BLCA data
svm_train_pred <- svm_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
svm_test_pred <- svm_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# svm_train_pred <- svm_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# svm_test_pred <- svm_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

svm_mod <- list()

for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future also
  set.seed(100+i)
  
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
  
  # tuning using the training set
  set.seed(i)
  tune <- tune.svm(x_train, as.factor(y_train), gamma = 10^(-3:3), cost = 10^(-3:3),
                   kernel="polynomial", scale=FALSE)
  par <- tune$best.parameters
  
  # fit svm using the training set
  svm_fit <- svm(x=x_train, y=as.factor(y_train), gamma = par[1], cost = par[2],
                 method = "C-classification", kernel = "polynomial", scale=FALSE, 
                 cross = 3, probability = TRUE)
  svm_mod[[i]] <- svm_fit
  
  # make predictions
  train_pred_svm <- predict(svm_fit, newdata=x_train)
  test_pred_svm <- predict(svm_fit, newdata=x_test)
  
  svm_train_pred[,i] <- train_pred_svm
  svm_test_pred[,i] <- test_pred_svm
  
  svm_train_prob[,i] <- attr(predict(svm_fit, x_train, probability=TRUE), "probabilities")[,1]
  svm_test_prob[,i] <- attr(predict(svm_fit, x_test, probability=TRUE), "probabilities")[,1]
  
  # train error
  train_err <- mean(y_train!=train_pred_svm)
  # test error
  test_err <- mean(y_test!=test_pred_svm)
  
  # store the minimum lambda, the train error, and the test error
  svm_pred_err[i,] <- c(par$gamma, par$cost, train_err, test_err)
}

##########################################################################################
# 6. Fit Random Forest
##########################################################################################
# create null matrices to store values from the loop
rf_pred_err <- matrix(0, ncol=3, nrow=nsim)
colnames(rf_pred_err) <- c("oob", "train_err", "test_err")

# For BLCA data
rf_train_pred <- rf_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d))-1)
rf_test_pred <- rf_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d))+1)

# For BRCA, LIHC, and LUAD data
# rf_train_pred <- rf_train_prob <- matrix(0, ncol=nsim, nrow=floor(.60*nrow(d)))
# rf_test_pred <- rf_test_prob <- matrix(0, ncol=nsim, nrow=nrow(d)-floor(.60*nrow(d)))

rf_top_var <- matrix(0, ncol=10, nrow=nsim)
rf_mod <- list()

for (i in 1:nsim) {
  # set seed so that same sample can be reproduced in future also
  set.seed(100+i)
  
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
  
  # fit random forest using the training set
  rf_fit <- randomForest(x_train, as.factor(y_train), importance = T, ntree = 500)
  plot(rf_fit)
  rf_mod[[i]] <- rf_fit
  
  oob_est <- rf_fit$err.rate[500,1]
  
  # make predictions
  train_pred_rf <- predict(rf_fit, x_train)
  test_pred_rf <- predict(rf_fit, x_test)
  
  rf_train_pred[,i] <- train_pred_rf
  rf_test_pred[,i] <- test_pred_rf
  
  rf_train_prob[,i] <- predict(rf_fit, x_train, type = "prob")[,2]
  rf_test_prob[,i] <- predict(rf_fit, x_test, type = "prob")[,2]
  
  var_imp <- data.frame(importance(rf_fit, type=2))
  var_imp$Variables <- row.names(var_imp)  
  var_imp <- var_imp[order(var_imp$MeanDecreaseGini, decreasing = T),]
  rf_top_var[i,] <- var_imp$Variables[1:10]
  
  # train error
  train_err <- mean(as.factor(y_train)!=train_pred_rf)
  # test error 
  test_err <- mean(as.factor(y_test)!=test_pred_rf)
  
  # store the minimum lambda, the train error, and the test error
  rf_pred_err[i,] <- c(oob_est, train_err, test_err)
}

save.image("/extra/akim127/Project1/blca.RData")
# save.image("/extra/akim127/Project1/brca.RData")
# save.image("/extra/akim127/Project1/lihc.RData")
# save.image("/extra/akim127/Project1/luad.RData")

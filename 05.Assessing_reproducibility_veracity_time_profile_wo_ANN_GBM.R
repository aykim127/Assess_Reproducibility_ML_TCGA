library(glmnet)
library(e1071)
library(randomForest)

dim.X = 100
numSamples = 100
nsim = 30

### with 100 features
lasso.100 <- enet.100 <- svm.100 <- rf.100 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)

  ##### Generate simulated data
  
  ## Lasso
  set.seed(i)
  lasso.100[i] <- system.time (cv.glmnet(X, y, nfold=3, family="binomial", standardize=FALSE))[3]

  ## Elastic net
  enet.cv <- list()
  alpha_seq <- seq(0.9, 0.1, -0.1)
  enet.100[i] <- system.time(for (j in 1:length(alpha_seq)){
      set.seed(j)
      enet.cv[[j]] <- cv.glmnet(X, y, nfold=3, family="binomial", alpha=alpha_seq[j],
                       standardize=FALSE)
      })[3]
  
  ## SVM
  set.seed(i)
  svm.100[i] <- system.time(tune.svm(X, as.factor(y), gamma=10^(-3:3), cost=10^(-3:3),
                                 kernel="polynomial", scale=FALSE))[3]
  
  ## RF
  set.seed(i)
  rf.100[i] <- system.time(randomForest(X, as.factor(y), importance=T, ntree=500))[3]
}


### with 1000 features
dim.X = 1000
lasso.1000 <- enet.1000 <- svm.1000 <- rf.1000 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)
  
  ##### Generate simulated data
  
  ## Lasso
  set.seed(i)
  lasso.1000[i] <- system.time (cv.glmnet(X, y, nfold=3, family="binomial", standardize=FALSE))[3]
  
  ## Elastic net
  enet.cv <- list()
  alpha_seq <- seq(0.9, 0.1, -0.1)
  enet.1000[i] <- system.time(for (j in 1:length(alpha_seq)){
    set.seed(j)
    enet.cv[[j]] <- cv.glmnet(X, y, nfold=3, family="binomial", alpha=alpha_seq[j],
                              standardize=FALSE)
  })[3]
  
  ## SVM
  set.seed(i)
  svm.1000[i] <- system.time(tune.svm(X, as.factor(y), gamma=10^(-3:3), cost=10^(-3:3),
                                     kernel="polynomial", scale=FALSE))[3]
  
  ## RF
  set.seed(i)
  rf.1000[i] <- system.time(randomForest(X, as.factor(y), importance=T, ntree=500))[3]
}


### with 10000 features
dim.X = 10000
lasso.10000 <- enet.10000 <- svm.10000 <- rf.10000 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)
  
  ##### Generate simulated data
  
  ## Lasso
  set.seed(i)
  lasso.10000[i] <- system.time (cv.glmnet(X, y, nfold=3, family="binomial", standardize=FALSE))[3]
  
  ## Elastic net
  enet.cv <- list()
  alpha_seq <- seq(0.9, 0.1, -0.1)
  enet.10000[i] <- system.time(for (j in 1:length(alpha_seq)){
    set.seed(j)
    enet.cv[[j]] <- cv.glmnet(X, y, nfold=3, family="binomial", alpha=alpha_seq[j],
                              standardize=FALSE)
  })[3]
  
  ## SVM
  set.seed(i)
  svm.10000[i] <- system.time(tune.svm(X, as.factor(y), gamma=10^(-3:3), cost=10^(-3:3),
                                      kernel="polynomial", scale=FALSE))[3]
  
  ## RF
  set.seed(i)
  rf.10000[i] <- system.time(randomForest(X, as.factor(y), importance=T, ntree=500))[3]
}

# Change the directory to save the RData file
save.image("/extra/akim127/Project1/time_profile_wo_ANN_GBM.RData")

dim.X = 100
numSamples = 100
nsim = 30

library(RCurl)
library(bitops)
library(rjson)
library(statmod)
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '100g')

### with 100 features
ann.100 <- gbm.100 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)
  train <- data.frame(X, y)
  p <- ncol(train)

  ##### Generate simulated data
  
  ## ANN
  ann.100[i] <- system.time(h2o.deeplearning(x = 1:p-1, y = p, 
                              training_frame = as.h2o(train), epochs = 10,
                              hidden = 500, adaptive_rate = T, rho = 0.8, epsilon = 1e-10, 
                              export_weights_and_biases = T, nfolds = 3, stopping_tolerance = 0.001))
  
  ## GBM
  gbm.100[i] <- system.time(h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                                     learn_rate=0.01, min_rows = 3, sample_rate=0.5,
                                     stopping_tolerance = 0.001, stopping_metric = "misclassification"))
}


### with 1000 features
dim.X = 1000
ann.1000 <- gbm.1000 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)
  train <- data.frame(X, y)
  p <- ncol(train)
  
  ##### Generate simulated data
  
  ## ANN
  ann.1000[i] <- system.time(h2o.deeplearning(x = 1:p-1, y = p, 
                                             training_frame = as.h2o(train), epochs = 10,
                                             hidden = 500, adaptive_rate = T, rho = 0.8, epsilon = 1e-10, 
                                             export_weights_and_biases = T, nfolds = 3, stopping_tolerance = 0.001))
  
  ## GBM
  gbm.1000[i] <- system.time(h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                                    learn_rate=0.01, min_rows = 3, sample_rate=0.5,
                                    stopping_tolerance = 0.001, stopping_metric = "misclassification"))
}

### with 10000 features
dim.X = 10000
ann.10000 <- gbm.10000 <- rep(0, nsim)

for (i in (1:nsim)) {
  set.seed(i)
  X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
  z = X%*% c(rep(3,5), rep(0,dim.X-5))
  y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
  y = as.matrix(y)
  train <- data.frame(X, y)
  p <- ncol(train)
  
  ##### Generate simulated data
  
  ## ANN
  ann.10000[i] <- system.time(h2o.deeplearning(x = 1:p-1, y = p, 
                                              training_frame = as.h2o(train), epochs = 10,
                                              hidden = 500, adaptive_rate = T, rho = 0.8, epsilon = 1e-10, 
                                              export_weights_and_biases = T, nfolds = 3, stopping_tolerance = 0.001))
  
  ## GBM
  gbm.10000[i] <- system.time(h2o.gbm(x = 1:p-1, y = p, training_frame = as.h2o(train), nfolds = 3, ntrees = 1000, 
                                     learn_rate=0.01, min_rows = 3, sample_rate=0.5,
                                     stopping_tolerance = 0.001, stopping_metric = "misclassification"))
}

# Change the directory to save the RData file
save.image("/extra/akim127/Project1/time_profile_ANN_GBM.RData")

library(precrec)
library(ggplot2)

# Load one data at a time to evaluate the classification and identifiability
# Run the rest of the code starting from line 25 for each of 12 results data listed below

# For the regression and RF models
load("C:/Users/akim127/Downloads/blca.RData")
# load("C:/Users/akim127/Downloads/brca.RData")
# load("C:/Users/akim127/Downloads/lihc.RData")
# load("C:/Users/akim127/Downloads/luad.RData")

# For ANN models
# load("C:/Users/akim127/Downloads/blca_ann.RData")
# load("C:/Users/akim127/Downloads/brca_ann.RData")
# load("C:/Users/akim127/Downloads/lihc_ann.RData")
# load("C:/Users/akim127/Downloads/luad_ann.RData")

# For GBM models
# load("C:/Users/akim127/Downloads/blca_gbm.RData")
# load("C:/Users/akim127/Downloads/brca_gbm.RData")
# load("C:/Users/akim127/Downloads/lihc_gbm.RData")
# load("C:/Users/akim127/Downloads/luad_gbm.RData")

# For ANN and GBM model set nsim to 1
nsim <- 30

##########################################################################################
# 1. Train/Test errors and some parameter values
##########################################################################################
lasso_err_mean <- apply(lasso_pred_err, 2, mean)
lasso_err_sd <- apply(lasso_pred_err, 2, sd)

enet_err_mean <- apply(enet_pred_err, 2, mean)
enet_err_sd <- apply(enet_pred_err, 2, sd)

svm_err_mean <- apply(svm_pred_err, 2, mean)
svm_err_sd <- apply(svm_pred_err, 2, sd)

rf_err_mean <- apply(rf_pred_err, 2, mean)
rf_err_sd <- apply(rf_pred_err, 2, sd)

# round(ann_pred_err1, 3)
# round(ann_pred_err2, 3)
# 
# round(gbm_pred_err1, 3)
# round(gbm_pred_err2, 3)
# round(gbm_pred_err3, 3)

##########################################################################################
# 2. sensitivity, specificity, precision, F1 scores
##########################################################################################
# change d1 to get the performance metrics for each method

d1 <- mmdata(lasso_test_prob, sample_test_y)
# d1 <- mmdata(enet_test_prob, sample_test_y)
# d1 <- mmdata(svm_test_prob, sample_test_y)
# d1 <- mmdata(rf_test_prob, sample_test_y)

# d1 <- mmdata(ann_test_prob1, y_test)
# d1 <- mmdata(ann_test_prob2, y_test)
# 
# d1 <- mmdata(gbm_test_prob1, y_test)
# d1 <- mmdata(gbm_test_prob2, y_test)
# d1 <- mmdata(gbm_test_prob3, y_test)

d1_eval <- evalmod(d1, mode = "basic")

# specificity and sensitivity plot
# autoplot(d1_eval, c("specificity", "sensitivity", "precision", "fscore"), show_legend = FALSE)
autoplot(d1_eval, "specificity", show_legend = FALSE) + geom_point(size=5) + theme(text = element_text(size=30), axis.title.x=element_blank())
autoplot(d1_eval, "sensitivity", show_legend = FALSE) + geom_point(size=5) + theme(text = element_text(size=30), axis.title.x=element_blank())
autoplot(d1_eval, "precision", show_legend = FALSE) + geom_point(size=5) + theme(text = element_text(size=30), axis.title.x=element_blank())
autoplot(d1_eval, "fscore", show_legend = FALSE) + geom_point(size=5) + theme(text = element_text(size=30), axis.title.x=element_blank())

# find median, mean, and SD of the each performance metric
d1_sp_median <- d1_sp_mean <- d1_sp_sd <- rep(0, nsim)
d1_sn_median <- d1_sn_mean <- d1_sn_sd <- rep(0, nsim)
d1_pre_median <- d1_pre_mean <- d1_pre_sd <- rep(0, nsim)
d1_f1_median <- d1_f1_mean <- d1_f1_sd <- rep(0, nsim)

for (i in 1:nsim) {
  d1_sp_median[i] <- median(d1_eval$sp[[i]]$y)
  d1_sp_mean[i] <- mean(d1_eval$sp[[i]]$y)
  d1_sp_sd[i] <- sd(d1_eval$sp[[i]]$y)
  
  d1_sn_median[i] <- median(d1_eval$sn[[i]]$y)
  d1_sn_mean[i] <- mean(d1_eval$sn[[i]]$y)
  d1_sn_sd[i] <- sd(d1_eval$sn[[i]]$y)
  
  d1_pre_median[i] <- median(d1_eval$prec[[i]]$y)
  d1_pre_mean[i] <- mean(d1_eval$prec[[i]]$y)
  d1_pre_sd[i] <- sd(d1_eval$prec[[i]]$y)
  
  d1_f1_median[i] <- median(d1_eval$fscore[[i]]$y)
  d1_f1_mean[i] <- mean(d1_eval$fscore[[i]]$y)
  d1_f1_sd[i] <- sd(d1_eval$fscore[[i]]$y)
}

round(mean(d1_sp_mean), 3)
round(mean(d1_sp_sd), 3)

round(mean(d1_sn_mean), 3)
round(mean(d1_sn_sd), 3)

round(mean(d1_pre_mean), 3)
round(mean(d1_pre_sd), 3)

round(mean(d1_f1_mean), 3)
round(mean(d1_f1_sd), 3)


##########################################################################################
# 3. model identifiability
##########################################################################################
lasso_gene <- enet_gene <- list()
lasso_gene_count <- enet_gene_count <- rep(0, nsim)

for (i in 1:nsim) {
  lasso_gene[[i]] <- which(lasso_coef[i,]!=0)
  lasso_gene_count[i] <- length(which(lasso_coef[i,]!=0))
  enet_gene[[i]] <- which(enet_coef[i,]!=0)
  enet_gene_count[i] <- length(which(enet_coef[i,]!=0))
}

lasso_gene2 <- unlist(lasso_gene)
lasso_genenames <- names(lasso_gene2)
lasso_gene_freq <- sort(table(lasso_genenames))
lasso_gene_freq

enet_gene2 <- unlist(enet_gene)
enet_genenames <- names(enet_gene2)
enet_gene_freq <- sort(table(enet_genenames))
enet_gene_freq

sort(table(rf_top_var))

# h2o.varimp(ann_mod1[[1]])
# h2o.varimp(ann_mod2[[1]])
# 
# h2o.varimp(gbm_mod1[[1]])
# h2o.varimp(gbm_mod2[[1]])
# h2o.varimp(gbm_mod3[[1]])
# 
# h2o.varimp_plot(ann_mod1[[1]])
# h2o.varimp_plot(ann_mod2[[1]])
# 
# h2o.varimp_plot(gbm_mod1[[1]])
# h2o.varimp_plot(gbm_mod2[[1]])
# h2o.varimp_plot(gbm_mod3[[1]])


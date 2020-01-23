# Assess_Reproducibility_ML_TCGA
Source code to reproduce study results in "Assessing Reproducibility and Veracity across Machine Learning Techniques in Biomedicine: A Case Study using TCGA Data."

How to run the code:
1) Obtain the TCGA data using the TCGA2STAT library in R
 - Run Assessing_reproducibility_veracity_data_all.R to prepare data

2) Perform the principal component analysis and plot the results to select datasets
 - Assessing_reproducibility_veracity_PCA.R

3) Train models
 - Run Assessing_reproducibility_veracity_wo_ANN_GBM.R to train the lasso, elastic net, SVM, and random forest(RF) models
 - Run Assessing_reproducibility_veracity_ANN.R to train the artificial neural networks(ANN)
 - Run Assessing_reproducibility_veracity_GBM.R to train the gradient boosting machines(GBM)

4) Evaluate the classification ability and model identifiability
 - Run Assessing_reproducibility_veracity_precrec.R
 
5) Perform simulations to estimate the computation times for all our models
 - For the regression and RF models: Assessing_reproducibility_veracity_time_profile_wo_ANN_GBM.R
 - For the ANN and GBM models: Assessing_reproducibility_veracity_time_profile_ANN_GBM.R
 - To combine and summarize the simulation results from all models: Assessing_reproducibility_veracity_time_profile_sum.R

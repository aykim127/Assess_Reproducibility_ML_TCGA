# Load the time profile simulation results
load("/Users/ahyoungkim/Downloads/time_profile_ANN_GBM.RData")
load("/Users/ahyoungkim/Downloads/time_profile_wo_ANN_GBM.RData")

library(reshape2)  # to use melt function to format the dataset proper for ggplot2
library(ggplot2)   # to plot

sim100 <- data.frame(lasso.100, enet.100, svm.100, rf.100, ann.100, gbm.100)
sim1000 <- data.frame(lasso.1000, enet.1000, svm.1000, rf.1000, ann.1000, gbm.1000)
sim10000 <- data.frame(lasso.10000, enet.10000, svm.10000, rf.10000, ann.10000, gbm.10000)
colnames(sim100) <- colnames(sim1000) <- 
  colnames(sim10000) <- c("Lasso", "Elastic Net", "SVM", "RF", "ANN", "GBM")

# If you want to plot the time profile results using boxplots, run lines 14-43
# In our case, since the scale difference was so large among different machine learning methods
# boxplots were not informative. 

sim100 <- melt(sim100)
sim1000 <- melt(sim1000)
sim10000 <- melt(sim10000)

ggplot(sim100) + geom_boxplot(aes(x=reorder(sim100$variable, sim100$value, FUN=median), 
                                  y=sim100$value))  +  
  theme_linedraw()+ labs(title=paste(100, 'Features'), x='',y='Time (seconds)')+
  theme(
    plot.title = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text.x = element_text(color="black", size=24, face="bold", angle = 90))

ggplot(sim1000) + geom_boxplot(aes(x=reorder(sim1000$variable, sim1000$value, FUN=median), 
                                  y=sim1000$value))  +  
  theme_linedraw()+ labs(title=paste(1000, 'Features'), x='',y='Time (seconds)')+
  theme(
    plot.title = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text.x = element_text(color="black", size=24, face="bold", angle = 90))

ggplot(sim10000) + geom_boxplot(aes(x=reorder(sim10000$variable, sim10000$value, FUN=median), 
                                   y=sim10000$value))  +  
  theme_linedraw()+ labs(title=paste(10000, 'Features'), x='',y='Time (seconds)')+
  theme(
    plot.title = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text.x = element_text(color="black", size=24, face="bold", angle = 90))

# To get the summary of the results (mean and SD of computation time for each method)
round(apply(sim100, 2, mean), 3)
round(apply(sim100, 2, sd), 3)

round(apply(sim1000, 2, mean), 3)
round(apply(sim1000, 2, sd), 3)

round(apply(sim10000, 2, mean), 3)
round(apply(sim10000, 2, sd), 3)




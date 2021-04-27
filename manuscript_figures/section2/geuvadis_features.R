#01-06-2020 # Diogo Ribeiro @ UNIL
# Script to perform a logistic regression between genes being co-expressed and several molecular features

pdf("geuvadis_features.pdf",10,10)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

# seed for the sampling for training/testing dataset
set.seed(666)

# number of digits after comma
options("scipen"=100, "digits"=2)

inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

######################
# Verify if there is missing data in dataset
######################
nas = nrow(mergedData) - nrow(na.omit(mergedData))
if (nas > 0){
  stop(paste("Rows with NAs:", nas, call.=FALSE)) }
# Verify that co-expressed and non-co-expressed gene pairs have same sample size. If not, there is a problem!
if (table(mergedData$significant)[1] != table(mergedData$significant)[2]) {
  stop(paste("Number of pairs is different", table(mergedData$significant)[1], table(mergedData$significant)[2]), call.=FALSE) }

#########
# Define training/test set
#########

data = get_train_test_set(mergedData)
trainData = data[1][[1]]
testData = data[2][[1]]

# cheat to get 2nd randomisation
data = get_train_test_set(mergedData)
trainData = data[1][[1]]
testData = data[2][[1]]

#########
# Perform logistic regressions and see how they perform
#########

discolored = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#8dd3c7","darkred","black")

# Naming models M1, M2 etc will ensure their nice ordering when plotting
models = list(
  get_perf(glm(significant ~ distance, data = trainData, family = "binomial"), "M1: Distance", testData)
  ,get_perf(glm(significant ~ totalCTCF + invertedCTCF, data = trainData, family = "binomial"), "M2: CTCF", testData)
  ,get_perf(glm(significant ~ tssContact, data = trainData, family = "binomial"), "M3: Hi-C", testData)
  ,get_perf(glm(significant ~ totalEnhancers + sharedEnhancers, data = trainData, family = "binomial"), "M4: Enhancer", testData)
  ,get_perf(glm(significant ~ totalTFBS + sharedTF, data = trainData, family = "binomial"), "M5: TF", testData)
  ,get_perf(glm(significant ~ diffExpr + diffCoef, data = trainData, family = "binomial"), "M6: Expression", testData)
  ,get_perf(glm(significant ~ goSharing, data = trainData, family = "binomial"), "M7: GO", testData)
  ,get_perf(glm(significant ~ LD_R2, data = trainData, family = "binomial"), "M8: LD", testData)
  ,get_perf(glm(significant ~ eqtlSharing + eGenes, data = trainData, family = "binomial"), "M9: eQTL", testData)
  ,get_perf(glm(significant ~ distance + totalCTCF + invertedCTCF + tssContact + totalEnhancers + sharedEnhancers + totalTFBS + sharedTF + LD_R2 + diffExpr + diffCoef + goSharing, data = trainData, family = "binomial"), "X: All except eQTL", testData)
  ,get_perf(glm(significant ~ distance + totalCTCF + invertedCTCF + tssContact + totalEnhancers + sharedEnhancers + totalTFBS + sharedTF + eqtlSharing + eGenes + LD_R2 + diffExpr + diffCoef + goSharing, data = trainData, family = "binomial"), "X: All together", testData)
)

# Function to plot ROC curve
plot_ROC = function(models){
  library(ggplot2)
  rocPlot = ggplot() + geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 3.5)
  count = 0
  for (m in models){ # add a line for each model
    if (count < length(models) - 1){
      rocPlot = rocPlot + geom_line(data = models[[count+1]], aes(x = fpr, y = tpr ), color = "black", size = 2) +
        geom_line(data = m, aes(x = fpr, y = tpr, color = paste(model,": ",round(auc,2), sep = "" )  ), size = 1.4)
      count = count + 1}
  }
  rocPlot = rocPlot + geom_line(data = models[[count+1]], aes(x = fpr, y = tpr), color = "black", size = 3.0) +
    geom_line(data = models[[count+1]], aes(x = fpr, y = tpr, color = paste(model,": ",round(auc,2), sep = "" )  ), size = 2.5)
  rocPlot + ylab("Sensitivity") + xlab("1-Specificity") +
    scale_color_manual(values = c(discolored[1:count], discolored[count+1]), name = paste("Model"),
                       labels = c("Distance", "CTCF", "Hi-C", "Enhancer", "TF", "Expression", "GO", "LD", "eQTL", "All except eQTL", "All together")) +
    annotate(geom = "text", x = 0.93, y = 0.46, label = " AUC", size = 8) +
    annotate(geom = "text", x = 0.93, y = 0.22, label = "0.50\n0.57\n0.51\n0.64\n0.61\n0.65\n0.54\n0.52\n0.74\n0.72\n0.82", size = 7.5, lineheight = .825) +
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
                             axis.text.x = element_text(angle = 0, vjust=0.6),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             legend.position = c(0.67, 0.27), legend.text=element_text(size=22), legend.title=element_text(size = 24), 
                             axis.line = element_line(colour = "black", size = 1),
                             panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)
}

plot_ROC(models)


lm = glm(significant ~ distance + totalCTCF + invertedCTCF + tssContact + totalEnhancers + sharedEnhancers + totalTFBS + sharedTF + diffExpr + diffCoef + LD_R2 + goSharing + eqtlSharing + eGenes, data = mergedData, family = "binomial")
summary(lm)

summary(aov(significant ~ distance + totalCTCF + invertedCTCF + tssContact + totalEnhancers + sharedEnhancers + totalTFBS + sharedTF + diffExpr + diffCoef + LD_R2 + goSharing + eqtlSharing + eGenes, data = mergedData))


dev.off()

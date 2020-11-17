#01-06-2020 # Diogo Ribeiro @ UNIL
# Script to plot distance AUC before distance-matching procedure

pdf("distance_roc.pdf",10,10)

library(data.table)
library(ggplot2)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

# seed for the sampling for training/testing dataset
set.seed(666)

# number of digits after comma
options("scipen"=100, "digits"=2)

########
# Load/process data
########

# inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/CODer_raw_results.bed"
# copFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/CODer_distance_controlled_null.bed_positive"
inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_raw_results.bed"
copFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_distance_controlled_null.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
copData = fread( copFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# COP and non-COP definition
copData = copData[significant == 1]
mergedData$significant = 0
mergedData[pairID %in% copData$pairID]$significant = 1

mergedData = unique(mergedData[,.(pairID,significant, distance)])

# Getting same sample size between positives and negatives
mergedData = rbind(mergedData[significant == 0][sample(nrow(mergedData[significant == 1]))], mergedData[significant == 1])

table(mergedData$significant)

######################
# Verify if there is missing data in your dataset (will be discussed later)
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

data = get_train_test_set_no_dist(mergedData)
trainData = data[1][[1]]
testData = data[2][[1]]

#########
# Perform logistic regressions and see how they perform
#########

discolored = c("#e41a1c")
colored = c("#e41a1c")

# Naming models M1, M2 etc will ensure their nice ordering when plotting
models = list(
  get_perf(glm(significant ~ distance, data = trainData, family = "binomial"), "Distance", testData)
)

# Function to plot ROC curve
plot_ROC = function(models){
  library(ggplot2)
  rocPlot = ggplot() + geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 3.5)
  count = 0
  for (m in models){ # add a line for each model
    if (count < length(models) - 1){
      rocPlot = rocPlot + geom_line(data = models[[count+1]], aes(x = fpr, y = tpr ), color = "black", size = 2.5) +
        geom_line(data = m, aes(x = fpr, y = tpr, color = paste(model,": ",round(auc,2), sep = "" )  ), size = 2)
      count = count + 1}
  }
  rocPlot = rocPlot + geom_line(data = models[[count+1]], aes(x = fpr, y = tpr), color = "black", size = 2.5) +
    annotate("text", label = "AUC = 0.81", x = 0.8, y = 0.3, size = 8) +
    geom_line(data = models[[count+1]], aes(x = fpr, y = tpr, color = paste(model,": ",round(auc,2), sep = "" )  ), size = 2.5)

  rocPlot + ylab("Sensitivity") + xlab("1-Specificity") +
    scale_color_manual(values = c(discolored[1:count], colored[count+1]), name = "AUC") +
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
                             axis.text.x = element_text(angle = 0, vjust=0.6),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             axis.line = element_line(colour = "black", size = 1),
                             panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, legend.position = "none")
}

plot_ROC(models)


dev.off()  

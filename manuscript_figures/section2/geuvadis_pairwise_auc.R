#01-06-2020 # Diogo Ribeiro @ UNIL
# Script to perform a logistic regression between genes being co-expressed and several molecular features

pdf("geuvadis_pairwise_auc.pdf",10,10)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

# seed for the sampling for training/testing dataset
set.seed(666)

# number of digits after comma
options("scipen"=100, "digits"=2)

########
# Load/process data
########

# wantedCorrSign = "+"
inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")


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

data = get_train_test_set(mergedData)
trainData = data[1][[1]]
testData = data[2][[1]]

#########
# Getting pairwise metric AUCs
#########

metrics = c("distance","tssContact","totalTFBS","sharedTF","totalCTCF","invertedCTCF","diffExpr",
            "diffCoef","LD_R2","sharedEnhancers","totalEnhancers","goSharing")
#"eqtlSharing","eGenes",

aucMatrix = data.table()
for (metric1 in metrics){
  for (metric2 in metrics){
    if (metric1 != metric2){
      print(paste(metric1,"&",metric2))
      resultDF = auc_rand(mergedData, metric1, metric2, nrands = 50)
      aucMatrix = rbind(aucMatrix, data.table(metric1 = metric1, metric2 = metric2, auc1 = resultDF$m1Auc, auc2 = resultDF$m2Auc, aucPair = resultDF$mPairAuc))
    }
  }
}


plot_pairwise_AUC_matrix = function(auxMatrix){
  
  # Get increase in AUC due to combination of metrics
  aucMatrix$maxAuc = apply(aucMatrix[,c("auc1","auc2")], 1, max)
  aucMatrix$aucDiff = aucMatrix$aucPair - aucMatrix$maxAuc
  
  # Plot with AUC values
  ggplot(aucMatrix, aes(x = metric1, y = metric2, fill = aucPair)) + 
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_text(aes(label = round(aucPair,2)), size = 4.5) +
    scale_fill_gradient(low = "#f7f7f7", high = "#3288bd", limits=c(0.49,max(aucMatrix$aucPair)), name = "AUC") + 
    scale_y_discrete(limits = c("distance","totalCTCF","invertedCTCF","tssContact","sharedEnhancers","totalEnhancers",
                                "sharedTF","totalTFBS","diffExpr","diffCoef","goSharing","LD_R2"),
                     labels = c("Distance","Total CTCF sites","Inverted CTCF sites","Hi-C contacts","Enhancer sharing","Total enhancers",
                                "Shared TFs","Total TFBS","Diff. expr. level","Diff. coef. var.","GO term sharing","Linkage disequilibrium")) +
    scale_x_discrete(limits = c("distance","totalCTCF","invertedCTCF","tssContact","sharedEnhancers","totalEnhancers",
                                "sharedTF","totalTFBS","diffExpr","diffCoef","goSharing","LD_R2"),
                     labels = c("Distance","Total CTCF sites","Inverted CTCF sites","Hi-C contacts","Enhancer sharing","Total enhancers",
                                "Shared TFs","Total TFBS","Diff. expr. level","Diff. coef. var.","GO term sharing","Linkage disequilibrium")) +
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust=1),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=14), 
                             aspect.ratio = 1)
}

plot_pairwise_AUC_matrix(aucMatrix)

dev.off()


#23-09-2019 # Diogo Ribeiro @ UNIL
# Script to perform a logistic regression between genes being co-expressed and several molecular features
library(data.table)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

# seed for the sampling for training/testing dataset
set.seed(666)

# number of digits after comma
options("scipen"=100, "digits"=2)

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
# Perform logistic regressions and see how they perform
# Getting AUC for many randomisations
#########
resultDF = data.table()
for (i in seq(1:50)){
  print(paste("Rand",i))
  data = get_train_test_set(mergedData)
  # data = get_train_test_set(mergedData, distance_below_cutoff = 200000)
  # data = get_train_test_set(mergedData, distance_above_cutoff = 200000)
  trainData = data[1][[1]]
  testData = data[2][[1]]
  
  # copy above models here
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
  
  
  aucVector = c()
  modelVector = c()
  for (m in models){
    modelVector = c(modelVector,m$model[1])
    aucVector = c(aucVector,m$auc[1]) }
  resultDF = rbind(resultDF, data.table(auc = aucVector, model = modelVector))
}

labels = c("1: Distance", "2: CTCF","3: Hi-C","4: Enhancer","5: TF","6: Expression", "7: GO sharing","8: LD","9: Genetic", "X: all together")
plot_rand_ROC(resultDF, labels)


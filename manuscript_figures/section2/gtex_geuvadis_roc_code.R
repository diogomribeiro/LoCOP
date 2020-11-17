#02-June-2020 # Diogo Ribeiro @ UNIL
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
# Define training/test set
#########

data = get_train_test_set(mergedData)
# data = get_train_test_set(mergedData, distance_below_cutoff = 500000)
# data = get_train_test_set(mergedData, distance_above_cutoff = 500000)
trainData = data[1][[1]]
testData = data[2][[1]]

#########
# Getting AUC for many randomisations
#########
resultDF = data.table()
for (i in seq(1:50)){
  print(paste("Rand",i))
  data = get_train_test_set(mergedData)
  # data = get_train_test_set(mergedData, distance_below_cutoff = 500000)
  # data = get_train_test_set(mergedData, distance_above_cutoff = 500000)
  trainData = data[1][[1]]
  testData = data[2][[1]]
  
  models = list(
    get_perf(glm(significant ~ distance, data = trainData, family = "binomial"), "1: Distance", testData)
    ,get_perf(glm(significant ~ totalCTCF + invertedCTCF, data = trainData, family = "binomial"), "2: CTCF", testData)
    ,get_perf(glm(significant ~ totalEnhancers + sharedEnhancers, data = trainData, family = "binomial"), "3: Enhancer", testData)
    ,get_perf(glm(significant ~ totalTFBS + sharedTF, data = trainData, family = "binomial"), "4: TF", testData)
    ,get_perf(glm(significant ~ diffExpr, data = trainData, family = "binomial"), "5: Expression", testData)
    ,get_perf(glm(significant ~ goSharing, data = trainData, family = "binomial"), "6: GO", testData)
    ,get_perf(glm(significant ~ eqtlSharing + eGenes, data = trainData, family = "binomial"), "7: eQTL", testData)
    ,get_perf(glm(significant ~ distance + totalTFBS + sharedTF + totalCTCF + invertedCTCF + totalEnhancers + sharedEnhancers + eqtlSharing + eGenes + goSharing + diffExpr, data = trainData, family = "binomial"), "X: All together", testData)
  )
  
  aucVector = c()
  modelVector = c()
  for (m in models){
    modelVector = c(modelVector,m$model[1])
    aucVector = c(aucVector,m$auc[1]) }
  resultDF = rbind(resultDF, data.table(auc = aucVector, model = modelVector))
}

labels = c("1: Distance","2: CTCF", "3: Enhancer","4: TF","5: Expression","6: GO","7: eQTL", "X: all together")
plot_rand_ROC(resultDF, labels)

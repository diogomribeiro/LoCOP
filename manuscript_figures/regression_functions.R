#!/usr/bin/env Rscript
# Useful R functions related to regression analysis / ROC curves

#######################
# Function to evaluate performance of a model
#######################
get_perf = function(model, name, testData){
  library(ROCR)
  predicted = predict(model, newdata=testData, type="response")
  pred = prediction(predicted, testData$significant)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  df = data.table(fpr = unlist(perf@x.values), tpr = unlist(perf@y.values))
  df$model = name
  df$auc = performance(pred, measure = "auc")@y.values[[1]]
  return(df)
}


#######################
# Function to plot ROC curve
#######################
plot_ROC = function(models, text){
  library(ggplot2)
  rocPlot = ggplot() + geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2, size = 1.5)
  count = 0
  for (m in models){ # add a line for each model
    rocPlot = rocPlot + geom_line(data = m, aes(x = fpr, y = tpr, color = paste(model,": ",round(auc,2), sep = "" )  ), size = 1.3) }
  rocPlot + annotate("text", label = text, x = 0.15, y = 0.9, size = 6) +
    ggtitle("ROC curve") + ylab("Sensitivity") + xlab("1-Specificity") + scale_color_brewer(palette = "Set3", name = "Model: AUC") + 
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), axis.text.x = element_text(angle = 0, vjust=0.6),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             legend.position = c(0.8, 0.3), legend.text=element_text(size=14))
}


#######################
# Function to plot AUC averages for each metric
#######################
plot_rand_ROC = function(resultDF,wantedLabels){
  # Boxplot + jitter
  meansDF = aggregate(resultDF$auc, list(resultDF$model), mean)
  colnames(meansDF) = c("model","auc")
  text = paste("Positives:",table(trainData$significant)[2]+table(testData$significant)[2], "\n Negatives:", table(trainData$significant)[1]+table(testData$significant)[1])
  # text = paste("Positives:",table(trainData$significant)[2], "\nNegatives:", table(trainData$significant)[1])
  ggplot( resultDF, aes(x = model, y = auc, fill = model ))  +
    geom_boxplot( width = 0.8) +
    geom_jitter( width = 0.2, height = 0, size = 1, alpha = 0.5) +
    geom_abline(slope = 0, intercept = 0.5, color = "red") +
    geom_text(data = meansDF, aes(x = model, y = auc-0.03, label = round(auc,2)), size = 7, fontface = "bold", color = "#525252" ) +
    annotate("text", label = text, x = 0.1, y = 0.8, size = 10, hjust = 0) +
    ylab("AUC") +
    xlab("Models") +
    scale_fill_brewer(palette="Set3",
      labels = wantedLabels) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank())

}


#######################
# Function to split training and test sets
#######################
get_train_test_set = function(dataset, percent_training = 0.8,  distance_above_cutoff = 0, distance_below_cutoff = 1000000) {
  
  # apply distance filtering
  dataset = dataset[distance >= distance_above_cutoff]  
  dataset = dataset[distance < distance_below_cutoff]
  
  # Split dataset into significant and non-significant
  datasetSig = dataset[significant == 1]
  datasetNonSig = dataset[significant == 0]
  
  # Split test and training sets
  trainRowsSig = sample(nrow(datasetSig),nrow(datasetSig)*percent_training)
  trainDataSig = datasetSig[trainRowsSig,]
  testDataSig = datasetSig[-trainRowsSig,]

  trainDataNonSig = datasetNonSig[nullId %in% trainDataSig$nullId]
  testDataNonSig = datasetNonSig[nullId %in% testDataSig$nullId]

  # merge data again
  trainData = rbind(trainDataSig, trainDataNonSig)
  testData = rbind(testDataSig, testDataNonSig)
  
  # verify same number of significant and non-significant
  if (table(trainData$significant)[1] != table(trainData$significant)[2]) {
      stop(paste("Number of pairs is different", table(trainData$significant)[1], table(trainData$significant)[2]), call.=FALSE) }
  # verify same number of significant and non-significant
  if (table(testData$significant)[1] != table(testData$significant)[2]) {
    stop(paste("Number of pairs is different", table(testData$significant)[1], table(testData$significant)[2]), call.=FALSE) }
  
  return(list(trainData,testData))
}


#######################
# Function to split training and test sets, without using distance and null ID
#######################
get_train_test_set_no_dist = function(dataset, percent_training = 0.8) {
    
  # Split dataset into significant and non-significant
  datasetSig = dataset[significant == 1]
  datasetNonSig = dataset[significant == 0]
  
  # Split test and training sets
  trainRowsSig = sample(nrow(datasetSig),nrow(datasetSig)*percent_training)
  trainDataSig = datasetSig[trainRowsSig,]
  testDataSig = datasetSig[-trainRowsSig,]

  trainRowsNonSig = sample(nrow(datasetNonSig),nrow(datasetNonSig)*percent_training)
  trainDataNonSig = datasetNonSig[trainRowsNonSig,]
  testDataNonSig = datasetNonSig[-trainRowsNonSig,]

  # merge data again
  trainData = rbind(trainDataSig, trainDataNonSig)
  testData = rbind(testDataSig, testDataNonSig)
  
  # verify same number of significant and non-significant
  if (table(trainData$significant)[1] != table(trainData$significant)[2]) {
      stop(paste("Number of pairs is different", table(trainData$significant)[1], table(trainData$significant)[2]), call.=FALSE) }
  # verify same number of significant and non-significant
  if (table(testData$significant)[1] != table(testData$significant)[2]) {
    stop(paste("Number of pairs is different", table(testData$significant)[1], table(testData$significant)[2]), call.=FALSE) }
  
  return(list(trainData,testData))
}

#######################
# Function to perform regression models with two metrics and their combination X times and report mean AUC
#######################
auc_rand = function(mergedData,metric1,metric2,nrands = 20){
  resultDF = data.table()
  metric1Auc = data.table()
  metric2Auc = data.table()
  metricPairAuc = data.table()
  for (i in seq(1:nrands)){
    data = get_train_test_set_no_dist(mergedData)
    trainData = data[1][[1]]
    testData = data[2][[1]]
    
    m1Auc = get_perf(glm(significant ~ get(metric1), data = trainData, family = "binomial"), "", testData)$auc[1]
    m2Auc = get_perf(glm(significant ~ get(metric2), data = trainData, family = "binomial"), "", testData)$auc[2]
    mPairAuc = get_perf(glm(significant ~ get(metric1) + get(metric2), data = trainData, family = "binomial"), "", testData)$auc[3]
    
    metric1Auc = rbind(metric1Auc, data.table(auc = m1Auc))
    metric2Auc = rbind(metric2Auc, data.table(auc = m2Auc))
    metricPairAuc = rbind(metricPairAuc, data.table(auc = mPairAuc))
    
  }    
  resultDF = rbind(resultDF, data.table(m1Auc = mean(metric1Auc$auc), m2Auc = mean(metric2Auc$auc), mPairAuc = mean(metricPairAuc$auc)))
  return(resultDF)
}


#######################
# Function to plot a matrix of logistic regression AUC results of pairwise metric combinations
#######################
plot_pairwise_AUC_matrix = function(auxMatrix){
  
  # Get increase in AUC due to combination of metrics
  aucMatrix$maxAuc = apply(aucMatrix[,c("auc1","auc2")], 1, max)
  aucMatrix$aucDiff = aucMatrix$aucPair - aucMatrix$maxAuc
  
  # Plot with AUC values
  ggplot(aucMatrix, aes(x = metric1, y = metric2, fill = aucPair)) + 
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_text(aes(label = round(aucPair,2))) +
    scale_fill_gradient(low = "#f7f7f7", high = "#3288bd", limits=c(0.49,max(aucMatrix$aucPair))) + 
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), axis.text.x = element_text(angle = 90, hjust=1),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=14))
}


#######################
# Function to plot a matrix of logistic regression AUC results of pairwise metric combinations
#######################
plot_pairwise_AUC_diff_matrix = function(auxMatrix){
  # Get increase in AUC due to combination of metrics
  aucMatrix$maxAuc = apply(aucMatrix[,c("auc1","auc2")], 1, max)
  aucMatrix$aucDiff = aucMatrix$aucPair - aucMatrix$maxAuc
  # Plot with increase/decrease AUC value due to combination
  ggplot(aucMatrix, aes(x = metric1, y = metric2, fill = aucDiff)) + 
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_text(aes(label = round(aucDiff,2))) +
    scale_fill_gradient2(low = "#d73027", mid = "#f7f7f7", high = "#1a9850") + 
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), axis.text.x = element_text(angle = 90, hjust=1),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=14))

}

#######################
# Function to plot a matrix of pairwise correlation between metrics
#######################
plot_pairwise_corr_matrix = function(corMatrix){
  ggplot(corMatrix, aes(x = metric1, y = metric2, fill = cor)) + 
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_text(aes(label = round(cor,2))) +
    scale_fill_gradient2(low = "#d73027", mid = "#ffffff", high = "#1a9850") + 
    theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), axis.text.x = element_text(angle = 90, hjust=1),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=14))
}


# Function to plot several ggplots together
grid_arrange_shared_legend <- function(plots) {
  require(grid)
  require(gridExtra)
  
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),legend,ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


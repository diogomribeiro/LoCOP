#02-June-2020 # Diogo Ribeiro @ UNIL
# Script to perform molecular feature result comparisons

pdf("gtex_geuvadis_auc_boxplots.pdf",13,6)

library(data.table)
library(ggplot2)

########
# INPUT FILES
########

filename = "multi_features_distance_matched.bed"
scriptname = "/users/dribeir1/code/cod/src/cod/paper_figures/section2/gtex_geuvadis_roc_code.R"

########
# DATASET 1
########

# Chose file and processing
tissue = "Lung"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Run and store
mergedData1 = mergedData
source(scriptname)
dataset1 = resultDF
dataset1Name = "1: Lung"
dataset1$dataset = dataset1Name

########
# DATASET 2
########

# Chose file and processing
tissue =  "Muscle_Skeletal"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Run and store
mergedData2 = mergedData
source(scriptname)
dataset2 = resultDF
dataset2Name = "2: Muscle Skeletal"
dataset2$dataset = dataset2Name

########
# DATASET 3
########

# Chose file and processing
tissue = "Cells_EBV-transformed_lymphocytes"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Run and store
mergedData3 = mergedData
source(scriptname)
dataset3 = resultDF
dataset3Name = "3: LCL (GTEx)"
dataset3$dataset = dataset3Name

########
# DATASET 4
########

# Chose file and processing
inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Run and store
mergedData4 = mergedData
source(scriptname)
dataset4 = resultDF
dataset4Name = "4: LCL (Geuvadis)"
dataset4$dataset = dataset4Name

# Merge datasets
# resultDF = rbind(dataset1,dataset2,dataset3)
resultDF = rbind(dataset1,dataset2,dataset3,dataset4)

##########
# Plot all AUC results
##########

## 4 datasets
meansDF = data.table(aggregate(resultDF$auc, list(resultDF$model,resultDF$dataset), mean))
colnames(meansDF) = c("model","dataset","auc")
ggplot( resultDF, aes(x = model, y = auc, fill = dataset ))  +
  geom_boxplot( width = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +
  geom_text(data = meansDF[dataset == dataset1Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = 2.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset2Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = 1.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset3Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = -0.0, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset4Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = -1.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  ylim(0.47,0.83) +
  ylab("AUC") +
  xlab("Models") +
  # coord_flip() +
  scale_fill_brewer(palette="Set3", labels = c("Lung (N = 4398)","Muscle Skeletal (N = 5401)", "LCL (GTEx) (N = 4702)","LCL (Geuvadis) (N = 6668)")) +
  scale_x_discrete(labels = c("Distance","CTCF","Enhancer","TF","Expression","GO","eQTL","All together")) + 
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=22), 
        axis.text.x = element_text(vjust=0.6),
        legend.position = c(0.20, 0.83), legend.text=element_text(size=18), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


dev.off()

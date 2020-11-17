#01-Jun-2020 # Diogo Ribeiro @ UNIL
# Script to perform molecular feature result comparisons

pdf("geuvadis_distance_split_auc_boxplots.pdf",16,8)

library(data.table)
library(ggplot2)

########
# INPUT FILES
########

filename = "multi_features_distance_matched.bed" #"current_multi_features_distance_matched_samefunction.bed" #"current_multi_features_distance_matched_clean_filtered.bed"
scriptname = "/users/dribeir1/code/cod/src/cod/paper_figures/section2/geuvadis_roc_code.R"

########
# DATASET 1
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/geuvadis/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Run and store
mergedData1 = mergedData
source(scriptname)
dataset1 = resultDF
dataset1Name = "1: All COPs"
dataset1$dataset = dataset1Name

########
# DATASET 2
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/geuvadis/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIds = mergedData[distance < 200000][significant == 1]$nullId
mergedData = mergedData[nullId %in% nullIds]

# Run and store
mergedData2 = mergedData
source(scriptname)
dataset2 = resultDF
dataset2Name = "2: [0-200]kb"
dataset2$dataset = dataset2Name

########
# DATASET 3
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/geuvadis/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIds = mergedData[distance >= 200000][significant == 1]$nullId
mergedData = mergedData[nullId %in% nullIds]

# Run and store
mergedData3 = mergedData
source(scriptname)
dataset3 = resultDF
dataset3Name = "3: [200-1000]kb"
dataset3$dataset = dataset3Name

##########
# Plot all AUC results
##########

resultDF = rbind(dataset1,dataset2,dataset3)

# 3 datasets
meansDF = data.table(aggregate(resultDF$auc, list(resultDF$model,resultDF$dataset), mean))
colnames(meansDF) = c("model","dataset","auc")
ggplot( resultDF, aes(x = model, y = auc, fill = dataset ))  +
  geom_boxplot( width = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +
  # annotate(geom = "text", x = 1, y = 0.8, label = paste("D1: ",nrow(mergedData1)/2, "\nD2: ",nrow(mergedData2)/2,"\nD3: ",nrow(mergedData3)/2,sep=""),hjust = 0, size = 6, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset1Name], aes(x = model, y = auc-0.02, label = round(auc,2)), hjust = 1.6, size = 6, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset2Name], aes(x = model, y = auc-0.02, label = round(auc,2)), hjust = 0.5, size = 6, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset3Name], aes(x = model, y = auc-0.02, label = round(auc,2)), hjust = -0.6, size = 6, fontface = "bold", color = "#525252" ) +
  # ylim(0.4,1.0) +
  ylab("AUC") +
  xlab("Models") +
  # coord_flip() +
  scale_fill_brewer(palette="Set3", labels = c("All COPs (N = 6668)","[0-200]kb (N = 4858)", "[200-1000]kb (N = 1810)")) +
  scale_x_discrete(labels = c("Distance","CTCF","Hi-C","Enhancer","TF","Expression","GO","LD","eQTL","All except eQTL","All together")) + 
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
                           axis.text.x = element_text(angle = 20, vjust=0.6),
                           legend.position = c(0.15, 0.9), legend.text=element_text(size=20), legend.title=element_blank(),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           axis.line = element_line(colour = "black", size = 1),
                           panel.border = element_rect(colour = "black", fill=NA, size=1))


dev.off()


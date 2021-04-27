#02-June-2020 # Diogo Ribeiro @ UNIL
# Script to perform molecular feature result comparisons

pdf("gtex_geuvadis_tf_enhancer_boxplots.pdf", 10, 10)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

########
# INPUT FILES
########

filename = "multi_features_distance_matched.bed"

tissue = "Lung"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData1 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset1Name = "1: Lung"

# Chose file and processing
tissue =  "Muscle_Skeletal"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData2 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2Name = "2: Muscle Skeletal"

tissue = "Cells_EBV-transformed_lymphocytes"
inFile = paste("cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData3 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset3Name = "3: LCL (GTEx)"

inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData4 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset4Name = "4: LCL (Geuvadis)"

##########
# Plot enhancer and TF details
##########

mergedData1 = mergedData1[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData2 = mergedData2[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData3 = mergedData3[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData4 = mergedData4[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]

mergedData1$dataset = dataset1Name
mergedData2$dataset = dataset2Name
mergedData3$dataset = dataset3Name
mergedData4$dataset = dataset4Name

mergedData = rbind(mergedData1,mergedData2,mergedData3,mergedData4)

##########
# Enhancers
##########

p1 = wilcox.test(mergedData[dataset == dataset1Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset1Name][significant == 0]$totalEnhancers)$p.value
p2 = wilcox.test(mergedData[dataset == dataset2Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset2Name][significant == 0]$totalEnhancers)$p.value
p3 = wilcox.test(mergedData[dataset == dataset3Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset3Name][significant == 0]$totalEnhancers)$p.value
p4 = wilcox.test(mergedData[dataset == dataset4Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset4Name][significant == 0]$totalEnhancers)$p.value
pvals = data.table(pvals = c(p1,p2,p3,p4), dataset = c(dataset1Name,dataset2Name,dataset3Name,dataset4Name))

meansDF = data.table(aggregate(mergedData$totalEnhancers, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )
g1 = ggplot( data = mergedData, aes(x = dataset, y = totalEnhancers, fill = as.factor(significant)) ) +
  # geom_violin( width = 0.6, size = 1) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 6.5, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 6.5, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 38.5, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 7 ) +
  ylim(0,40) +
  # ggtitle("Enhancers") +
  ylab("Total enhancers") +
  xlab("Dataset") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  scale_x_discrete(labels = c("Lung","Muscle Skeletal","LCL (GTEx)", "LCL (Geuvadis)")) + 
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30),
        axis.text.x = element_text(angle = 20, vjust=0.6), 
        legend.title=element_blank(), legend.key = element_rect(size = 6), 
        legend.key.size = unit(2, 'lines'), legend.position = c(0.3,0.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
  )

##########
# TFs
##########
p1 = wilcox.test(mergedData[dataset == dataset1Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset1Name][significant == 0]$totalEnhancers)$p.value
p2 = wilcox.test(mergedData[dataset == dataset2Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset2Name][significant == 0]$totalEnhancers)$p.value
p3 = wilcox.test(mergedData[dataset == dataset3Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset3Name][significant == 0]$totalEnhancers)$p.value
p4 = wilcox.test(mergedData[dataset == dataset4Name][significant == 1]$totalEnhancers, mergedData[dataset == dataset4Name][significant == 0]$totalEnhancers)$p.value
pvals = data.table(pvals = c(p1,p2,p3,p4), dataset = c(dataset1Name,dataset2Name,dataset3Name,dataset4Name))

meansDF = data.table(aggregate(mergedData$totalTFBS, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )
g2 = ggplot( data = mergedData, aes(x = dataset, y = totalTFBS, fill = as.factor(significant)) ) +
  # geom_violin( width = 0.6, size = 1) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 6.5, fontface = "bold", color = "#525252", nudge_x = 0.42) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 6.5, fontface = "bold", color = "#525252", nudge_x = -0.42) +
  geom_text(data = pvals, aes(x = dataset, y = 1200, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 7 ) +
  # ggtitle("Transcription factors") +
  ylab("Total TFBS") +
  xlab("Dataset") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  scale_x_discrete(labels = c("Lung","Muscle Skeletal","LCL (GTEx)", "LCL (Geuvadis)")) + 
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30),
        axis.text.x = element_text(angle = 20, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.title=element_blank(), legend.key = element_rect(size = 6), 
        legend.key.size = unit(2, 'lines'), legend.position = c(0.3,0.85),
  )

# grid.arrange(g1,g2,g3,g4, nrow = 2)
# grid.arrange(g1,g2, nrow = 2)
g1
g2

dev.off()

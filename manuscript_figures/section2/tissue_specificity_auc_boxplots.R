#17-June-2020 # Diogo Ribeiro @ UNIL
# Script to perform molecular feature result comparisons

args = commandArgs(trailingOnly=TRUE)

tissue = args[1] #"Cells_EBV-transformed_lymphocytes"  #"Lung"  #"Muscle_Skeletal"

pdf(paste("tissue_specificity_auc_boxplots_",tissue,".pdf",sep=""),15,15)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

########
# INPUT FILES
########

filename = "multi_features_distance_matched.bed"
scriptname = "/users/dribeir1/code/cod/src/cod/paper_figures/section2/gtex_geuvadis_roc_code.R"

uniqueCOPs = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/unique_COPs_mintissue5.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")
specificCOPs = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/specific_COPs_mintissue5_cutoff0.15.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")
prevalentCOPs = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/prevalent_COPs_mintissue5_cutoff0.15_0.5.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")
conservedCOPs = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/conserved_COPs_mintissue5_cutoff0.5.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")


########
# DATASET 1
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIDs = mergedData[significant == 1][pairID %in% uniqueCOPs$pairID]$nullId
mergedData = mergedData[nullId %in% nullIDs]

# Run and store
mergedData1 = mergedData
source(scriptname)
dataset1 = resultDF
dataset1Name = "1: Unique" # COPs"
dataset1$dataset = dataset1Name

########
# DATASET 2
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIDs = mergedData[significant == 1][pairID %in% specificCOPs$pairID]$nullId
mergedData = mergedData[nullId %in% nullIDs]

# Run and store
mergedData2 = mergedData
source(scriptname)
dataset2 = resultDF
dataset2Name = "2: Specific"# COPs"
dataset2$dataset = dataset2Name

########
# DATASET 3
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIDs = mergedData[significant == 1][pairID %in% prevalentCOPs$pairID]$nullId
mergedData = mergedData[nullId %in% nullIDs]

# Run and store
mergedData3 = mergedData
source(scriptname)
dataset3 = resultDF
dataset3Name = "3: Prevalent"# COPs"
dataset3$dataset = dataset3Name

########
# DATASET 4
########

# Chose file and processing
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nullIDs = mergedData[significant == 1][pairID %in% conservedCOPs$pairID]$nullId
mergedData = mergedData[nullId %in% nullIDs]

# Run and store
mergedData4 = mergedData
source(scriptname)
dataset4 = resultDF
dataset4Name = "4: Conserved"# COPs"
dataset4$dataset = dataset4Name

##########
# Plot all AUC results
##########
# Merge datasets
resultDF = rbind(dataset1,dataset2,dataset3,dataset4)

## 4 datasets
meansDF = data.table(aggregate(resultDF$auc, list(resultDF$model,resultDF$dataset), mean))
colnames(meansDF) = c("model","dataset","auc")
g1 = ggplot( resultDF, aes(x = model, y = auc, fill = dataset ))  +
  geom_boxplot( width = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", size = 1) +
  geom_text(data = meansDF[dataset == dataset1Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = 2.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset2Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = 1.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset3Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = -0.0, size = 4.5, fontface = "bold", color = "#525252" ) +
  geom_text(data = meansDF[dataset == dataset4Name], aes(x = model, y = auc-0.022, label = round(auc,2)), hjust = -1.1, size = 4.5, fontface = "bold", color = "#525252" ) +
  ylim(0.47,0.85) +
  ylab("AUC") +
  xlab("Models")
  ## GTEx LCL
  if (tissue == "Cells_EBV-transformed_lymphocytes"){
    g1 = g1 + scale_fill_brewer(palette="Set3", labels = c("Unique (N = 1429)","Specific (N = 1217)", "Prevalent (N = 968)","Conserved (N = 836)"))
  }
  if (tissue == "Lung"){
    g1 = g1 + scale_fill_brewer(palette="Set3", labels = c("Unique (N = 517)","Specific (N = 847)", "Prevalent (N = 1397)","Conserved (N = 1580)"))
  }
  if (tissue == "Muscle_Skeletal"){
    g1 = g1 + scale_fill_brewer(palette="Set3", labels = c("Unique (N = 1279)","Specific (N = 1606)", "Prevalent (N = 1395)","Conserved (N = 1064)"))
  }
  g1 = g1 + scale_x_discrete(labels = c("Distance","CTCF","Enhancer","TF","Expression","GO","eQTL","All together")) + 
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        axis.text.x = element_text(angle = 20, vjust=0.6),
        legend.position = c(0.15, 0.85), legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


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
g2 = ggplot( data = mergedData, aes(x = dataset, y = totalEnhancers, fill = as.factor(significant)) ) +
  # geom_violin( width = 0.6, size = 1) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 6, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 6, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 38, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 6 ) +
  ylim(0,40) +
  # ggtitle("Enhancers") +
  ylab("Total enhancers") +
  xlab("Dataset") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  scale_x_discrete(labels = c("Unique","Specific","Prevalent", "Conserved")) + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30),
        axis.text.x = element_text(angle = 20, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines')  
  )

##########
# TFs
##########
p1 = wilcox.test(mergedData[dataset == dataset1Name][significant == 1]$totalTFBS, mergedData[dataset == dataset1Name][significant == 0]$totalTFBS)$p.value
p2 = wilcox.test(mergedData[dataset == dataset2Name][significant == 1]$totalTFBS, mergedData[dataset == dataset2Name][significant == 0]$totalTFBS)$p.value
p3 = wilcox.test(mergedData[dataset == dataset3Name][significant == 1]$totalTFBS, mergedData[dataset == dataset3Name][significant == 0]$totalTFBS)$p.value
p4 = wilcox.test(mergedData[dataset == dataset4Name][significant == 1]$totalTFBS, mergedData[dataset == dataset4Name][significant == 0]$totalTFBS)$p.value
pvals = data.table(pvals = c(p1,p2,p3,p4), dataset = c(dataset1Name,dataset2Name,dataset3Name,dataset4Name))

meansDF = data.table(aggregate(mergedData$totalTFBS, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )
g3 = ggplot( data = mergedData, aes(x = dataset, y = totalTFBS, fill = as.factor(significant)) ) +
  # geom_violin( width = 0.6, size = 1) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 6, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 6, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 1150, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 6 ) +
  # ylim(0,40) +
  # ggtitle("Enhancers") +
  ylab("Total TFBS") +
  xlab("Dataset") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  scale_x_discrete(labels = c("Unique","Specific","Prevalent", "Conserved")) + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30),
        axis.text.x = element_text(angle = 20, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines')  
  )

lay <- rbind(c(1,1), c(2,3))
grid.arrange(g1,g2,g3, layout_matrix = lay)

dev.off()

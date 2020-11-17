#02-June-2020 # Diogo Ribeiro @ UNIL
# Script to perform molecular feature result comparisons

pdf("geuvadis_expression_control.pdf", 20, 20)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

########
# INPUT FILES
########
filename = "multi_features_distance_matched.bed"

tissue = "Lung"
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData1 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset1Name = "1: Lung"
tissue =  "Muscle_Skeletal"
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData2 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2Name = "2: Muscle Skeletal"
tissue = "Cells_EBV-transformed_lymphocytes"
inFile = paste("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/GTEx/",tissue,"/",filename,sep="")
mergedData3 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset3Name = "3: LCL (GTEx)"
inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData4 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset4Name = "4: LCL (Geuvadis)"

mergedData1 = mergedData1[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData2 = mergedData2[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData3 = mergedData3[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]
mergedData4 = mergedData4[,c("corr","distance","eGenes","eqtlSharing","totalTFBS","sharedTF","totalCTCF","invertedCTCF","sharedEnhancers","totalEnhancers","goSharing","diffExpr","significant")]

mergedData1$dataset = dataset1Name
mergedData2$dataset = dataset2Name
mergedData3$dataset = dataset3Name
mergedData4$dataset = dataset4Name

mergedData = rbind(mergedData1,mergedData2,mergedData3,mergedData4)

p1 = wilcox.test(mergedData[dataset == dataset1Name][significant == 1]$diffExpr, mergedData[dataset == dataset1Name][significant == 0]$diffExpr)$p.value
p2 = wilcox.test(mergedData[dataset == dataset2Name][significant == 1]$diffExpr, mergedData[dataset == dataset2Name][significant == 0]$diffExpr)$p.value
p3 = wilcox.test(mergedData[dataset == dataset3Name][significant == 1]$diffExpr, mergedData[dataset == dataset3Name][significant == 0]$diffExpr)$p.value
p4 = wilcox.test(mergedData[dataset == dataset4Name][significant == 1]$diffExpr, mergedData[dataset == dataset4Name][significant == 0]$diffExpr)$p.value
pvals = data.table(pvals = c(p1,p2,p3,p4), dataset = c(dataset1Name,dataset2Name,dataset3Name,dataset4Name))

meansDF = data.table(aggregate(mergedData$diffExpr, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )
g1 = ggplot( data = mergedData, aes(x = dataset, y = totalEnhancers, fill = as.factor(significant)) ) +
  # geom_violin( width = 0.6, size = 1) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 32, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 8 ) +
  ylim(c(0,35)) +
  ylab("Expression difference") +
  xlab("Dataset") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  scale_x_discrete(labels = c("Lung","Muscle Skeletal","LCL (GTEx)", "LCL (Geuvadis)")) + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30),
        axis.text.x = element_text(angle = 20, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")  
  )


########
# Other plots
########

inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/CODer_distance_and_variable_controlled_null.bed_positive_expr_metrics"
beforeData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/controlling_expression/CODer_distance_and_variable_controlled_null.bed_positive"
afterData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

### Small plots with mean Expr 
g2 = ggplot(beforeData, aes(x = log2(meanExpr), color = as.factor(significant)) ) +
  geom_density( aes(linetype = as.factor(significant)), alpha = 0.5, size = 3) +
  scale_color_manual( values = c("#d9d9d9","#66C2A5")) +
  ggtitle("Distance control only") +
  xlab("mean expression (log2)") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=24),
        axis.text.x = element_text(angle = 25, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")
)

g3 = ggplot(afterData, aes(x = log2(meanExpr), color = as.factor(significant)) ) +
  geom_density( aes(linetype = as.factor(significant)), alpha = 0.5, size = 3) +
  scale_color_manual( values = c("#d9d9d9","#66C2A5")) +
  ggtitle("Distance and expression control") +
  xlab("mean expression (log2)") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=24),
      axis.text.x = element_text(angle = 25, vjust=0.6),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
      legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm") 
)

# report how many cops
nrow(beforeData[significant == 1])
# report how many cops
nrow(afterData[significant == 1])

beforeData = beforeData[,.(diffExpr,diffCoef,significant)]
afterData = afterData[,.(diffExpr,diffCoef,significant)]
beforeData$dataset = "1: Distance control"
afterData$dataset = "2: Distance + Expression control"
mergedData = rbind(beforeData,afterData)

################
# Expr diff
################
p1 = wilcox.test(mergedData[dataset == "1: Distance control"][significant == 1]$diffExpr, mergedData[dataset == "1: Distance control"][significant == 0]$diffExpr)$p.value
p2 = wilcox.test(mergedData[dataset == "2: Distance + Expression control"][significant == 1]$diffExpr, mergedData[dataset == "2: Distance + Expression control"][significant == 0]$diffExpr)$p.value
pvals = data.table(pvals = c(p1,p2), dataset = c("1: Distance control","2: Distance + Expression control"))
meansDF = data.table(aggregate(mergedData$diffExpr, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )

g4 = ggplot( data = mergedData, aes(x = dataset, y = diffExpr, fill = as.factor(significant)) ) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 2.2, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 8 ) +
  # ggtitle("Gene pair expression level difference") +
  ylab("Expression difference") +
  xlab("") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=35),
        axis.text.x = element_text(vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")  )

################
# Coef var diff
################
p1 = wilcox.test(mergedData[dataset == "1: Distance control"][significant == 1]$diffCoef, mergedData[dataset == "1: Distance control"][significant == 0]$diffCoef)$p.value
p2 = wilcox.test(mergedData[dataset == "2: Distance + Expression control"][significant == 1]$diffCoef, mergedData[dataset == "2: Distance + Expression control"][significant == 0]$diffCoef)$p.value
pvals = data.table(pvals = c(p1,p2), dataset = c("1: Distance control","2: Distance + Expression control"))
meansDF = data.table(aggregate(mergedData$diffCoef, list(significant = mergedData$significant, dataset = mergedData$dataset), mean ) )

g5 = ggplot( data = mergedData, aes(x = dataset, y = diffCoef, fill = as.factor(significant)) ) +
  geom_boxplot( outlier.shape = NA, width = 0.5, size = 1) +
  geom_text(data = meansDF[significant == 1], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  geom_text(data = meansDF[significant == 0], aes(x = dataset, y = x, label = round(x,1)), size = 8, fontface = "bold", color = "#525252", nudge_x = -0.35) +
  geom_text(data = pvals, aes(x = dataset, y = 2.2, fill = "1", label = paste("P =",format(pvals, scientific = TRUE, digits = 1)) ), size = 8 ) +
  # ggtitle("Gene pair expression coefficient of variation difference") +
  ylab("Coef. variation difference") +
  xlab("") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=35),
        axis.text.x = element_text(vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")  )

pdf("/users/dribeir1/code/cod/src/cod/paper_figures/section2/geuvadis_expression_control.pdf", 26, 16)

lay <- rbind(c(1,1,1,1,2,2,2,3,3,3), c(4,4,4,4,4,5,5,5,5,5))
grid.arrange(g1,g2,g3,g4,g5, layout_matrix = lay)

dev.off()


#01-Jul-2020 Diogo Ribeiro @ UNIL
# Script to plot functional enrichment of shared and unshared eQTLs

pdf(paste("lead_eqtl_enrichment_main.pdf",sep=""),11,13)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)

inFile = "cod_analysis/geuvadis/functional_enrichments/eQTL/meta_results.out"
data = fread( inFile, stringsAsFactors = FALSE, header = F, sep=" ")

data$annotation = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[1])))$V1
data$category = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[2])))$V1

##############
# Geuvadis LCL Overlap against background
##############

g1 = ggplot() +
  geom_segment(data = data[category == "shared"], aes(x = V8, xend = V6, y = annotation, yend = annotation, color = annotation), position = position_nudge(y = 0.2), size = 1 ) +
  geom_point(data = data[category == "shared"], aes(x = V7, y = annotation, fill = annotation), color = "black", position = position_nudge(y = 0.2), size = 3, shape = 21) +
  geom_text(data = data[category == "shared"], aes(x = V7+0.12, y = annotation, color = annotation, label = paste0("P=",format.pval(V5, digits = 1), " OR=",format.pval(V7, digits = 1) ) ), position = position_nudge(y = 0.35), size = 4, fontface = "bold") +
  geom_segment(data = data[category == "unshared"], aes(x = V8, xend = V6, y = annotation, yend = annotation, color = annotation), position = position_nudge(y = -0.2), size = 1, alpha = 0.5 ) +
  geom_point(data = data[category == "unshared"], aes(x = V7, y = annotation, fill = annotation), color = "black", position = position_nudge(y = -0.2), size = 3, shape = 24, alpha = 0.5) +
  geom_text(data = data[category == "unshared"], aes(x = V7+0.12, y = annotation, color = annotation, label = paste0("P=",format.pval(V5, digits = 1), " OR=",format.pval(V7, digits = 1) )), position = position_nudge(y = -0.05), size = 4, fontface = "bold", alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = "dashed") +
  # ggtitle("Odds ratio against background") +
  xlab("Odds ratio") +
  xlim(c(0,4.25) ) +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  # scale_color_manual(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), 
        axis.title.y = element_blank(), legend.position = "None", 
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"))

g2 = ggplot() +
  geom_bar( data = data[category == "shared"], stat = "identity", aes(y = V1*100/V2, x = annotation, fill = annotation), position = position_nudge(x = 0.2), color = "black",  width = 0.4) +
  geom_bar( data = data[category == "unshared"], stat = "identity", aes(y = V1*100/V2, x = annotation, fill = annotation), position = position_nudge(x = -0.2), color = "black", width = 0.4, alpha = 0.5) +
  # ggtitle("Observed and expected (sd)") +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
    # scale_fill_manual(values = c("#7fc97f","#fdc086"), name = "eQTL") +
  ylab("% overlap") +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none", legend.key = element_rect(size = 1), legend.key.size = unit(1, 'lines'), legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.2,0.8,0.2,0), "cm"))

##############
# Geuvadis LCL
##############
mergedData = merge(data[category == "shared"], data[category == "unshared"], by = "annotation")

mergedData$odds = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$estimate )
mergedData$confmin = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$conf.int[1] )
mergedData$confmax = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$conf.int[2] )
mergedData$pval = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$p.value )

g3 = ggplot(mergedData, aes(x = odds, y = annotation, fill = annotation)) +
  geom_segment(data = mergedData, aes(x = confmin, xend = confmax, y = annotation, yend = annotation, color = annotation), size = 1 ) +
  geom_point(size = 4, shape = 21) +
  geom_text(data = mergedData, aes(x = odds, y = annotation, label = paste("OR=",round(odds,1), sep = "" ), color = annotation ), position = position_nudge(y = 0.25), size = 4, fontface = "bold") +
  geom_text(data = mergedData, aes(x = odds, y = annotation, label = paste("P=",format.pval(pval, digits = 1), sep = "" ), color = annotation ), position = position_nudge(y = -0.25), size = 4, fontface = "bold") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # ggtitle("Shared lead eQTLs vs not shared") +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  xlim(c(0,min(max(mergedData$confmax), 3.5)) ) +
  xlab("Odds ratio") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "None", panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), , plot.margin = unit(c(0.2,0.8,0.2,0), "cm"))

##############
# GTEx LCL
##############
inFile = "cod_analysis/GTEx/functional_enrichments/LCLs/eQTL/meta_results.out"
data = fread( inFile, stringsAsFactors = FALSE, header = F, sep=" ")
data$annotation = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[1])))$V1
data$category = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[2])))$V1

mergedData = merge(data[category == "shared"], data[category == "unshared"], by = "annotation")

mergedData$odds = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$estimate )
mergedData$confmin = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$conf.int[1] )
mergedData$confmax = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$conf.int[2] )
mergedData$pval = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3])-as.numeric(x[2]),as.numeric(x[12]),as.numeric(x[13])-as.numeric(x[12]) ),nrow = 2), conf.level = 0.95)$p.value )

g4 = ggplot(mergedData, aes(x = odds, y = annotation, fill = annotation)) +
  geom_segment(data = mergedData, aes(x = confmin, xend = confmax, y = annotation, yend = annotation, color = annotation), size = 1 ) +
  geom_point(size = 4, shape = 21) +
  geom_text(data = mergedData, aes(x = odds, y = annotation, label = paste("OR=",round(odds,1), sep = "" ), color = annotation ), position = position_nudge(y = 0.25), size = 4, fontface = "bold") +
  geom_text(data = mergedData, aes(x = odds, y = annotation, label = paste("P=",format.pval(pval, digits = 1), sep = "" ), color = annotation ), position = position_nudge(y = -0.25), size = 4, fontface = "bold") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # ggtitle("Shared lead eQTLs vs not shared") +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  xlim(c(0,min(max(mergedData$confmax), 3.5)) ) +
  xlab("Odds ratio") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "None", panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank())

#####
statsFile = "cod_identification/GTEx/coder_stats.txt"
copFile = "zcat cod_identification/GTEx/CODer_final_dataset_cops_merged.bed.gz"
sharedEQTLFile = "cod_analysis/GTEx/eqtl_sharing/meta_results_shared.eqtl"
eQTLFile = "eqtl/GTEx/eQTLs/permutation_pass/results/meta_results.eqtl"
colorFile = "raw_input/GTEx/v8/other/color_code_computer.txt"

statsData = fread( statsFile, stringsAsFactors = FALSE, header = FALSE, sep="\t")
copData = fread( copFile, stringsAsFactors = FALSE, header = T, sep="\t")
sharedEQTLData = fread( sharedEQTLFile, stringsAsFactors = FALSE, header = F, sep="\t")
eQTLData = fread( eQTLFile, stringsAsFactors = FALSE, header = F, sep=" ")
colorCode = fread( colorFile, stringsAsFactors = FALSE, header = FALSE, sep="\t")

paste("Total unique lead eQTLs:",length(unique(eQTLData$V8)) )
paste("Total unique eGenes:",length(unique(eQTLData$V1)) )
paste("Total unique shared lead eQTLs:",length(unique(sharedEQTLData$V1)) )
paste("Total unique genes:",length(unique(sharedEQTLData$V2)) )

colnames(statsData) = c("Tissue","Metric","Value")
copData = unique(copData)

eGenesData = unique(eQTLData[,.(V1,V21)])
eGenesPerTissue = data.table(table(eGenesData$V21))
colnames(eGenesPerTissue) = c("Tissue","eGenes")
mergedData = merge(eGenesPerTissue, statsData[Metric == "Total Genes tested"], by = "Tissue")
mergedData$perc = mergedData$eGenes * 100 / mergedData$Value
sharedEQTLAll = unique(sharedEQTLData[,.(V1,V2,V4)])
eGenesData = unique(sharedEQTLAll[,.(V2,V4)])
eGenesPerTissue = data.table(table(eGenesData$V4))
colnames(eGenesPerTissue) = c("Tissue","shared_eGenes")
mergedData = merge(eGenesPerTissue, mergedData, by = "Tissue")
mergedData$perc = mergedData$shared_eGenes * 100 / mergedData$eGenes

##################
# Plot % COPs that have a shared eQTL per tissue
##################

sharedEQTLData$V2 = data.table(unlist(lapply(sharedEQTLData$V2, function(x) gsub("[.].*","",x) )))
sharedEQTLData$V3 = data.table(unlist(lapply(sharedEQTLData$V3, function(x) gsub("[.].*","",x) )))
p1 = unite(sharedEQTLData[V2 <= V3], pairID, c(V2, V3), remove=FALSE, sep = "|")
p2 = unite(sharedEQTLData[V2 > V3], pairID, c(V3, V2), remove=FALSE, sep = "|")
sharedEQTLData = rbind(p1,p2)
# nrow(copData[pairID %in% sharedEQTLData$pairID]) / nrow(copData)
#85460/135662 # 63% of COPs share an eQTL (across tissues) 
d1 = merge(sharedEQTLData, copData, by.x = c("pairID","V4"), by.y = c("pairID","tissue"), all.x = T)
shareCopPerTissue = merge(data.table(table(d1$V4)), data.table(table(copData$tissue)), by = "V1")
colnames(shareCopPerTissue) = c("Tissue","sharedCOP","allCOP")
shareCopPerTissue$perc = shareCopPerTissue$sharedCOP * 100 / shareCopPerTissue$allCOP
summary(shareCopPerTissue$perc)

shareCopPerTissue$tissueName = data.table(unlist(lapply(shareCopPerTissue$Tissue, function(x) unlist(gsub("_"," ",x[1]))[1])))$V1

g5 = ggplot( shareCopPerTissue, aes(x = reorder(tissueName, perc), y = perc, fill = tissueName, label = sharedCOP ) ) +
  geom_bar( stat = "identity", color = "black", show.legend = FALSE, alpha = 0.5) +
  geom_text(hjust = 1, size = 4.5, angle = 90) +
  ylab("% COPs with shared lead eQTL") +
  xlab("Tissue") +
  scale_fill_manual(values = colorCode$V2) +
  # coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(hjust = 1, angle = 45), axis.title.x = element_text(size = 20))

lay <- rbind(c(1,1,1,2,3,3,4,4), c(5,5,5,5,5,5,5,5))
grid.arrange(g1,g2,g3,g4,g5, layout_matrix = lay)

dev.off()


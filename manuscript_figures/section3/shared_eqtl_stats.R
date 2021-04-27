#18-Aug-2020 Diogo Ribeiro @ UNIL
# Script to plot shared lead eQTLs across GTEx tissues

pdf("shared_eqtl_stats.pdf",13,12)

library(data.table)
library(ggplot2)
library(tidyr)

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

##################
# Plot % eGenes per tissue
##################

eGenesData = unique(eQTLData[,.(V1,V21)])
eGenesPerTissue = data.table(table(eGenesData$V21))
colnames(eGenesPerTissue) = c("Tissue","eGenes")

mergedData = merge(eGenesPerTissue, statsData[Metric == "Total Genes tested"], by = "Tissue")
mergedData$perc = mergedData$eGenes * 100 / mergedData$Value

median(mergedData$perc)
  
# ggplot( mergedData, aes(x = reorder(Tissue, perc), y = perc, fill = Tissue, label = eGenes ) ) +
#   geom_bar( stat = "identity", color = "black", show.legend = FALSE, alpha = 0.5) +
#   geom_text(hjust = 1, size = 6) +
#   ylab("% genes") +
#   xlab("Tissue") +
#   ggtitle("% genes as eGenes per tissue") +
#   scale_fill_manual(values = colorCode$V2) +
#   coord_flip() +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1))

##################
# Plot % eGenes that have a shared eQTL per tissue
##################

sharedEQTLAll = unique(sharedEQTLData[,.(V1,V2,V4)])
eGenesData = unique(sharedEQTLAll[,.(V2,V4)])
eGenesPerTissue = data.table(table(eGenesData$V4))
colnames(eGenesPerTissue) = c("Tissue","shared_eGenes")

mergedData = merge(eGenesPerTissue, mergedData, by = "Tissue")
mergedData$perc = mergedData$shared_eGenes * 100 / mergedData$eGenes

summary(mergedData$perc)

# ggplot( mergedData, aes(x = reorder(Tissue, perc), y = perc, fill = Tissue, label = shared_eGenes ) ) +
#   geom_bar( stat = "identity", color = "black", show.legend = FALSE, alpha = 0.5) +
#   geom_text(hjust = 1, size = 6) +
#   ylab("% eGenes") +
#   xlab("Tissue") +
#   ggtitle("% eGenes with shared lead eQTL per tissue") +
#   scale_fill_manual(values = colorCode$V2) +
#   coord_flip() +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1))

##################
# Plot % COPs that have a shared eQTL per tissue
##################

sharedEQTLData$V2 = data.table(unlist(lapply(sharedEQTLData$V2, function(x) gsub("[.].*","",x) )))
sharedEQTLData$V3 = data.table(unlist(lapply(sharedEQTLData$V3, function(x) gsub("[.].*","",x) )))

p1 = unite(sharedEQTLData[V2 <= V3], pairID, c(V2, V3), remove=FALSE, sep = "|")
p2 = unite(sharedEQTLData[V2 > V3], pairID, c(V3, V2), remove=FALSE, sep = "|")
sharedEQTLData = rbind(p1,p2)

nrow(copData[pairID %in% sharedEQTLData$pairID]) / nrow(copData)
#85460/135662 # 63% of COPs share an eQTL (across tissues) 

d1 = merge(sharedEQTLData, copData, by.x = c("pairID","V4"), by.y = c("pairID","tissue"), all.x = T)

shareCopPerTissue = merge(data.table(table(d1$V4)), data.table(table(copData$tissue)), by = "V1")
colnames(shareCopPerTissue) = c("Tissue","sharedCOP","allCOP")
shareCopPerTissue$perc = shareCopPerTissue$sharedCOP * 100 / shareCopPerTissue$allCOP
summary(shareCopPerTissue$perc)

# ggplot( shareCopPerTissue, aes(x = reorder(Tissue, perc), y = perc, fill = Tissue, label = sharedCOP ) ) +
#   geom_bar( stat = "identity", color = "black", show.legend = FALSE, alpha = 0.5) +
#   geom_text(hjust = 1, size = 6) +
#   ylab("% COPs") +
#   xlab("Tissue") +
#   ggtitle("% COPs with shared lead eQTL per tissue") +
#   scale_fill_manual(values = colorCode$V2) +
#   coord_flip() +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1))


##################
# COPs with eQTL sharing per tissue sample size
##################
mergedData = unique(merge(shareCopPerTissue,statsData[Metric == "Total samples"],by="Tissue"))

d = lm(mergedData$perc ~ mergedData$Value)

correlation = cor.test(mergedData$perc, mergedData$Value, method = "spearman")
correlationText = paste("Spearman R:",round(correlation$estimate,2), "P-value:",format.pval(correlation$p.value,2), sep = " ")
ggplot( mergedData, aes(x = Value, y = perc, color = Tissue, fill = Tissue ) ) +
  geom_abline(slope = d$coefficients[2], intercept = d$coefficients[1], color = "#de2d26") +
  geom_point( show.legend = FALSE, alpha = 0.9, size = 3, shape = 21, color = "black") +
  geom_text( aes(label = Tissue), size = 7, vjust = 1.5, check_overlap = TRUE, show.legend = FALSE) +
  annotate("text", x = 420, y = 10, label = correlationText, hjust = 0, vjust =1, size = 8.5  ) +
  # ggtitle("% COPs with eQTL sharing per tissue sample size") +
  ylab("% COPs with eQTL sharing") +
  xlab("Tissue sample size") +
  scale_fill_manual(values = colorCode$V2) +
  scale_color_manual(values = colorCode$V2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), axis.text.x = element_text(angle = 0, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2))

dev.off()


#23-June-2020 # Diogo Ribeiro @ UNIL
# Script to perform a logistic regression between genes being co-expressed and several molecular features

pdf("eqtl_sharing.pdf",10,14)

library(data.table)
library(ggplot2)
library(tidyr)
library(grid)
library(gridExtra)

################
# eGenes only
################
inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

mergedData$group = "COP"
mergedData[significant == 0]$group = "Non-COP"
mergedData[eqtlSharing == 2]$eqtlSharing = 1

mergedData = mergedData[eGenes == 2]
d = data.table(table(mergedData$group,mergedData$eqtlSharing))
colnames(d) = c("group","eqtlSharing","N")
d$fill = c(2,0,3,1)
perc1 = d[group == "COP"][eqtlSharing == 1]$N * 100.0 / (d[group == "COP"][eqtlSharing == 1]$N + d[group == "COP"][eqtlSharing == 0]$N)
perc2 = d[group == "Non-COP"][eqtlSharing == 1]$N * 100.0 / (d[group == "Non-COP"][eqtlSharing == 1]$N + d[group == "Non-COP"][eqtlSharing == 0]$N)

g1 = ggplot( data = d, aes(x = group, y = N, fill = as.factor(fill) ) ) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c("#d9d9d9","#969696","#ccece6","#66C2A5"), label = c("Non-COP | Not Sharing","Non-COP | Sharing", "COP | Not Sharing", "COP | Sharing")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "COP", y = 1000, size = 8) + 
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "COP", y = 2200, size = 8) + 
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "Non-COP", y = 140, size = 8) + 
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "Non-COP", y = 700, size = 8) + 
  ylab("# gene pairs") +
  ggtitle("eQTL sharing | eGene pairs") +
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white",
                                        size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 5), legend.title = element_blank(),
        legend.key.size = unit(1.7, 'lines')  )

################
# All genes pairs, also trans
################

inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
mergedData$group = "COP"
mergedData[significant == 0]$group = "Non-COP"
mergedData[eqtlSharing == 2]$eqtlSharing = 1
mergedData = mergedData[,.(group, eqtlSharing, pairID)]

inFile = "cod_analysis/geuvadis/eQTLs/eqtl_sharing/trans_vs_cis/4_groups_eqtl_sharing.bed"
mergedData2 = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
mergedData2 = mergedData2[group == "Trans-COP"]
mergedData2[eqtlSharing == 2]$eqtlSharing = 1
mergedData2 = unique(mergedData2[,.(group,eqtlSharing,pairID)])

mergedData = rbind(mergedData,mergedData2)

d = data.table(table(mergedData$group,mergedData$eqtlSharing))
colnames(d) = c("group","eqtlSharing","N")
d = d[order(group)]
d$fill = c(1,2,3,4,5,6)
perc1 = d[group == "COP"][eqtlSharing == 1]$N * 100.0 / (d[group == "COP"][eqtlSharing == 1]$N + d[group == "COP"][eqtlSharing == 0]$N)
perc2 = d[group == "Non-COP"][eqtlSharing == 1]$N * 100.0 / (d[group == "Non-COP"][eqtlSharing == 1]$N + d[group == "Non-COP"][eqtlSharing == 0]$N)
perc3 = d[group == "Trans-COP"][eqtlSharing == 1]$N * 100.0 / (d[group == "Trans-COP"][eqtlSharing == 1]$N + d[group == "Trans-COP"][eqtlSharing == 0]$N)

g2 = ggplot( data = d, aes(x = group, y = N, fill = as.factor(fill) ) ) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c("#ccece6","#66C2A5","#d9d9d9","#969696","#ffeda0", "#feb24c"), 
                     label = c("COP | No Sharing", "COP | Sharing", "Non-COP | No Sharing","Non-COP | Sharing", "Trans-COP | No Sharing", "Trans-COP | Sharing")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "COP", y = 1500, size = 8) + 
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "COP", y = 4800, size = 8) + 
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "Non-COP", y = 220, size = 8) + 
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "Non-COP", y = 3500, size = 8) + 
  annotate(geom = "text", label = paste0(round(perc3,1),"%"), x = "Trans-COP", y = 220, size = 8) + 
  annotate(geom = "text", label = paste0(round(100-perc3,1),"%"), x = "Trans-COP", y = 3500, size = 8) + 
  ylab("# gene pairs") +
  ggtitle("eQTL sharing | all gene pairs") +
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white",
                                        size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 5), legend.title = element_blank(),
        legend.key.size = unit(1.7, 'lines')  )

grid.arrange(g2,g1, nrow = 2)


dev.off()
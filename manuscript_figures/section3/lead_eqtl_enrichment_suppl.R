#01-Jul-2020 Diogo Ribeiro @ UNIL
# Script to plot functional enrichment of shared and unshared eQTLs

args = commandArgs(trailingOnly=TRUE)

inFile = args[1]
outTag = args[2]

# inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/geuvadis/functional_enrichments/eQTL/roadmap/meta_results.out"
# inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/functional_enrichments/LCLs/eQTL/roadmap/meta_results.out"
# inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/functional_enrichments/Lung/eQTL/roadmap/meta_results.out"
# inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/functional_enrichments/Muscle/eQTL/roadmap/meta_results.out"
# outTag = "gtex_lcl"

pdf(paste("/users/dribeir1/code/cod/src/cod/paper_figures/section3/lead_eqtl_enrichment_suppl_",outTag,".pdf",sep=""),14,10)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)

data = fread( inFile, stringsAsFactors = FALSE, header = F, sep=" ")

data$annotation = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[1])))$V1
data$category = data.table(unlist(lapply(data$V9, function(x) unlist(strsplit(x,".", fixed = T))[2])))$V1

##############
# Geuvadis LCL Overlap against background
##############

unique(data$annotation)

g1 = ggplot() +
  geom_segment(data = data[category == "shared"], aes(x = V8, xend = V6, y = annotation, yend = annotation, color = annotation), position = position_nudge(y = 0.2), size = 1 ) +
  geom_point(data = data[category == "shared"], aes(x = V7, y = annotation, fill = annotation), color = "black", position = position_nudge(y = 0.2), size = 3, shape = 21) +
  geom_text(data = data[category == "shared"], aes(x = V7, y = annotation, color = annotation, label = paste0("P=",format.pval(V5, digits = 1), " OR=",format.pval(V7, digits = 1) ) ), position = position_nudge(y = 0.35), size = 4, fontface = "bold") +
  geom_segment(data = data[category == "unshared"], aes(x = V8, xend = V6, y = annotation, yend = annotation, color = annotation), position = position_nudge(y = -0.2), size = 1, alpha = 0.5 ) +
  geom_point(data = data[category == "unshared"], aes(x = V7, y = annotation, fill = annotation), color = "black", position = position_nudge(y = -0.2), size = 3, shape = 24, alpha = 0.5) +
  geom_text(data = data[category == "unshared"], aes(x = V7, y = annotation, color = annotation, label = paste0("P=",format.pval(V5, digits = 1), " OR=",format.pval(V7, digits = 1)  )), position = position_nudge(y = -0.05), size = 4, fontface = "bold", alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = "dashed") +
  scale_y_discrete(limits = c("15_Quies","14_ReprPCWk","13_ReprPC","12_EnhBiv","11_BivFlnk","10_TssBiv","9_Het","8_ZNF","7_Enh","6_EnhG","5_TxWk","4_Tx","3_TxFlnk","2_TssAFlnk","1_TssA")) +
  xlab("Odds ratio") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), 
        axis.title.y = element_blank(), legend.position = "None", 
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"))

g2 = ggplot() +
  geom_bar( data = data[category == "shared"], stat = "identity", aes(y = V1*100/V2, x = annotation, fill = annotation), position = position_nudge(x = 0.2), color = "black",  width = 0.4) +
  geom_bar( data = data[category == "unshared"], stat = "identity", aes(y = V1*100/V2, x = annotation, fill = annotation), position = position_nudge(x = -0.2), color = "black", width = 0.4, alpha = 0.5) +
  ylab("% overlap") +
  scale_x_discrete(limits = c("15_Quies","14_ReprPCWk","13_ReprPC","12_EnhBiv","11_BivFlnk","10_TssBiv","9_Het","8_ZNF","7_Enh","6_EnhG","5_TxWk","4_Tx","3_TxFlnk","2_TssAFlnk","1_TssA")) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none", legend.key = element_rect(size = 1), legend.key.size = unit(1, 'lines'), legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.2,1,0.2,0), "cm"))

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
  scale_y_discrete(limits = c("15_Quies","14_ReprPCWk","13_ReprPC","12_EnhBiv","11_BivFlnk","10_TssBiv","9_Het","8_ZNF","7_Enh","6_EnhG","5_TxWk","4_Tx","3_TxFlnk","2_TssAFlnk","1_TssA")) +
  xlim(c(0,min(max(mergedData$confmax), 3.5)) ) +
  xlab("Odds ratio") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "None", panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank())


lay <- rbind(c(1,1,1,2,3,3))
grid.arrange(g1,g2,g3, layout_matrix = lay)

dev.off()


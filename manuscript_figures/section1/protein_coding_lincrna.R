#21-May-2020 Diogo Ribeiro @ UNIL
# Plotting protein coding vs lincRNA COPs

pdf("protein_coding_lincrna.pdf",18,9)

library(data.table)
require(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

baseFolder = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset"
resFile = paste(baseFolder,"/CODer_raw_results.bed",sep="")
copFile = paste(baseFolder,"/CODer_cod_identification_cops.bed",sep="")

resData = fread( resFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
copData = fread( copFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# get only central phenotype info
centralData = unique(resData[,.(centralPhenotype,centralInfo)])
centralData$inCops = 0
listGenesInCops = unique(c(copData$centralPhenotype, copData$cisPheno))
centralData[centralData$centralPhenotype %in% listGenesInCops]$inCops = 1

##############
# Bar plot between coding and lincRNA
##############

coding = centralData[grepl("protein_coding", centralInfo)]
lincrna = centralData[grepl("lincRNA",centralInfo)]
mergedData = rbind(data.table(table(coding$inCops)), data.table(table(lincrna$inCops)))
mergedData$gene_type = c("Protein coding","Protein coding","LincRNA","LincRNA")
colnames(mergedData) = c("COP","N","gene_type")

mergedData$yvalue = c(mergedData[COP == 0][gene_type == "Protein coding"]$N / 2 + mergedData[COP == 1][gene_type == "Protein coding"]$N, 
                      mergedData[COP == 1][gene_type == "Protein coding"]$N / 2, 
                      mergedData[COP == 0][gene_type == "LincRNA"]$N / 2 + mergedData[COP == 1][gene_type == "LincRNA"]$N, 
                      mergedData[COP == 1][gene_type == "LincRNA"]$N / 2)

mergedData$total = c(mergedData[COP == 0][gene_type == "Protein coding"]$N + mergedData[COP == 1][gene_type == "Protein coding"]$N,
                     mergedData[COP == 0][gene_type == "Protein coding"]$N + mergedData[COP == 1][gene_type == "Protein coding"]$N,
                     mergedData[COP == 0][gene_type == "LincRNA"]$N + mergedData[COP == 1][gene_type == "LincRNA"]$N,
                     mergedData[COP == 0][gene_type == "LincRNA"]$N + mergedData[COP == 1][gene_type == "LincRNA"]$N)

g1 = ggplot( mergedData, aes(x = gene_type, y = N, fill = COP ))  +
  geom_bar( stat = "identity", width = 0.7, size = 2, color = "black") +
  geom_text(aes(x = gene_type, y = yvalue, label = paste(N," (",round(N * 100 / total,2),"%)", sep = "") ), size = 7, fontface = "bold", color = "white" ) +
  # ggtitle("COPs in protein coding and lincRNAs") + 
  ylab("# Genes") +
  xlab("Gene type") +
  # coord_flip() +
  # scale_fill_brewer(palette="Set1") +
  scale_fill_manual(values=c("#377eb8", "#e41a1c"), 
                    labels=c("Non-COP", "COP")) + 
  theme_minimal() +
  theme_linedraw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
                           axis.text.x = element_text(angle = 0, vjust=0.6),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           axis.line = element_line(colour = "black", size = 1),
                           panel.border = element_rect(colour = "black", fill=NA, size=1), 
                           legend.title=element_blank(), legend.position = c(0.25,0.85), legend.text = element_text(size = 24))

##############
# Pie plot between coding and lincRNA co-expression pairs
##############
pieData = data.table(gene_pair_type = c("coding-coding","coding-lincRNA","lincRNA-lincRNA"), N = c(nrow(copData[grepl("protein_coding", centralInfo)][grepl("protein_coding", cisInfo)]), nrow(copData[grepl("protein_coding", centralInfo)][grepl("lincRNA", cisInfo)]) + nrow(copData[grepl("lincRNA", centralInfo)][grepl("protein_coding", cisInfo)]), nrow(copData[grepl("lincRNA", centralInfo)][grepl("lincRNA", cisInfo)])) )
tot = sum(pieData$N)
pieData$perc = round(pieData$N * 100.0 / tot,2)

pieData = pieData %>%
  arrange(desc(gene_pair_type)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)

g2 = ggplot( pieData, aes(x = "", y = perc, fill = gene_pair_type))  +
  geom_bar( stat = "identity", width = 0.7, size = 2, color = "black") +
  geom_text(aes(y = lab.ypos, 
                label = paste(N, "(", perc, "%)", sep = "")), size=7, color = "white", fontface = "bold", nudge_y = 5) +
  # ggtitle("COP gene type proportion") + 
  coord_polar("y", start=0) +
  xlab("") +
  scale_fill_brewer(palette="Set2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.ticks = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(),
                           legend.title=element_blank(), legend.position = c(0.67,0.44), legend.text = element_text(size = 24))

grid.arrange(g1,g2,nrow = 1)

dev.off()


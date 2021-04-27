#21-May-2020 Diogo Ribeiro @ UNIL
# Plotting protein coding vs lincRNA COPs

pdf("protein_coding_lincrna.pdf",18,9)

library(data.table)
require(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

baseFolder = "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset"
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

###################
# Gene adjacency
###################
allData = fread( "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_raw_results.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")
inFile = "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_distance_controlled_null.bed"
copData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
table(copData$significant)
allData$tss1 = allData$centralStart
allData[centralStrand == "-"]$tss1 = allData[centralStrand == "-"]$centralEnd
allData$tss2 = allData$cisStart
allData[cisStrand == "-"]$tss2 = allData[cisStrand == "-"]$cisEnd
d1 = data.table(gene = allData$centralPhenotype, tss = allData$tss1, chr = allData$`#chr`)
d2 = data.table(gene = allData$cisPheno, tss = allData$tss2, chr = allData$`#chr`)
allGenes = unique(rbind(d1,d2))
allGenes = allGenes[order(chr,tss)]
allGenes$idx = seq(1,nrow(allGenes))
copData = unique(copData[,.(centralPhenotype,cisPheno,pairID,distance,significant,nullId)])
allData = unique(allData[,.(centralPhenotype,cisPheno,pairID,distance)])
allGenes$tss = NULL
allGenes$chr = NULL
colnames(allGenes) = c("centralPhenotype","centralIdx")
d1 = merge(allData,allGenes, by = "centralPhenotype")
colnames(allGenes) = c("cisPheno","cisIdx")
mergedData = merge(d1,allGenes, by = "cisPheno")
mergedData$idxDiff = abs(mergedData$centralIdx - mergedData$cisIdx)
mergedData = unique(mergedData[,.(pairID,idxDiff)])
mergedData$significant = -1
# mergedData[pairID %in% copData[significant == 0]$pairID]$significant = 0
mergedData[pairID %in% copData[significant == 1]$pairID]$significant = 1
table(mergedData$significant)
summary(mergedData$idxDiff)
summary(mergedData[significant == 1]$idxDiff)
summary(mergedData[significant == 0]$idxDiff)
length(unique(mergedData[significant == 1][idxDiff == 1]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == 1][idxDiff == 2]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == 1][idxDiff > 2]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == -1][idxDiff == 1]$pairID))*100/nrow(mergedData[significant == -1])
length(unique(mergedData[significant == -1][idxDiff == 2]$pairID))*100/nrow(mergedData[significant == -1])
length(unique(mergedData[significant == -1][idxDiff > 2]$pairID))*100/nrow(mergedData[significant == -1])

g3 = ggplot(mergedData[significant != 0], aes(x = idxDiff, color = as.factor(significant), fill = as.factor(significant) )) +
  geom_density( alpha = 0.5, size = 1, adjust = 2) +
  scale_color_manual( values = c("#fec44f","#66C2A5"), label = c("Not co-expressed","COPs"), name = "Gene pairs") +
  scale_fill_manual( values = c("#fec44f","#66C2A5"), label = c("Not co-expressed","COPs"), name = "Gene pairs") +
  xlab("# genes in the middle") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=20), aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")
  )

## Same but for a single strand
allData = fread( "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_raw_results.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")
inFile = "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_distance_controlled_null.bed"
copData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

allData = allData[centralStrand == "+"][cisStrand == "+"]

allData$tss1 = allData$centralStart
allData[centralStrand == "-"]$tss1 = allData[centralStrand == "-"]$centralEnd
allData$tss2 = allData$cisStart
allData[cisStrand == "-"]$tss2 = allData[cisStrand == "-"]$cisEnd
d1 = data.table(gene = allData$centralPhenotype, tss = allData$tss1, chr = allData$`#chr`)
d2 = data.table(gene = allData$cisPheno, tss = allData$tss2, chr = allData$`#chr`)
allGenes = unique(rbind(d1,d2))
allGenes = allGenes[order(chr,tss)]
allGenes$idx = seq(1,nrow(allGenes))
copData = unique(copData[,.(centralPhenotype,cisPheno,pairID,distance,significant,nullId)])
allData = unique(allData[,.(centralPhenotype,cisPheno,pairID,distance)])
allGenes$tss = NULL
allGenes$chr = NULL
colnames(allGenes) = c("centralPhenotype","centralIdx")
d1 = merge(allData,allGenes, by = "centralPhenotype")
colnames(allGenes) = c("cisPheno","cisIdx")
mergedData = merge(d1,allGenes, by = "cisPheno")
mergedData$idxDiff = abs(mergedData$centralIdx - mergedData$cisIdx)
mergedData = unique(mergedData[,.(pairID,idxDiff)])
mergedData$significant = -1
# mergedData[pairID %in% copData[significant == 0]$pairID]$significant = 0
mergedData[pairID %in% copData[significant == 1]$pairID]$significant = 1
table(mergedData$significant)
summary(mergedData$idxDiff)
summary(mergedData[significant == 1]$idxDiff)
summary(mergedData[significant == 0]$idxDiff)
length(unique(mergedData[significant == 1][idxDiff == 1]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == 1][idxDiff == 2]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == 1][idxDiff > 2]$pairID))*100/nrow(mergedData[significant == 1])
length(unique(mergedData[significant == -1][idxDiff == 1]$pairID))*100/nrow(mergedData[significant == -1])
length(unique(mergedData[significant == -1][idxDiff == 2]$pairID))*100/nrow(mergedData[significant == -1])
length(unique(mergedData[significant == -1][idxDiff > 2]$pairID))*100/nrow(mergedData[significant == -1])

g4 = ggplot(mergedData[significant != 0], aes(x = idxDiff, color = as.factor(significant), fill = as.factor(significant) )) +
  geom_density( alpha = 0.5, size = 1, adjust = 2) +
  scale_color_manual( values = c("#fec44f","#66C2A5"), label = c("Not co-expressed","COPs"), name = "Gene pairs") +
  scale_fill_manual( values = c("#fec44f","#66C2A5"), label = c("Not co-expressed","COPs"), name = "Gene pairs") +
  xlab("# genes in the middle") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=20), aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 6),legend.key.size = unit(2, 'lines'), plot.margin = unit(c(1,1,1,1), "cm")
  )


layout = rbind( c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))
grid.arrange(g3,g4,nrow = 1, layout_matrix = layout)


dev.off()


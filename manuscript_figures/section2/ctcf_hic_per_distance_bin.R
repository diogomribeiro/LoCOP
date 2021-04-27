#23-09-2019 # Diogo Ribeiro @ UNIL
# Script to perform a logistic regression between genes being co-expressed and several molecular features

pdf("ctcf_hic_per_distance_bin.pdf",14,14)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

# seed for the sampling for training/testing dataset
set.seed(666)

# number of digits after comma
options("scipen"=100, "digits"=2)

########
# Load/process data
########

# wantedCorrSign = "+"
inFile = "cod_analysis/multiple_features/geuvadis/multi_features_distance_matched.bed"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

###############
# Create distance bins
###############

# Get results for bins of distance
step = 100000 #bp
distanceBinned = data.frame()
for (i in seq(1,1000000,step)){  
  filtDF = mergedData[mergedData$distance >= i]
  tempDF = filtDF[filtDF$distance < i+step]
  tempDF$bin = (i-1.0)/1000 # paste( (i-1)/1000,"-",(i+step-1)/1000,sep="")
  distanceBinned = rbind(distanceBinned, tempDF)
}
distanceBinned$bin = as.factor(distanceBinned$bin)


################
# Inverted CTCF plot
################
distanceBinnedSizes = as.data.table(table(distanceBinned[,c("bin","significant")]))
# add the % mean increase between significant and non-significant
m = data.table(aggregate(distanceBinned$invertedCTCF, list(significant = distanceBinned$significant, Bin = distanceBinned$bin), mean ) )
distanceBinMeans = round((m[m$significant == 1]$x * 100.0 / m[m$significant == 0]$x) - 100, 0)
v = c()
for (i in seq(1,nrow(distanceBinnedSizes) ) ){
  if (i <= length(distanceBinMeans)){ v = c(v, paste( as.factor(distanceBinMeans[i]), "%", sep = "" ) ) }
  else{ v = c(v,"0")  }
}
distanceBinnedSizes$means = v

g1 = ggplot( data = distanceBinned, aes(x = bin, y = invertedCTCF, fill = as.factor(significant)) ) +
  geom_boxplot( outlier.shape = NA, width = 0.9, size = 1) +
  geom_text(data = distanceBinnedSizes[distanceBinnedSizes$significant == 0], aes(bin, -1, label = means), size = 5) +
  ylim(-1,25) +
  ggtitle("Inverted CTCF motifs between gene TSSs") +
  ylab("Inverted CTCF motifs") +
  xlab("Gene TSS distance bin (Kb)") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.title=element_blank(), legend.key = element_rect(size = 6), 
        legend.key.size = unit(2, 'lines'), legend.position = c(0.3,0.85), aspect.ratio = 1  )


################
# Total CTCF plot
################
distanceBinnedSizes = as.data.table(table(distanceBinned[,c("bin","significant")]))
# add the % mean increase between significant and non-significant
m = data.table(aggregate(distanceBinned$totalCTCF, list(significant = distanceBinned$significant, Bin = distanceBinned$bin), mean ) )
distanceBinMeans = round((m[m$significant == 1]$x * 100.0 / m[m$significant == 0]$x) - 100, 0)
v = c()
for (i in seq(1,nrow(distanceBinnedSizes) ) ){
  if (i <= length(distanceBinMeans)){ v = c(v, paste( as.factor(distanceBinMeans[i]), "%", sep = "" ) ) }
  else{ v = c(v,"0")  }
}
distanceBinnedSizes$means = v

g2 = ggplot( data = distanceBinned, aes(x = bin, y = totalCTCF, fill = as.factor(significant)) ) +
  geom_boxplot( outlier.shape = NA, width = 0.9, size = 1) +
  geom_text(data = distanceBinnedSizes[distanceBinnedSizes$significant == 0], aes(bin, -10, label = means), size = 5) +
  ylim(-10,350) +
  ggtitle("Total CTCF sites between gene TSSs") +
  ylab("Total CTCF sites") +
  xlab("Gene TSS distance bin (Kb)") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.title=element_blank(), legend.key = element_rect(size = 6), 
        legend.key.size = unit(2, 'lines'), legend.position = c(0.3,0.85), aspect.ratio = 1  )


################
# Hi-C plot
################

distanceBinnedSizes = as.data.table(table(distanceBinned[,c("bin","significant")]))
# add the % mean increase between significant and non-significant
m = data.table(aggregate(distanceBinned$tssContact, list(significant = distanceBinned$significant, Bin = distanceBinned$bin), mean ) )
distanceBinMeans = round((m[m$significant == 1]$x * 100.0 / m[m$significant == 0]$x) - 100, 0)
v = c()
for (i in seq(1,nrow(distanceBinnedSizes) ) ){
  if (i <= length(distanceBinMeans)){ v = c(v, paste( as.factor(distanceBinMeans[i]), "%", sep = "" ) ) }
  else{ v = c(v,"0")  }
}
distanceBinnedSizes$means = v

g3 = ggplot( data = distanceBinned, aes(x = bin, y = log(tssContact), fill = as.factor(significant)) ) +
  geom_boxplot( outlier.shape = NA, width = 0.9, size = 1) +
  geom_text(data = distanceBinnedSizes[distanceBinnedSizes$significant == 0], aes(bin, -1, label = means), size = 5) +
  ggtitle("Hi-C contacts between gene TSSs") +
  ylab("Hi-C contacts (log-scale)") +
  xlab("Gene TSS distance bin (Kb)") +
  scale_fill_manual( values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        legend.title=element_blank(), legend.key = element_rect(size = 6), 
        legend.key.size = unit(2, 'lines'), legend.position = c(0.7,0.85), aspect.ratio = 1  )


grid.arrange(g1,g2,g3,ncol = 2)


dev.off()

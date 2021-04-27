#26-Nov-2020 Originally from Olivier Delaneau @ UNIL
# Script to compare centimorgan distance between TSSs of COPs and non-COPs

pdf("geuvadis_centimorgan_tss.pdf", 14,7)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

dist0=c()
dist1=c()

D1 = fread("cod_analysis/multiple_features/geuvadis/centimorgan/geuvadis_cops_tss.tsv", head=TRUE, stringsAsFactors=FALSE)
D0 = fread("cod_analysis/multiple_features/geuvadis/centimorgan/geuvadis_non_cops_tss.tsv", head=TRUE, stringsAsFactors=FALSE)

for (CHR in 1:22) {
  M = read.table(paste("cod_analysis/multiple_features/geuvadis/centimorgan/chr", CHR, ".b37.gmap.gz", sep=""), head=TRUE, stringsAsFactors=FALSE)
  
  H0 = D0[chr == CHR]
  H1 = D1[chr == CHR]
  
  H0$cm1 = approx(M$pos, M$cM, xout = H0$tss1, yleft = min(M$cM), yright=max(M$cM))$y
  H0$cm2 = approx(M$pos, M$cM, xout = H0$tss2, yleft = min(M$cM), yright=max(M$cM))$y
  
  H1$cm1 = approx(M$pos, M$cM, xout = H1$tss1, yleft = min(M$cM), yright=max(M$cM))$y
  H1$cm2 = approx(M$pos, M$cM, xout = H1$tss2, yleft = min(M$cM), yright=max(M$cM))$y
  
  dist0 = c(dist0, abs(H0$cm2-H0$cm1))
  dist1 = c(dist1, abs(H1$cm2-H1$cm1))
}

D1$dist = dist1
D0$dist = dist0

mergedData = data.table(centiDist0 = sort(D0$dist), centiDist1 = sort(D1$dist))
mergedDataLog = mergedData
mergedDataLog$centiDist0 = -log10(mergedDataLog$centiDist0)
mergedDataLog$centiDist1 = -log10(mergedDataLog$centiDist1)

g1 = ggplot( data = mergedData, aes(x = centiDist0, y = centiDist1 ) ) +
  geom_point( shape = 1, size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("Non-COP TSS distance (cM)") +
  ylab("COP TSS distance (cM)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=22), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, 
        legend.title=element_blank(), legend.position = c(0.7,0.9), legend.text = element_text(size = 22))

summary(mergedDataLog[!is.infinite(centiDist0)][!is.infinite(centiDist1)]$centiDist0)
summary(mergedDataLog[!is.infinite(centiDist0)][!is.infinite(centiDist1)]$centiDist1)

g2 = ggplot( data = mergedDataLog[!is.infinite(centiDist0)][!is.infinite(centiDist1)], aes(x = centiDist0, y = centiDist1 ) ) +
  geom_point( shape = 1, size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("Non-COP TSS distance (cM)") +
  ylab("COP TSS distance (cM)") +
  scale_x_continuous(labels = c("1","0.1","0.01","1e-3","1e-4","1e-5","1e-6"), breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(labels = c("1","0.1","0.01","1e-3","1e-4","1e-5","1e-6","1e-7"), breaks = c(0,1,2,3,4,5,6,7)) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=22), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, 
        legend.title=element_blank(), legend.position = c(0.7,0.9), legend.text = element_text(size = 22))

layout = rbind( c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2), c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2))
grid.arrange(g1, g2, layout_matrix = layout)

#############
# LD blocks
#############
copFile = "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/CODer_distance_controlled_null.bed_positive"
allFile = "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/post_filter/CODer_raw_results.bed"
copData <- fread(copFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
allData <- fread(allFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

allData$tss1 = apply(allData, 1, function(x) {ifelse(x[6] == "+", x[2], x[3])})
allData$tss2 = apply(allData, 1, function(x) {ifelse(x[11] == "+", x[7], x[8])})

ldFile = "cod_analysis/multiple_features/geuvadis/LD_blocks/fourier_ls-all.bed"
ldData <- fread(ldFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

finalDataset = data.table()
for (c in seq(1,22)) {
  data = allData[`#chr` == c]
  ld = ldData[chr == paste0("chr",c)]
  
  data$tss1bin = cut(as.numeric(data$tss1), breaks=ld$stop)
  data$tss2bin = cut(as.numeric(data$tss2), breaks=ld$stop)
  
  finalDataset = rbind(finalDataset, data)
}

finalDataset$sameBlock = "across"
finalDataset[tss1bin == tss2bin]$sameBlock = "same"

# finalDataset$group = "COP"
# finalDataset[significant == 1]$group = "COP"
# finalDataset[significant == 0]$group = "Non-COP"

finalDataset$group = "All pairs"
finalDataset[pairID %in% copData[significant == 1]$pairID]$group = "COP"
finalDataset[pairID %in% copData[significant == 0]$pairID]$group = "Non-COP"

finalDataset = unique(finalDataset[,.(pairID,group,sameBlock)])

res = data.table(table(finalDataset$sameBlock, finalDataset$group))
colnames(res) = c("LD_block","group","N")

r = data.table(table(finalDataset$group))

percs =  c(
  res[group == "All pairs"][LD_block == "across"]$N / r[V1 == "All pairs"]$N,
  res[group == "All pairs"][LD_block == "same"]$N / r[V1 == "All pairs"]$N,
  
  res[group == "COP"][LD_block == "across"]$N / r[V1 == "COP"]$N,
  res[group == "COP"][LD_block == "same"]$N / r[V1 == "COP"]$N,
  
  res[group == "Non-COP"][LD_block == "across"]$N / r[V1 == "Non-COP"]$N,
  res[group == "Non-COP"][LD_block == "same"]$N / r[V1 == "Non-COP"]$N
)

res$perc = round(percs * 100,1)

res$fill = c(4,5,2,3,0,1)
ggplot(res, aes(x = group, y = perc, fill = as.factor(fill) )) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c("#d9d9d9","#969696","#ccece6","#66C2A5","#ffeda0","#feb24c"), 
                     label = c("Non-COP | Across LD blocks","Non-COP | Same LD block", "COP | Across LD blocks", "COP | Same LD block", "All pairs | Across LD blocks", "All pairs | Same LD block")) +
  annotate(geom = "text", label = paste0(res[LD_block == "across"][group == "COP"]$perc, "%"), x = "COP", y = max(res$perc)*1.08, size = 8) +
  annotate(geom = "text", label = paste0(res[LD_block == "same"][group == "COP"]$perc, "%"), x = "COP", y = max(res$perc)/1.8, size = 8) +
  annotate(geom = "text", label = paste0(res[LD_block == "across"][group == "Non-COP"]$perc, "%"), x = "Non-COP", y = max(res$perc)*1.08, size = 8) +
  annotate(geom = "text", label = paste0(res[LD_block == "same"][group == "Non-COP"]$perc,"%"), x = "Non-COP", y = max(res$perc)/1.8, size = 8) +
  annotate(geom = "text", label = paste0(res[LD_block == "across"][group == "All pairs"]$perc, "%"), x = "All pairs", y = max(res$perc)*0.98, size = 8) +
  annotate(geom = "text", label = paste0(res[LD_block == "same"][group == "All pairs"]$perc, "%"), x = "All pairs", y = max(res$perc)/1.8, size = 8) +
  ylab("% gene pairs") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1)



dev.off()

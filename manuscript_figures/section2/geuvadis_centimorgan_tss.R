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

# lm(mergedData$centiDist0 ~ mergedData$centiDist1)
# lm(mergedDataLog[!is.infinite(centiDist0)][!is.infinite(centiDist1)]$centiDist0 ~ mergedDataLog[!is.infinite(centiDist0)][!is.infinite(centiDist1)]$centiDist1)
#TODO add slope of R/lm in legend

# par(mfrow=c(1,2))
# plot(sort(dist0), sort(dist1), xlab="dist(tss1, tss2) in cM for non-COPs", ylab="dist(tss1, tss2) in cM for COPs", main="Normal scale")
# abline(0,1,col="red")
# plot(sort(-log10(dist0)), sort(-log10(dist1)), xlab="-log10(dist(tss1, tss2)) in cM for non-COPs", ylab="-log10(dist(tss1, tss2)) in cM for COPs", main="Logarithmic scale")
# abline(0,1,col="red")

dev.off()

#12-Aug-2020 Diogo Ribeiro @ UNIL
# Script to analyse results from GWAS eQTL analysis in terms of pleiotropy

pdf("pleiotropy.pdf",9,9)

library(data.table)
library(ggplot2)

inFile = "cod_analysis/GTEx/gwas/pheWAS/phewas_cutoff_0.00000005.tsv"

results = fread( inFile, stringsAsFactors = FALSE, header = F, sep="\t")
colnames(results) = c("variant","freq","trait")
results$freq = as.numeric(results$freq)

paste("Number of GWAS variants passing p-value cutoff:",nrow(results))
paste("% of GWAS variants associated with >1 trait/disease:",round(nrow(results[freq > 1])/nrow(results),2)*100)

#################
# eQTL sharing vs not sharing
#################

eqtlSharedFile = "cod_analysis/GTEx/eqtl_sharing/meta_results_shared.eqtl"
eqtlSharedData = fread( eqtlSharedFile, stringsAsFactors = FALSE, header = F, sep="\t")
eqtlSharedData$V1 = data.table(unlist(lapply(eqtlSharedData$V1, function(x) gsub("_b38","",x) )))
eqtlSharedData = unique(eqtlSharedData[,.(V1)])
eqtlSharedData$shared = 1

eqtlUnsharedFile = "cod_analysis/GTEx/eqtl_sharing/unshared/meta_results_unshared.eqtl"
eqtlUnsharedData = fread( eqtlUnsharedFile, stringsAsFactors = FALSE, header = F, sep=" ")
eqtlUnsharedData$V1 = data.table(unlist(lapply(eqtlUnsharedData$V1, function(x) gsub("_b38","",x) )))
eqtlUnsharedData = unique(eqtlUnsharedData[,.(V1)])
eqtlUnsharedData$shared = 0

ambigouous = eqtlSharedData[V1 %in% eqtlUnsharedData$V1]$V1

boundData = rbind(eqtlSharedData, eqtlUnsharedData)
# boundData = boundData[!V1 %in% ambigouous] ## excluding variants that are shared in one tissue but not shared in other tissues
boundData[V1 %in% ambigouous]$shared = 1 ## if a variant is shared in at least one tissue, it is considered shared
boundData = unique(boundData)

eqtlSharing = merge(results, boundData, by.x = "variant", by.y = "V1")

paste("# shared eQTLs:",length(unique(boundData[shared == 1]$V1)))
paste("# shared eQTLs with a trait",length(unique(eqtlSharing[shared == 1]$variant)))
paste("# shared eQTLs with >1 trait",length(unique(eqtlSharing[shared == 1][freq >1]$variant)))

paste("# unshared eQTLs:",length(unique(boundData[shared == 0]$V1)))
paste("# unshared eQTLs with a trait",length(unique(eqtlSharing[shared == 0]$variant)))
paste("# unshared eQTLs with >1 trait",length(unique(eqtlSharing[shared == 0][freq >1]$variant)))

TP = length(unique(eqtlSharing[shared == 1]$variant))
FP = length(unique(boundData[shared == 1]$V1)) - TP
FN = length(unique(eqtlSharing[shared == 0]$variant)) 
TN = length(unique(boundData[shared == 0]$V1)) - FN
fisher.test(matrix(c( TP, FN, FP, TN),nrow = 2))

TP = length(unique(eqtlSharing[freq > 1][shared == 1]$variant))
FP = length(unique(boundData[shared == 1]$V1)) - TP
FN = length(unique(eqtlSharing[freq > 1][shared == 0]$variant)) 
TN = length(unique(boundData[shared == 0]$V1)) - FN
fisher.test(matrix(c( TP, FN, FP, TN),nrow = 2))


m1 = mean(eqtlSharing[shared == 1]$freq)
m2 = mean(eqtlSharing[shared == 0]$freq)
t1 = wilcox.test(eqtlSharing[shared == 0]$freq,eqtlSharing[shared == 1]$freq) #independent 2-group Mann-Whitney U Test
n1 = nrow(eqtlSharing[shared == 1])
n2 = nrow(eqtlSharing[shared == 0])
ggplot( eqtlSharing, aes(x = freq, fill = as.factor(shared) ) ) +
  geom_density(alpha = 0.8,  outline.type = "full", adjust = 5, size = 1.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("Shared eQTL mean:",round(m1,2)), hjust = 1.07, vjust = 1.5, size = 8.5, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("Other eQTL mean:",round(m2,2)), hjust = 1.07, vjust = 3, size = 8.5, fontface = "bold"  ) +
  # annotate("text", x = Inf, y = Inf, label = paste("Wilcox p-val:",format.pval(t1$p.value,2)), hjust = 1, vjust = 4.5, size = 5, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("# shared:",n1), hjust = 1.1, vjust = 4.5, size = 8.5, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("# other:", n2), hjust = 1.1, vjust = 6, size = 8.5, fontface = "bold"  ) +
  scale_fill_manual(values = c("#ff7f00","#4daf4a"), name = "lead eQTL", labels = c("Other", "Shared")) +
  ggtitle("Lead eQTL trait pleiotropy") +
  xlab("Trait pleiotropy") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), text = element_text(size=28), legend.position = c(0.80, 0.5), 
        legend.text = element_text(size = 24), legend.title = element_text(size = 26),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1 ) 

#################
# Multiple tissues
#################

eqtlFile = "eQTLs/permutation_pass/results/meta_results.eqtl"
eqtlData = fread( eqtlFile, stringsAsFactors = FALSE, header = F, sep=" ")
eqtlData$V8 = data.table(unlist(lapply(eqtlData$V8, function(x) gsub("_b38","",x) )))
paste("Number of eQTLs:",length(unique(eqtlData$V8)))

nTissues = data.table(table(unique(eqtlData[,.(V8,V21)])$V8))
resultsEQTL = results[variant %in% nTissues$V1]

nTissueCutoff = 1

# Multiple tissues
resultsEQTL$multipleTissues = 0
resultsEQTL[variant %in% nTissues[N>nTissueCutoff]$V1]$multipleTissues = 1

# Fisher test
TP = length(unique(resultsEQTL[multipleTissues == 1]$variant))
FP = length(unique(nTissues[N>nTissueCutoff]$V1)) - TP
FN = length(unique(resultsEQTL[multipleTissues == 0]$variant))
TN = length(unique(nTissues[N<=nTissueCutoff]$V1)) - FN
fisher.test(matrix(c( TP, FN, FP, TN),nrow = 2))


# Distribution
m1 = mean(resultsEQTL[multipleTissues == 0]$freq)
m2 = mean(resultsEQTL[multipleTissues == 1]$freq)
wilcox.test(resultsEQTL[multipleTissues == 0]$freq,resultsEQTL[multipleTissues == 1]$freq)
n1 = nrow(resultsEQTL[multipleTissues == 0])
n2 = nrow(resultsEQTL[multipleTissues == 1])
ggplot( resultsEQTL, aes(x = as.numeric(freq), fill = as.factor(multipleTissues) ) ) +
  geom_density(alpha = 0.8,  outline.type = "full", adjust = 5, size = 1.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste(">1 tissue mean:",round(m2,2)), hjust = 1, vjust = 1.5, size = 8, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("one tissue mean:",round(m1,2)), hjust = 1, vjust = 3, size = 8, fontface = "bold"  ) +
  # annotate("text", x = Inf, y = Inf, label = paste("Wilcox p-val:",format.pval(t1$p.value,2)), hjust = 1, vjust = 4.5, size = 5, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("# more:", n2), hjust = 1, vjust = 4.5, size = 8, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = paste("# one:",n1), hjust = 1, vjust = 6, size = 8, fontface = "bold"  ) +
  scale_fill_manual(values = c("#ff7f00","#4daf4a"), name = "Tissues", labels = c("One", "More")) +
  # ggtitle("eQTLs: one tissue vs more tissues ") +
  xlab("Trait pleiotropy") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), legend.position = c(0.75, 0.5), aspect.ratio = 1 )


############
# Tissue & shared together
############
# Many shared eQTLs are present in multiple tissues

mergedData = merge(eqtlSharing,resultsEQTL, by = "variant")

data.table(table(mergedData$shared,mergedData$multipleTissues))

2231/(2231+8144)
2070/(2070+2489)

m1 = mean(mergedData[multipleTissues == 0][shared == 0]$freq.x)
m2 = mean(mergedData[multipleTissues == 1][shared == 0]$freq.x)
m3 = mean(mergedData[multipleTissues == 0][shared == 1]$freq.x)
m4 = mean(mergedData[multipleTissues == 1][shared == 1]$freq.x)

wilcox.test(mergedData[multipleTissues == 0][shared == 0]$freq.x,mergedData[multipleTissues == 0][shared == 1]$freq.x)
wilcox.test(mergedData[multipleTissues == 1][shared == 0]$freq.x,mergedData[multipleTissues == 1][shared == 1]$freq.x)

summary(aov(freq.x ~ shared * multipleTissues, data = mergedData))
summary(aov(freq.x ~ multipleTissues * shared, data = mergedData))
summary(lm(mergedData$freq.x ~ mergedData$shared + mergedData$multipleTissues))

dev.off()

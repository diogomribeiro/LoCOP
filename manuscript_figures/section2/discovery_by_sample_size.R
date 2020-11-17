
pdf("discovery_by_sample_size.pdf",13,12)

library(data.table)
library(ggplot2)

statsData = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/COP_specificity_stats.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")
copData = fread("zcat /scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/GTEx/CODer_final_dataset_cops_merged.bed.gz", stringsAsFactors = FALSE, header = TRUE, sep="\t")

coderStats = fread("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/GTEx/coder_stats.txt", stringsAsFactors = FALSE, header = F, sep="\t")
coderStats = coderStats[V2 == "Total samples"][,.(V1,V3)]
colnames(coderStats) = c("tissue","sample_size")

colorCode = fread( "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/GTEx/v8/other/color_code_computer.txt", stringsAsFactors = FALSE, header = FALSE, sep="\t")


mergedData = unique(merge(statsData,copData,by="pairID"))
mergedData = unique(merge(mergedData,coderStats,by="tissue"))

data = data.table(table(mergedData[copFreq == 1]$tissue)) #ratio > 0.2][ratio < 0.5
others = data.table(table(mergedData$tissue))
data = merge(data, others, by = "V1")
colnames(data) = c("tissue","uniques","total")
data$perc = data$uniques * 100.0 / data$total
data = merge(data,coderStats,by = "tissue")

correlation = cor.test(data$perc, data$sample_size, method = "spearman")
correlationText = paste("Spearman R:",round(correlation$estimate,2), "P-value:",format.pval(pv = c(correlation$p.value), digits = 2), sep = " ")
d = lm(data$perc ~ data$sample_size)
ggplot( data, aes(x = sample_size, y = perc, color = tissue, fill = tissue ) ) +
  geom_abline(slope = d$coefficients[2], intercept = d$coefficients[1], color = "#de2d26") +
  geom_point( show.legend = FALSE, alpha = 0.9, size = 4, shape = 21, color = "black") +
  geom_text( aes(label = tissue), size = 7, vjust = 1.5, check_overlap = TRUE, show.legend = FALSE) +
  annotate("text", x = 420, y = 38, label = correlationText, hjust = 0, vjust =1, size = 8.5  ) +
  # ggtitle("% unique COPs per tissue sample size") +
  ylab("% unique COPs") +
  xlab("Tissue sample size") +
  scale_color_manual(values = colorCode$V2) +
  scale_fill_manual(values = colorCode$V2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), axis.text.x = element_text(angle = 0, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2))

dev.off()

pdf(paste("pleiotropy_enrichment.pdf",sep=""),7,7)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

inFile = "cod_analysis/GTEx/candidates/enrichment_ensembl_v101/comparison/meta_results.out"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")

mergedData$odds = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[3]),as.numeric(x[5])-as.numeric(x[3]),as.numeric(x[4]),as.numeric(x[6])-as.numeric(x[4]) ),nrow = 2), conf.level = 0.95)$estimate )
mergedData$confmin = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[3]),as.numeric(x[5])-as.numeric(x[3]),as.numeric(x[4]),as.numeric(x[6])-as.numeric(x[4]) ),nrow = 2), conf.level = 0.95)$conf.int[1] )
mergedData$confmax = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[3]),as.numeric(x[5])-as.numeric(x[3]),as.numeric(x[4]),as.numeric(x[6])-as.numeric(x[4]) ),nrow = 2), conf.level = 0.95)$conf.int[2] )
mergedData$pval = apply(mergedData, 1, function(x) fisher.test(matrix(c(as.numeric(x[3]),as.numeric(x[5])-as.numeric(x[3]),as.numeric(x[4]),as.numeric(x[6])-as.numeric(x[4]) ),nrow = 2), conf.level = 0.95)$p.value )

ggplot(mergedData, aes(x = odds, y = tag, fill = tag)) +
  geom_segment(data = mergedData, aes(x = confmin, xend = confmax, y = tag, yend = tag, color = tag), size = 1 ) +
  geom_point(size = 4, shape = 21) +
  geom_text(data = mergedData, aes(x = odds, y = tag, label = paste("OR=",round(odds,1), sep = "" ), color = tag ), position = position_nudge(y = 0.25), size = 4, fontface = "bold") +
  geom_text(data = mergedData, aes(x = odds, y = tag, label = paste("P=",format.pval(pval, digits = 1), sep = "" ), color = tag ), position = position_nudge(y = -0.25), size = 4, fontface = "bold") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # ggtitle("Shared lead eQTLs vs not shared") +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf","black")) +
  scale_color_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf","black")) +
  # xlim(c(0,min(max(mergedData$confmax), 4.5)) ) +
  xlab("Odds ratio") +
  ylab("Annotation") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "None", panel.background = element_rect(colour = "black", fill = "white", size = 1),
        #axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank()
        )

## TODO: order annotations

dev.off()

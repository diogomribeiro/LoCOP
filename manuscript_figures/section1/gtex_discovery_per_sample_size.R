#27-05-2020 # Diogo Ribeiro @ UNIL
# Script to plot CODer.py summary results from several runs

pdf("gtex_discovery_per_sample_size.pdf",13,12)

library(data.table)
require(ggplot2)
library(tidyr)

statsFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/GTEx/coder_stats.txt"
copsFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/GTEx/CODer_cops_merged.bed"
colorFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/GTEx/v8/other/color_code_computer.txt"

statsData = fread( statsFile, stringsAsFactors = FALSE, header = FALSE, sep="\t")
copsData = fread( copsFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
colorCode = fread( colorFile, stringsAsFactors = FALSE, header = FALSE, sep="\t")

colnames(statsData) = c("Tissue","Metric","Value")

percs = c(statsData[statsData$Metric == "Total Genes in CODs"]$Value * 100.0 / statsData[statsData$Metric == "Total Genes tested"]$Value)

statsData$Value = as.numeric(statsData$Value)
statsData[statsData$Metric == "Percentage genes in CODs"]$Value = percs
statsData$tissueName = data.table(unlist(lapply(statsData$Tissue, function(x) unlist(gsub("_"," ",x[1]))[1])))$V1

##################
# Scatterplot # samples vs # COPs
##################
df1 = statsData[statsData$Metric == "Total samples"]
df2 = statsData[statsData$Metric == "Percentage genes in CODs"]
currentDF = data.table( Tissue = df1$tissueName, samples = df1$Value, perc_cods = df2$Value)
correlation = cor.test(currentDF$samples, currentDF$perc_cods, method = "spearman")
correlationText = paste("Spearman correlation:",round(correlation$estimate,2), "\nP-value:",format.pval(pv = c(correlation$p.value), digits = 2), sep = " ")

d = lm(currentDF$perc_cods ~ currentDF$samples)

ggplot( currentDF, aes(x = samples, y = perc_cods, color = Tissue, fill = Tissue ) ) +
  geom_abline(slope = d$coefficients[2], intercept = d$coefficients[1], color = "#de2d26") +
  geom_point( show.legend = FALSE, alpha = 0.9, size = 4, shape = 21, color = "black") +
  geom_text( aes(label = Tissue), size = 7, vjust = 1.5, check_overlap = TRUE, show.legend = FALSE) +
  annotate("text", x = 480, y = 20, label = correlationText, hjust = 0, vjust =1, size = 8.5  ) +
  # ggtitle("% co-expressed genes vs tissue sample size ") +
  ylab("% co-expressed genes") +
  xlab("Tissue sample size") +
  scale_fill_manual(values = colorCode$V2) +
  scale_color_manual(values = colorCode$V2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), axis.text.x = element_text(angle = 0, vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2))

## Overlap between GTEx and Geuvadis LCL gene pairs
# fisher.test(matrix(c(3069,6378,7202,181735),nrow = 2))
fisher.test(matrix(c(3069, 4133, 3309, 181735-(3069+4133+3309)), nrow = 2))

dev.off( )



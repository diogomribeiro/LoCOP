#27-05-2020 # Diogo Ribeiro @ UNIL
# Script to plot CODer.py summary results from several runs

pdf("gtex_cop_discovery.pdf",13,12)

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
# COPs per tissue
##################
ggplot( statsData[Metric == "Total Genes in CODs"], aes(x = reorder(tissueName, Value), y = Value, fill = tissueName, label = Value ) ) + #Percentage genes in CODs
  geom_bar( stat = "identity", color = "black", show.legend = FALSE, alpha = 0.5) +
  geom_text(hjust = 1, size = 6) +
  scale_fill_manual(values = colorCode$V2) +
  coord_flip() +
  ylab("Number of co-expressed genes") +
  xlab("Tissue") +
  # ggtitle("Number of co-expressed genes per tissue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1))


dev.off()


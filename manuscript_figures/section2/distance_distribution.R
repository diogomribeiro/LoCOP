#28-05-2020 # Diogo Ribeiro @ UNIL
# Script to plot distance distribution

pdf("distance_distribution.pdf",9,7)

library(data.table)
# require(ggplot2)
library(ggplot2, lib = "/home/dribeiro/Software/R/lib") ## newer ggplot2 version with geom_density outline.type = "full" and size bugfix

resData = fread( "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_raw_results.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")
null = fread( "cod_identification/geuvadis/all_chr/final_fdr0.01/final_dataset/CODer_distance_controlled_null.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")

# resData = fread( "cod_identification/GTEx/all_tissues/Lung/CODer_raw_results.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")
# null = fread( "cod_identification/GTEx/all_tissues/Lung/final_dataset/CODer_distance_controlled_null.bed", stringsAsFactors = FALSE, header = TRUE, sep="\t")

resData$cop = 0
resData[pairID %in% null[significant == 1]$pairID]$cop = 1

summary(resData$distance)

ggplot( data = resData, aes(x = distance, fill = as.factor(cop) ) ) +
  geom_density( alpha = 0.8, size = 0.8, outline.type = "full") +
  # ggtitle("Distance between gene TSS") +
  xlab("Absolute distance (bp)") +
  scale_fill_manual(limits = rev(levels(resData$cop)), values = c("#d9d9d9","#66C2A5"), label = c("Not co-expressed", "Co-expressed pairs") ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, 
        legend.title=element_blank(), legend.position = c(0.7,0.9), legend.text = element_text(size = 22))

################
# Gene pair calculations
################
dist = 200000
allCOPs = length(unique((resData[cop == 1]$pairID)))
allCOPsWithin = length(unique((resData[cop == 1][distance <= dist]$pairID)))
allNonCOPs = length(unique((resData[cop == 0]$pairID)))
allNonCOPsWithin = length(unique((resData[cop == 0][distance <= dist]$pairID)))

# proportion of co-expressed gene pairs within 200kb
allCOPsWithin/allCOPs
# proportion of non-co-expressed gene pairs within 200kb
allNonCOPsWithin/allNonCOPs
# How many gene pairs within are co-expressed
allCOPsWithin / (allCOPsWithin + allNonCOPsWithin)

summary(resData[cop == 1]$distance)
summary(resData[cop == 0]$distance)


dev.off( )

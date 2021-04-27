#07-Aug-2020 Diogo Ribeiro @ UNIL
# Script to read GWAS results and plot qq plots

pdf("qqplot_single_tissues.pdf",16,9)

library(data.table)
library(ggplot2)
# install.packages("remotes")
# remotes::install_github("sinarueeger/ggGWAS")
library(ggGWAS)
library(grid)
library(gridExtra)
# source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

inFile = "zcat cod_analysis/GTEx/gwas/all_gwas/35_traits/results/meta.non_redundant_reformatted.tsv.gz"

results = fread( inFile, stringsAsFactors = FALSE, header = F, sep="\t")
colnames(results) = c("variant","value","tissue","category","trait")

############
# Plots for the combinations of gwas and tissue
############

results = results[category == "shared_eqtl" | category == "unshared_eqtl" | category == "gwas_sample"]
results[category == "shared_eqtl"]$category = "1"
results[category == "unshared_eqtl"]$category = "2"
results[category == "gwas_sample"]$category = "3"


lambda1 = median(qchisq(1 - results[trait == "hypothyroidism"][tissue == "Thyroid"][category == "1"]$value,1))/qchisq(0.5,1)
lambda2 = median(qchisq(1 - results[trait == "hypothyroidism"][tissue == "Thyroid"][category == "2"]$value,1))/qchisq(0.5,1)
lambda3 = median(  chisq <- qchisq(1 - results[trait == "hypothyroidism"][category == "3"]$value,1))/qchisq(0.5,1)

g1 = ggplot( results[trait == "hypothyroidism"][tissue == "Thyroid"], aes(y = value, color = category ) ) +
  geom_gwas_qq(alpha = 0.8, size = 3.5, shape = 16) +
  geom_gwas_qq(data = results[trait == "hypothyroidism"][tissue == "gwas_sample"], aes(color = tissue), alpha = 0.8, size = 3.5, shape = 16) +
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0, y = Inf, label = paste("Shared eQTL inflation:", round(lambda1,2)), hjust = 0, vjust = 2.5, size = 7.5, fontface = "bold"  ) +
  annotate("text", x = 0, y = Inf, label = paste("Other eQTL inflation:", round(lambda2,2)), hjust = 0, vjust = 4, size = 7.5, fontface = "bold"  ) +
  annotate("text", x = 0, y = Inf, label = paste("GWAS (sample) inflation:", round(lambda3,2)), hjust = 0, vjust = 5.5, size = 7.5, fontface = "bold"  ) +
  scale_color_brewer(palette = "Set2", name = "Variants", labels = c("Shared eQTL", "Other eQTL","GWAS (sample)")) +
  ggtitle("Hypothyroidism in Thyroid") +
  xlab("Expected -log10(p-value)") +
  ylab("Observed -log10(p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        legend.position = c(0.22,0.72), legend.title = element_blank(), legend.text = element_text(size = 22), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1)

lambda1 = median(qchisq(1 - results[trait == "systolic_blood_pressure"][tissue == "Artery_Aorta"][category == "1"]$value,1))/qchisq(0.5,1)
lambda2 = median(qchisq(1 - results[trait == "systolic_blood_pressure"][tissue == "Artery_Aorta"][category == "2"]$value,1))/qchisq(0.5,1)
lambda3 = median(  chisq <- qchisq(1 - results[trait == "systolic_blood_pressure"][category == "3"]$value,1))/qchisq(0.5,1)

g2 = ggplot( results[trait == "systolic_blood_pressure"][tissue == "Artery_Aorta"], aes(y = value, color = category ) ) +
  geom_gwas_qq(alpha = 0.8, size = 3.5, shape = 16) +
  geom_gwas_qq(data = results[trait == "systolic_blood_pressure"][tissue == "gwas_sample"], aes(color = tissue), alpha = 0.8, size = 3.5, shape = 16) +
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0, y = Inf, label = paste("Shared eQTL inflation:", round(lambda1,2)), hjust = 0, vjust = 2.5, size = 7.5, fontface = "bold"  ) +
  annotate("text", x = 0, y = Inf, label = paste("Other eQTL inflation:", round(lambda2,2)), hjust = 0, vjust = 4, size = 7.5, fontface = "bold"  ) +
  annotate("text", x = 0, y = Inf, label = paste("GWAS (sample) inflation:", round(lambda3,2)), hjust = 0, vjust = 5.5, size = 7.5, fontface = "bold"  ) +
  scale_color_brewer(palette = "Set2", name = "Variants", labels = c("Shared eQTL", "Other eQTL","GWAS (sample)")) +
  ggtitle("Systolic Blood Pressure in Artery Aorta") +
  xlab("Expected -log10(p-value)") +
  ylab("Observed -log10(p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        legend.position = c(0.22,0.72), legend.title = element_blank(), legend.text = element_text(size = 22), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1)


grid.arrange(g1,g2,nrow = 1)

dev.off()

# subset[category == "shared_eqtl"]$category = "1: shared_eqtl"
# subset[category == "unshared_eqtl"]$category = "2: not_shared_eqtl"
# subset[category == "gwas_sample"]$category = "3: gwas_sample"
# 
# table(subset$category)
# table(subset$trait)
# table(subset$tissue)
# 
# wantedPairs = c(
#   "hypothyroidism|Thyroid",
#   # "white_blood_cell|Whole_Blood",
#   "systolic_blood_pressure|Artery_Aorta",
#   # "myocardial_infarction|Artery_Aorta"
# )
# 
# listPlot = list()
# i=0
# for (pair in wantedPairs){
#   print(pair)
#   wantedTrait = unlist(strsplit(pair,"[|]"))[1]
#   wantedTissue = unlist(strsplit(pair,"[|]"))[2]
# 
#   lambda1 = median(qchisq(1 - subset[tissue == wantedTissue][trait == wantedTrait][category == "1: shared_eqtl"]$value,1))/qchisq(0.5,1)
#   lambda2 = median(qchisq(1 - subset[tissue == wantedTissue][trait == wantedTrait][category == "2: not_shared_eqtl"]$value,1))/qchisq(0.5,1)
#   lambda3 = median(  chisq <- qchisq(1 - subset[trait == wantedTrait][category == "3: gwas_sample"]$value,1))/qchisq(0.5,1)
#   
#   i=i+1
#   listPlot[[i]] = ggplot( subset[tissue == wantedTissue][trait == wantedTrait], aes(y = value, color = category ) ) +
#       geom_gwas_qq(alpha = 0.8, size = 2, shape = 16) +
#       geom_gwas_qq(data = subset[tissue == "gwas_sample"][trait == wantedTrait], aes(color = tissue), alpha = 0.8, size = 2, shape = 16) +
#       geom_abline(slope = 1, intercept = 0) +
#       annotate("text", x = 0, y = Inf, label = paste("Shared inflation:", round(lambda1,2)), hjust = 0, vjust = 1.5, size = 4.5, fontface = "bold"  ) +
#       annotate("text", x = 0, y = Inf, label = paste("Not shared inflation:", round(lambda2,2)), hjust = 0, vjust = 3, size = 4.5, fontface = "bold"  ) +
#       annotate("text", x = 0, y = Inf, label = paste("GWAS inflation:", round(lambda3,2)), hjust = 0, vjust = 4.5, size = 4.5, fontface = "bold"  ) +
#       scale_color_brewer(palette = "Set2", name = "Variants", labels = c("Shared eQTL", "Not shared eQTL","GWAS (sample)")) +
#       ggtitle(pair) + 
#       xlab("Expected -log10(p-value)") +
#       ylab("Observed -log10(p-value)") +
#       theme_minimal() +
#       theme(plot.title = element_text(hjust = 0.5, size = 22), text = element_text(size=16),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1)
#   
# }
# grid_arrange_shared_legend(listPlot)




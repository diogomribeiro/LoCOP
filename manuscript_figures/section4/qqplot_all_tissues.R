#14-Aug-2020 Diogo Ribeiro @ UNIL
# Script to read GWAS results and plot qq plots

args = commandArgs(trailingOnly=TRUE)

# run = "part1"
run = args[1]

pdf(paste("qqplot_all_tissues_",run,".pdf",sep = ""),16,7)

library(data.table)
library(ggplot2)
# install.packages("remotes")
# remotes::install_github("sinarueeger/ggGWAS")
library(ggGWAS)
library(grid)
library(gridExtra)
source("/users/dribeir1/code/cod/src/cod/util/r_utils/regression_functions.R")

inFile = "zcat /scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/gwas/all_gwas/35_traits_all_tissues/results/meta.non_redundant.tsv.gz"
results = fread( inFile, stringsAsFactors = FALSE, header = F, sep="\t")
colnames(results) = c("variant","value","category","trait")

subset = results[category == "all_tissues:shared_eqtl" | category == "all_tissues:unshared_eqtl" | category == "gwas_sample"]

# Removing unshared cases which are shared
d1 = subset[category == "all_tissues:shared_eqtl"]
d2 = subset[category == "all_tissues:unshared_eqtl"][!variant %in% d1$variant]
subset = rbind(d1,d2,subset[category == "gwas_sample"])

############
# Plots for the combinations of gwas and tissue
############

# quantiles = data.table(aggregate(subset[,pvalue_log10],list(subset$category, subset$trait), FUN = quantile, probs=c(0.5, 0.75, 0.9, 0.95)))
traitFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/GWAS/uk_biobank/list_traits.txt"
traits = fread( traitFile, stringsAsFactors = FALSE, header = F, sep="\t")
traits = traits[order(V1)]

if (run == "part1"){
  traits = traits[1:10]
}
if (run == "part2"){
  traits = traits[11:20]
}
if (run == "part3"){
  traits = traits[21:30]
}
if (run == "part4"){
  traits = traits[31:35]
}

listPlot = list()
for (i in seq(nrow(traits) ) ){
  print(paste(i,traits[i]))

  wantedTrait = traits[i]
  
  displayName = wantedTrait
  titleSize = 20

  if (wantedTrait == "age_hay_fever_rhinitis_eczema_diagnosed"){
    displayName = paste("age_hay_fever","rhinitis_eczema_diagnosed",sep = "\n")
    titleSize = titleSize - 3
  }
  if (wantedTrait == "age_started_wearing_glasses_or_contact_lenses"){
    displayName = paste("age_started_wearing_glasses","_or_contact_lenses",sep="\n")
    titleSize = titleSize - 3
  }
  if (wantedTrait == "allergy_adverse_effect_of_penicillin"){
    displayName = paste("allergy_adverse","effect_of_penicillin",sep="\n")
    titleSize = titleSize - 3
  }
  if (wantedTrait == "esophagitis_gerd_and_related_diseases"){
    displayName = paste("esophagitis_gerd","_and_related_diseases",sep="\n")
    titleSize = titleSize - 3
  }

  displayName = gsub("_"," ",displayName)

  lambda1 = median(qchisq(1 - subset[trait == wantedTrait][category == "all_tissues:shared_eqtl"]$value,1))/qchisq(0.5,1)
  lambda2 = median(qchisq(1 - subset[trait == wantedTrait][category == "all_tissues:unshared_eqtl"]$value,1))/qchisq(0.5,1)
  lambda3 = median(  chisq <- qchisq(1 - subset[trait == wantedTrait][category == "gwas_sample"]$value,1))/qchisq(0.5,1)

  listPlot[[i]] = ggplot( subset[trait == wantedTrait], aes(y = value, color = category ) ) +
    geom_gwas_qq(alpha = 0.8, size = 2, shape = 16) +
    geom_abline(slope = 1, intercept = 0) +
    annotate("text", x = 0, y = Inf, label = paste("Shared:", round(lambda1,2)), hjust = 0, vjust = 1.5, size = 5, fontface = "bold", color = "#4daf4a"  ) +
    annotate("text", x = 0, y = Inf, label = paste("Other:", round(lambda2,2)), hjust = 0, vjust = 3, size = 5, fontface = "bold", color = "#ff7f00"  ) +
    annotate("text", x = 0, y = Inf, label = paste("GWAS:", round(lambda3,2)), hjust = 0, vjust = 4.5, size = 5, fontface = "bold", color = "#377eb8"  ) +
    scale_color_brewer(palette = "Set2", name = "Variants", labels = c("Shared eQTL", "Not shared eQTL","GWAS (sample)")) +
    ggtitle(displayName) +
    # xlab("Expected") +
    # ylab("Observed") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = titleSize), text = element_text(size=16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
          axis.title = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1)

}

if (length(traits$V1) > 5){
  grid.arrange(grobs = listPlot, nrow = 2)
} else {
  grid.arrange(grobs = listPlot, nrow = 1)
}



# grid_arrange_shared_legend(listPlot)



dev.off()


# ggplot( subset, aes(y = value, color = category ) ) +
#   geom_gwas_qq(alpha = 0.8, size = 1.5, shape = 16) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_color_brewer(palette = "Set2", name = "Variants", labels = c("Shared eQTL", "Not shared eQTL","GWAS (sample)")) +
#   xlab("Expected -log10(p-value)") +
#   ylab("Observed -log10(p-value)") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5, size = 14), text = element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
#         panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1) +
#   facet_wrap(~trait, scales = "free")

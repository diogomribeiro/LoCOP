# 25-May-2020 Diogo Ribeiro
# Plotting results from fisher exact test and its odds ratio on groups

pdf("functional_related_odds_ratio.pdf",12,10)

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)

######################
# Read input files
######################

# args = commandArgs(trailingOnly=TRUE)
# if (length(args)<3) {
#   stop("All arguments need to be provided", call.=FALSE)
# }

inputFile = "cod_analysis/geuvadis/same_function/results.txt"
# inputFile = "cod_analysis/geuvadis/same_function/after_paralog_filter/results.txt"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetWanted = dataset[ItemGroup != "High corr"]
datasetWanted = datasetWanted[ExternalList != "functional_related"]

colnames(datasetWanted) = c("ExternalList","ExtSize","ItemGroup","GrpSize","Overlap","OddsRatio","Pvalue")

# number of digits after comma
options("scipen"=100, "digits"=2)

# Python float point limit is 2.225e-308
datasetWanted[Pvalue == 0]$Pvalue = 2.225e-308
datasetWanted$log10pval = -log10(datasetWanted$Pvalue)

datasetWanted[ExternalList == "go_sharing"]$ExternalList = "GO sharing"
datasetWanted[ExternalList == "evolutionarily_related"]$ExternalList = "Paralogs"
datasetWanted[ExternalList == "same_pathway"]$ExternalList = "Same pathway"
datasetWanted[ExternalList == "same_complex"]$ExternalList = "Same complex"
datasetWanted[ItemGroup == "Neg corr"]$ItemGroup = "Neg. corr."
datasetWanted[ItemGroup == "Pos corr"]$ItemGroup = "Pos. corr."

# Odds ratio plot
ggplot( datasetWanted, aes(x = ExternalList, y = ItemGroup, fill = OddsRatio)) +
  geom_tile( colour = "black", size = 1) + 
  # scale_fill_continuous( low = "white", high = "#de2d26", name = "Fisher's Exact Test \nOdds ratio", na.value = "#de2d26") +
  scale_fill_gradient2(low = "#31a354", mid = "white", high = "#4daf4a", midpoint = 1,  name = "Fisher's Exact Test \nOdds ratio", na.value = "white") +
  xlab("Functional group") +
  ylab("COP category") +
  coord_flip() +
  geom_text( label = paste("N=",datasetWanted$GrpSize, "\nNO=",datasetWanted$Overlap, "\nPv=",round(datasetWanted$log10pval,1), "\nOR=",datasetWanted$OddsRatio, sep=""), size = 6) +
  theme_minimal() + 
  theme(text = element_text(size=20),axis.title.x=element_text(vjust=-0.6), axis.text.x=element_text(angle=45, vjust = 0.5), aspect.ratio = 1)

dev.off( )

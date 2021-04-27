#24-05-2019 # Diogo Ribeiro @ UNIL
# Script to plot GTEx COP matrix

pdf("tissue_clustering.pdf")
# png("tissue_clustering.png")

library(data.table)
library(gplots)
require(ggplot2)
library(tidyr)

# matrixFile = "cod_analysis/GTEx/tissue_clustering/v8/tissue_pair_perc_share.tsv"
matrixFile = "cod_analysis/GTEx/tissue_clustering/v8/tissue_pair_perc_share_before_filters.tsv"
colorFile = "raw_input/GTEx/v8/other/color_code_computer.txt"

matrixData = fread( matrixFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
colorCode = fread( colorFile, stringsAsFactors = FALSE, header = FALSE, sep="\t")

matrixData$Tissue1 = data.table(unlist(lapply(matrixData$Tissue1, function(x) unlist(gsub("_"," ",x[1]))[1])))$V1
matrixData$Tissue2 = data.table(unlist(lapply(matrixData$Tissue2, function(x) unlist(gsub("_"," ",x[1]))[1])))$V1

matrix = spread(matrixData, Tissue2, pi1)
matrix[is.na(matrix)] = 0
matrix = as.matrix(matrix, rownames = 1)

dist.pear <- function(x) as.dist(1-cor(t(x)))
# hclust.ave <- function(x) hclust(x, method="average")
my_palette <- colorRampPalette(c("#ece7f2", "#023858"))(n = 50)

heatmap.2(matrix,
          ColSideColors=colorCode$V2,
          RowSideColors=colorCode$V2,
          # density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          cexRow=0.42,
          cexCol=0.42,
          # sepwidth=c(0.001,0.001),
          # sepcolor="black",
          # colsep=NA,
          # rowsep=NA,
          # dendrogram="both",
          dendrogram="none",
          # main = "", #"Enrichment of lincRNA interactions on network modules",
          margins=c(23,23),
          # Colv="NA",    # turn off column clustering
          # Rowv="NA",
          offsetRow = 0,
          offsetCol = 0,
          distfun = dist.pear,
          hclustfun = hclust,
          key = F, # remove color key
          key.title = "COPs shared",
          key.xlab = "%",
          col=my_palette,       # use on color palette defined earlier
          # lmat = rbind( c(1,1,1), c(1,1,1) ),
          lhei = c(1.5, 4),
          lwid = c(1.5, 4, 0.75)
)


dev.off()

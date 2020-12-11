#20-Nov-2020 Diogo Ribeiro @ UNIL
# Script to analyse results from coloc across all COPs

pdf("coloc_results.pdf",16,9)

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

genePairInfo = fread( "cod_analysis/geuvadis/eQTLs/coloc/prepare_input_files/geuvadis_cops_noncop_info.tsv", stringsAsFactors = FALSE, header = T, sep="\t")

prob = 0.5

inFile = "cod_analysis/geuvadis/eQTLs/coloc/full_run/meta_results.out"
mergedData = fread( inFile, stringsAsFactors = FALSE, header = F, sep=" ")
mergedData$V1 = NULL
colnames(mergedData) = c("H0","H1","H2","H3","H4","pairID")

mergedData = merge(mergedData, genePairInfo, by = "pairID", all.x = T)

#############
# H4 and eQTL sharing
mergedDataCOP = mergedData[significant == 1]
mergedDataNonCOP = mergedData[significant == 0]

mergedDataCOP = melt(mergedDataCOP[,.(pairID,H0,H1,H2,H3,H4,eqtlSharing)], id.vars = c("pairID","eqtlSharing"))
mergedDataCOP = mergedDataCOP[variable == "H4"]
mergedDataCOP$call = "No"
mergedDataCOP[value > prob]$call = "Yes"

mergedDataCOP$sharing = ""
mergedDataCOP[eqtlSharing == 2]$sharing = "eQTL sharing"
mergedDataCOP[eqtlSharing == 1]$sharing = "eQTL sharing"
mergedDataCOP[eqtlSharing == 0]$sharing = "No sharing"

dt = data.table(table(mergedDataCOP$call, mergedDataCOP$sharing))
colnames(dt) = c("H4","sharing","N")
perc1 = dt[sharing == "eQTL sharing"][H4 == "Yes"]$N * 100.0 / (dt[sharing == "eQTL sharing"][H4 == "Yes"]$N + dt[sharing == "eQTL sharing"][H4 == "No"]$N)
perc2 = dt[sharing == "No sharing"][H4 == "Yes"]$N * 100.0 / (dt[sharing == "No sharing"][H4 == "Yes"]$N + dt[sharing == "No sharing"][H4 == "No"]$N)
g1 = ggplot( data = dt, aes(x = as.factor(sharing), y = N, fill = as.factor(H4) ) ) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c("#ccece6","#66C2A5")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "eQTL sharing", y = 600, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "eQTL sharing", y = 2000, size = 8) +
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "No sharing", y = 50, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "No sharing", y = 1900, size = 8) +
  ylab("# gene pairs") +
  xlab("COPs and eQTL sharing") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white",
                                        size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 5), legend.title = element_blank(),
        legend.key.size = unit(1.7, 'lines'), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")  )


#####
# COPs vs non-COPs
mergedData = mergedData[eGenes == 2]

mergedData$group = "COP"
mergedData[significant == 0]$group = "Non-COP"
meltedData = melt(mergedData[,.(pairID,H0,H1,H2,H3,H4,group,eGenes)], id.vars = c("pairID","group", "eGenes"))
meltedData$call = 0
meltedData[value > prob]$call = 1
meltedData = meltedData[call == 1]

dt = data.table(table(meltedData$variable,meltedData$group))
colnames(dt) = c("Hypothesis","group","N")

dt = dt[Hypothesis == "H4"]
dt = rbind(dt, data.table(Hypothesis = "Other", group = "COP", N = nrow(mergedData[significant == 1]) - sum(dt[group == "COP"]$N)))
dt = rbind(dt, data.table(Hypothesis = "Other", group = "Non-COP", N = nrow(mergedData[significant == 0]) - sum(dt[group == "Non-COP"]$N)))
dt$fill = c(3,1,2,0)
perc1 = dt[group == "COP"][Hypothesis == "H4"]$N * 100.0 / (dt[group == "COP"][Hypothesis == "H4"]$N + dt[group == "COP"][Hypothesis == "Other"]$N)
perc2 = dt[group == "Non-COP"][Hypothesis == "H4"]$N * 100.0 / (dt[group == "Non-COP"][Hypothesis == "H4"]$N + dt[group == "Non-COP"][Hypothesis == "Other"]$N)
g2 = ggplot( data = dt, aes(x = group, y = N, fill = as.factor(fill) ) ) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c("#d9d9d9","#969696","#ccece6","#66C2A5"), label = c("Non-COP | Not colocalized","Non-COP | Colocalized", "COP | Not colocalized", "COP | Colocalized")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "COP", y = 500, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "COP", y = 1700, size = 8) +
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "Non-COP", y = 50, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "Non-COP", y = 650, size = 8) +
  ylab("# gene pairs") +
  xlab("Gene pair group") +
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 28), text = element_text(size=26),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white",
                                        size = 1.5, linetype = "solid"),
        legend.key = element_rect(size = 5), legend.title = element_blank(),
        legend.key.size = unit(1.7, 'lines'), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")  )

lay = rbind(c(1,1,1,1,2,2,2,2,2,2,2), c(1,1,1,1,2,2,2,2,2,2,2))
grid.arrange(g1,g2, layout_matrix = lay)

dev.off()

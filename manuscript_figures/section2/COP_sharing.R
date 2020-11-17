#16-Jun-2020 # Diogo Ribeiro @ UNIL

pdf("COP_sharing.pdf",20,10)

library(data.table)
require(ggplot2)
library(grid)
library(gridExtra)

inFile = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/GTEx/tissue_conservation/COP_specificity_stats.txt"

data = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

summary(data$copFreq)
summary(data$allFreq)

mean(data$copFreq)
mean(data$allFreq)

data[copFreq == 1]
data[copFreq > 1]

22621/40999

data[copFreq > 48]
#"#377eb8""#e41a1c"
g1 = ggplot( ) +
  geom_histogram( data = data, aes(x = allFreq, color = "Co-presence" ), binwidth = 1, fill = "white", alpha = 1) +
  geom_histogram( data = data, aes(x = copFreq, color = "Co-expression" ), binwidth = 1, fill = "white", alpha = 0.2) +
  # ggtitle("COP and gene pair presence across tissues") +
  xlab("# tissues") +
  scale_color_brewer(palette = "Set1") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        axis.text.x = element_text(vjust=0.6),
        legend.position = c(0.5, 0.8), legend.text=element_text(size=24), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

g2 = ggplot( ) +
  geom_histogram( data = data[allFreq > 5][copFreq == 1], aes(x = ratio*100), binwidth = 5, color = "#e41a1c", fill = "#e41a1c", alpha = 0.5) +
  geom_histogram( data = data[allFreq > 5][copFreq > 1], aes(x = ratio*100), binwidth = 5, color = "#4daf4a", fill = "#4daf4a", alpha = 0.5) +
  # geom_vline(xintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  geom_vline(xintercept = 50, linetype = "dashed") +
  annotate(geom = "text", label = "Unique (1 tissue)", x = 0, y = 7000, angle = 90, size = 7, color = "#e41a1c", fontface = "bold") +
  annotate(geom = "text", label = "Specific COPs (<15%)", x = 12, y = 7000, angle = 90, size = 7, color = "#4daf4a", fontface = "bold") +
  annotate(geom = "text", label = "Prevalent COPs (15-50%)", x = 35, y = 7000, angle = 90, size = 7, color = "#4daf4a", fontface = "bold") +
  annotate(geom = "text", label = "Conserved COPs (50-100%)", x = 80, y = 7000, angle = 90, size = 7, color = "#4daf4a", fontface = "bold") +
  # ggtitle("% tissues where COP is found") +
  xlab("% tissues") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        axis.text.x = element_text(vjust=0.6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

grid.arrange(g1,g2,nrow = 1)

dev.off()

# figure_2.R
# makes a heatmap to visualize the stoichiometric matrix of a network

library(pheatmap)

data <- read.csv("data/ab_3_2ins_5outs_full_S.csv", header = F)
invisible(png('data/figure_4b.png'))
pheatmap(t(data), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

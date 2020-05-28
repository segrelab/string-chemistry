# figure_4b.R
# make heatmap to visualize S matrix

suppressMessages(library(tidyverse))

S <- read.csv("data/ab_3_2ins_3outs_full_S.csv", header = F)
invisible(png("data/figure_4b.png"))
pheatmap(t(S), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

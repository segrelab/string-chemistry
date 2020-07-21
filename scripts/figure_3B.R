# figure_3B.R
# make heatmap to visualize S matrix

library(pheatmap)

file <- commandArgs(trailingOnly = T)[1]
S <- read.csv(file, header = F)
png("data/figure_3B.png", height = 2000, width = 2000, res = 600)
pheatmap(
  t(S), # opinions differ on whether reactions should be the columns or rows
  # don't cluster anything
  cluster_rows = F, cluster_cols = F,
  # row and column names are deafaults so they're meaningless
  show_rownames = F, show_colnames = F,
  # set legend breaks
  legend_breaks = c(-1, 0, 1, 2)
)
invisible(dev.off())

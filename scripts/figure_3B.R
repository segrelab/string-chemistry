# figure_2B.R
# make heatmap to visualize the S matrix of the network shown in figure 2A

library(pheatmap)

S <- read.csv("data/figure_2B_data.csv", header = F)
png("data/figure_2B.png", height = 2000, width = 2000, res = 600)
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

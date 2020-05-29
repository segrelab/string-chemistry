# figure_4b.R
# make heatmap to visualize S matrix

suppressMessages(library(tidyverse))

file <- commandArgs(trailingOnly = T)[1]
S <- read.csv(file, header = F)
png("data/figure_4b.png")
pheatmap(t(S), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

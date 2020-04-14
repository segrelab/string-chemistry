# figure_4bc.R
# makes heatmaps to visualize S matrices and reaction inclusion in pruned networks

library(pheatmap)
library(tidyverse)

# S matrix
S <- read.csv("data/ab_3_2ins_2outs_full_S.csv", header = F)
invisible(png("data/figure_4b.png"))
pheatmap(t(S), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

# reaction inclusion
rxn_incl <- read.csv(
  "data/ab_3_2ins_2outs_1000reps.csv",
  header = T,
  row.names = 1,
  # read the bitstring/binary vector as a character otherwise R condenses it
  # with scientific notation and loses some of the information
  colClasses = c("numeric", "character", "numeric")
) %>%
  # only need the reaction inclusion vector
  select(bitstring) %>%
  # make one column for each reaction
  separate(bitstring, into = sep = "")

# make the heatmap
invisible(png("data/figure_4c.png"))
pheatmap(rxn_incl, show_rownames = F, show_colnames = F)
invisible(dev.off())

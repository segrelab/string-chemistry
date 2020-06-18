# random_vs_min_prune.R
# makes a heatmap to visualize the differences between randomly-pruned networks
# and the min flux-pruned network for a given large network

library(pheatmap)
suppressMessages(library(tidyverse))

file <- commandArgs(trailingOnly = T)[1]

bitstring_df <- read.csv(
  file,
  header = T,
  row.names = 1,
  # read the bitstring/binary vector as a character otherwise R condenses it
  # with scientific notation and loses some of the information
  colClasses = c("numeric", "character", "numeric")
)

# make a dataframe that has one column for each reaction and is 0 or 1
rxn_incl <- as.data.frame(t(as.data.frame(strsplit(bitstring_df$bitstring, "")))) %>%
  mutate_all(function(col) as.numeric(as.character(col)))

# make the heatmap
png("data/random_vs_min_prune.png")
pheatmap(rxn_incl, show_rownames = F, show_colnames = F)
invisible(dev.off())

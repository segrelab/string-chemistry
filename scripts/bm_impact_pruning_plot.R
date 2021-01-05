# bm_impact_pruning_plot.R
# visualizes the differences in the two pruning algorithms using the data in
# bm_impact_pruning_data.csv

library(tidyverse)

theme_set(theme_bw())

data <- read.csv("data/bm_impact_pruning_data.csv")

min_data <- filter(data, type == "min")
bm_data <- filter(data, type == "bm")

ggplot() +
  geom_line(
    data = min_data, aes(x = step, y = rxn_count, alpha = trial), color = "blue"
  ) +
  geom_line(
    data = bm_data, aes(x = step, y = rxn_count, alpha = trial), color = "red"
  )

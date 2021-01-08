# bm_impact_pruning_plot.R
# visualizes the differences in the two pruning algorithms using the data in
# bm_impact_pruning_data.csv

library(tidyverse)

# set the theme once
theme_set(theme_bw())

# read in data
data <- read.csv("data/bm_impact_pruning_data.csv")

# plot number of reactions vs step for each trial and have each trial be a
# different color and each pruner be a different linetype
data %>%
  filter(trial < 20) %>%
  mutate(trial = as.factor(trial)) %>%
  ggplot(aes(x = step, y = rxn_count, color = trial, linetype = type)) +
    geom_line() +
    # remove legend for which color is which trial
    guides(color = FALSE) + 
    labs(x = "Pruning step", y = "Number of Reactions", linetype = "Pruner")

data %>%
  mutate(trial = as.factor(trial)) %>%
  ggplot(aes(x = step, y = jaccard, color = trial)) + geom_line() +
    guides(color = FALSE)

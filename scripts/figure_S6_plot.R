# figure_S6_plot.R
# visualizes the differences in the two pruning algorithms using the data in
# figure_S6_data.csv

suppressMessages(library(tidyverse))

# set the theme once
theme_set(theme_bw())

# read in data
data <- read.csv("data/figure_S6_data.csv")

# plot number of reactions vs step for each trial and have each trial be a
# different color and each pruner be a different linetype
png("data/figure_S6.png", height = 6000, width = 6700, res = 600)
data %>%
  mutate(type = ifelse(type == "bm", "Biomass-Focused", "Flux-Focused")) %>%
  ggplot(aes(x = step, y = rxn_count, col = type)) +
    geom_line() + facet_wrap(.~trial, ncol = 10, nrow = 10) +
    labs(x = "Pruning step", y = "Number of Reactions", col = "Pruner") +
    # remove labels for individual facets
    theme(strip.background = element_blank(), strip.text.x = element_blank())
invisible(dev.off())

# plot distribution of Jaccard similarities
png("data/pruning_jaccards.png", height = 4000, width = 6000, res = 600)
data %>%
  ggplot(aes(x = jaccard)) + 
  geom_histogram(fill = "royalblue", binwidth = 0.01) +
  labs(
    x = "Jaccard Similarity Between Pruned Networks",
    y = "Number of Network Pairs"
  )
invisible(dev.off())

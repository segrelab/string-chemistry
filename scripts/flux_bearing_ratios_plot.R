# flux_bearing_ratios_plot.R
# make a boxplot showing the portion of reactions with fluxes after running FBA
# on a complete string chemistry network

suppressMessages(library(tidyverse))

data <- read.csv("data/flux_bearing_ratios.csv", row.names = 1) %>%
  mutate_at(vars(monos, max_len), as.factor)

png("data/flux_bearing_ratios.png")
ggplot(data, aes(x = max_len, y = nonzero_ratio)) + 
  geom_boxplot(aes(col = monos)) +
  theme_bw()
invisible(dev.off())

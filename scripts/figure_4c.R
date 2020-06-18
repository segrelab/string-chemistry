# figure_4c.R
# make a boxplot showing the portion of reactions with fluxes after running FBA
# on a complete string chemistry network

suppressMessages(library(tidyverse))

data <- read.csv("data/figure_4c_data.csv", row.names = 1) %>%
  mutate_at(vars(monos, max_len), as.factor)

png("data/figure_4c.png")
ggplot(data, aes(x = max_len, y = nonzero_ratio)) + 
  geom_boxplot(aes(col = monos)) +
  theme_bw()
invisible(dev.off())

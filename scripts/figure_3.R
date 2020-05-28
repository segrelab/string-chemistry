# figure_3.R
# make the plots for figure 3 that show sizes of string chemistry networks
# as functions of the number of unique monomers and the maximum string length
# used to generate them

# load package(s)
suppressMessages(library(tidyverse))

# make all plots be nice
theme_set(theme_bw())

# read in data
counts <- read.csv(
  "data/figure_3_data.csv", header = F,
  col.names = c("monos", "max_len", "mets", "rxns")
) %>%
  mutate(monos = as.factor(monos)) %>%
  filter(rxns > 0) %>%
  mutate(ratio = rxns/mets)

# iJO1336
ecoli_mets <- 1805
ecoli_rxns <- 2583
# Recon3D
human_mets <- 5835
human_rxns <- 10600
# Yeast8
yeast_mets <- 2666
yeast_rxns <- 3895

# metabolite count as a function of monomers and max length
invisible(png("data/figure_3_met.png", width = 700, height = 500))
counts %>%
  ggplot(aes(x = max_len, y = mets)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_mets)) +
    annotate("text", x = 3, y = 1300, label = "E. coli") +
    geom_hline(aes(yintercept = human_mets)) +
    annotate("text", x = 3, y = 8000, label = "Human") +
    geom_hline(aes(yintercept = yeast_mets)) +
    annotate("text", x = 3, y = 3500, label = "S. cerevisiae") +
    labs(
      x = "Maximum String Length", y = "Number of Metabolites", 
      color = "Types of Monomers", 
      title = "Network Sizes By Metabolite Count"
    )
invisible(dev.off())

# reaction count as a function of monomers and max length
invisible(png("data/figure_3_rxn.png", width = 700, height = 500))
counts %>%
  ggplot(aes(x = max_len, y = rxns)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_rxns)) +
    annotate("text", x = 2.5, y = 1900, label = "E. coli") +
    geom_hline(aes(yintercept = human_rxns)) +
    annotate("text", x = 2.5, y = 15000, label = "Human") +
    geom_hline(aes(yintercept = yeast_rxns)) +
    annotate("text", x = 2.5, y = 5500, label = "S. cerevisiae") +
    labs(
      x = "Maximum String Length", y = "Number of Reactions", 
      color = "Types of Monomers", 
      title = "Network Sizes By Reaction Count"
    )
invisible(dev.off())

# reaction-to-metabolite ratio as a function of monomers and max length
invisible(png("data/figure_3_ratio.png", width = 700, height = 500))
ggplot() +
  # plot ratios from universal networks
  geom_line(data = counts, aes(x = max_len, y = ratio, color = monos)) +
  # set x-axis tick spacing
  scale_x_continuous(breaks = c(2,4,6,8,10)) +
  # add horizontal lines for real organisms
  geom_hline(aes(yintercept = ecoli_rxns/ecoli_mets)) +
  annotate("text", x = 8, y = 1.3, label = "E. coli") +
  geom_hline(aes(yintercept = human_rxns/human_mets)) +
  annotate("text", x = 8, y = 2, label = "Human") +
  geom_hline(aes(yintercept = yeast_rxns/yeast_mets)) +
  annotate("text", x = 8, y = 1.6, label = "S. cerevisiae") +
  labs(
    x = "Maximum String Length", 
    y = "Number of Reactions / Number of Metabolites", 
    color = "Types of Monomers", 
    title = "Network Sizes By Ratio of Reaction Count to Metabolite Count"
  )
invisible(dev.off())

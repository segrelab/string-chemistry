# figure_2.R
'''
Make a lineplot showing the reaction-to-metabolite ratios for several string
chemistry networks and compare them to the ratios for real metabolic networks
and pruned string chemistry networks
'''

# load package(s)
library(ggpubr)
suppressMessages(library(tidyverse))

# make all plots be nice
theme_set(theme_bw())

# read in data
counts <- read.csv(
  "data/figure_2_data.csv", header = F,
  col.names = c("monos", "max_len", "mets", "rxns")
) %>%
  mutate(monos = as.factor(monos)) %>%
  filter(rxns > 0) %>%
  mutate(ratio = rxns/mets)

pruned_counts <- read.csv(
  "data/varied_ab_5_2to5ins_2to10outs_exp.csv", 
  header = F,
  col.names = c("ins", "outs", "ratio", "pruned_pct")
)

# iJO1336
ecoli_mets <- 1805
ecoli_rxns <- 2583
# Recon3D
human_mets <- 5835
human_rxns <- 10600
# Yeast8
yeast_mets <- 2666
yeast_rxns <- 3895


png("data/figure_S2.png", height = 8000, width = 5000, res = 600)
# reaction-to-metabolite ratio as a function of monomers and max length
ggplot() +
  # plot ratios from universal networks
  geom_line(data = counts, aes(x = max_len, y = ratio, color = monos)) +
  # set x-axis tick spacing
  scale_x_continuous(breaks = c(2,4,6,8,10)) +
  # add horizontal lines for real organisms
  geom_hline(aes(yintercept = ecoli_rxns/ecoli_mets)) +
  annotate("text", x = 8, y = 1.2, label = "E. coli") +
  geom_hline(aes(yintercept = human_rxns/human_mets)) +
  annotate("text", x = 8, y = 2, label = "Human") +
  geom_hline(aes(yintercept = yeast_rxns/yeast_mets)) +
  annotate("text", x = 8, y = 1.62, label = "S. cerevisiae") +
  # add points for pruned networks
  geom_boxplot(data = pruned_counts, aes(x = 5, y = ratio, color = "2")) +
  annotate("text", x = 6.2, y = 1, label = "Pruned networks") +
  # add an arrow to this boxplot from the A = 2 line
  geom_segment(
    aes(x = 5, y = 3.1, xend = 5, yend = 1.4, color = "2"), 
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  # add labels
  labs(
    x = "Maximum String Length", 
    y = "Reactions / Metabolites", 
    color = "Types of Monomers", 
    title = "Network Sizes By Ratio of Reaction Count to Metabolite Count"
  )
invisible(dev.off())

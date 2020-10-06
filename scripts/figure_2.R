# figure_2.R
# make the plots for figure 3 that show sizes of string chemistry networks
# as functions of the number of unique monomers and the maximum string length
# used to generate them

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

# metabolite count as a function of monomers and max length
panel_a <- counts %>%
  ggplot(aes(x = max_len, y = mets)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_mets)) +
    annotate("text", x = 3, y = 1200, label = "E. coli") +
    geom_hline(aes(yintercept = human_mets)) +
    annotate("text", x = 3, y = 8500, label = "Human") +
    geom_hline(aes(yintercept = yeast_mets)) +
    annotate("text", x = 3, y = 3700, label = "S. cerevisiae") +
    labs(
      x = "Maximum String Length", y = "Metabolites", 
      color = "Types of Monomers", 
      title = "Network Sizes By Metabolite Count"
    )


# reaction count as a function of monomers and max length
panel_b <- counts %>%
  ggplot(aes(x = max_len, y = rxns)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_rxns)) +
    annotate("text", x = 2.5, y = 1700, label = "E. coli") +
    geom_hline(aes(yintercept = human_rxns)) +
    annotate("text", x = 2.5, y = 16000, label = "Human") +
    geom_hline(aes(yintercept = yeast_rxns)) +
    annotate("text", x = 2.5, y = 6000, label = "S. cerevisiae") +
    labs(
      x = "Maximum String Length", y = "Reactions", 
      color = "Types of Monomers", 
      title = "Network Sizes By Reaction Count"
    )

# reaction-to-metabolite ratio as a function of monomers and max length
panel_c <- ggplot() +
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

png("data/figure_2.png", height = 8000, width = 5000, res = 600)
ggarrange(panel_a, panel_b, panel_c, nrow = 3, labels = c("a", "b", "c"))
invisible(dev.off())

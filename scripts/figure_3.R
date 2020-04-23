# figure_3.R
# make the plots for figure 3 that show sizes of string chemistry networks
# as functions of the number of unique monomers and the maximum string length
# used to generate them

# load package(s)
library(tidyverse)

# make all plots be nice
theme_set(theme_bw())

# read in data
counts <- read.csv(
  "../data/figure_3_data.csv", header = F,
  col.names = c("monos", "max_len", "mets", "rxns")
) %>%
  mutate(monos = as.factor(monos)) %>%
  filter(rxns > 0) %>%
  mutate(ratio = rxns/mets)

# find all ratios files
all_files <- list.files("../data")
ratio_files <- all_files[grepl("ratios.csv", all_files)]

# parse information from all of those files into a nice dataframe
parse_ratio_file <- function(ratio_file) {
  # get mean and standard deviation of reaction-to-metabolite ratios
  ratio_dist <- read.csv(paste("../data/", ratio_file, sep = ""), header = F)
  mean_ratio <- mean(ratio_dist$V1)
  ratio_sd <- sd(ratio_dist$V1)
  # get number of monomers and max length used to generate the universal
  # network all these ratios correspond to
  filename_bits <- strsplit(ratio_file, "_")[[1]]
  monos <- filename_bits[1]
  max_len <- filename_bits[2]
  return(c(monos, max_len, mean_ratio, ratio_sd))
}

ratio_df <- as.data.frame(t(sapply(ratio_files, parse_ratio_file)))
# make the dataframe a bit easier to work with
colnames(ratio_df) <- c("monos", "max_len", "mean_ratio", "ratio_sd")
ratio_df <- ratio_df %>%
  mutate_at(vars(-monos), function(col) as.numeric(as.character(col))) %>%
  # for plotting purposes we only need to know how many monomers there were,
  # not exactly what they were
  mutate(monos = nchar(as.character(monos))) %>%
  # since we're using monos to separate lines in the plots, make it a factor
  mutate(monos = as.factor(monos))

# googled these so they're hard-coded
ecoli_mets <- 1805
ecoli_rxns <- 2583
kegg_mets <- 18686
kegg_rxns <- 11395

# metabolite count as a function of monomers and max length
png("data/figure_3_met.png")
counts %>%
  ggplot(aes(x = max_len, y = mets)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_mets)) +
    annotate("text", x = 3, y = 2700, label = "E. coli Metabolites") +
    geom_hline(aes(yintercept = kegg_mets)) +
    annotate("text", x = 3, y = 30000, label = "KEGG Metabolites") +
    labs(
      x = "Maximum String Length", y = "Number of Metabolites", 
      color = "Types of Monomers", 
      title = "Network Sizes By Metabolite Count"
    )
dev.off()

# reaction count as a function of monomers and max length
png("data/figure_3_rxn.png")
counts %>%
  ggplot(aes(x = max_len, y = rxns)) + 
    geom_line(aes(color = monos)) + 
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    scale_y_continuous(trans = "log10") +
    geom_hline(aes(yintercept = ecoli_rxns)) +
    annotate("text", x = 3, y = 4000, label = "E. coli Reactions") +
    geom_hline(aes(yintercept = kegg_rxns)) +
    annotate("text", x = 3, y = 20000, label = "KEGG Reactions") +
    labs(
      x = "Maximum String Length", y = "Number of Reactions", 
      color = "Types of Monomers", 
      title = "Network Sizes By Reaction Count"
    )
dev.off()

# reaction-to-metabolite ratio as a function of monomers and max length
png("data/figure_3_ratio.png")
ggplot() +
  # plot ratios from universal networks
  geom_line(data = counts, aes(x = max_len, y = ratio, color = monos)) +
  # plot ratios from pruned networks
  geom_point(data = ratio_df, aes(x = max_len, y = mean_ratio, color = monos)) +
  geom_errorbar(
    data = ratio_df, 
    aes(
      x = max_len, y = mean_ratio, color = monos,
      ymin = mean_ratio-ratio_sd, ymax = mean_ratio+ratio_sd
    )
  ) +
  # set x-axis tick spacing
  scale_x_continuous(breaks = c(2,4,6,8,10)) +
  # add horizontal line for E. coli ratio
  geom_hline(aes(yintercept = ecoli_rxns/ecoli_mets)) +
  annotate("text", x = 8, y = 1.75, label = "E. coli Ratio") +
  # add horizontal line for KEGG ratio
  geom_hline(aes(yintercept = kegg_rxns/kegg_mets)) +
  annotate("text", x = 8, y = 0.9, label = "KEGG Ratio") +
  labs(
    x = "Maximum String Length", 
    y = "Number of Reactions / Number of Metabolites", 
    color = "Types of Monomers", 
    title = "Network Sizes By Ratio of Reaction Count to Metabolite Count"
  )
dev.off()

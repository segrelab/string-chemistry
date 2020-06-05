# figure_4cd.R
# make reaction-inclusion heatmaps for the pruned networks produced by figure_4a
suppressMessages(library(tidyverse))

file <- commandArgs(trailingOnly = T)[1]
# R makes the reaction-inclusion vectors into scientific notation if you don't
# read them as character vectors
raw_data <- read.csv(file, colClasses = "character")

# have to do networks with export reactions separately from networks without
raw_exp_data <- raw_data %>%
  filter(export == "export")
raw_no_exp_data <- raw_data %>%
  filter(export == "no export")

exp_data <- raw_exp_data %>%
  # make a column for each reaction
  as.data.frame(t(as.data.frame(strsplit(.$rxn_incl, "")))) %>%
  mutate_all(function(col) as.numeric(as.character(col))) %>%
  cbind(., select(raw_exp_data, -X, -rxn_incl))

no_exp_data <- raw_no_exp_data %>%
  # make a column for each reaction
  as.data.frame(t(as.data.frame(strsplit(.$rxn_incl, "")))) %>%
  mutate_all(function(col) as.numeric(as.character(col))) %>%
  cbind(., select(raw_no_exp_data, -X, -rxn_incl))

png("data/figure_4c.png")
pheatmap(exp_data, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

png("data/figure_4d.png")
pheatmap(no_exp_data, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
invisible(dev.off())

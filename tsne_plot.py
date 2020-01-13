# tsne_plot.py
# make a t-SNE plot to compare bitstrings from many different networks

import sys
import pandas as pd
from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# read in dataframe of bitstrings
filename = sys.argv[1]
bitstring_df = pd.read_csv(filename) #, sep = '\t')
# want each reaction bit in its own column and need to not have the rxn_count
# column in the dataframe we pass to the PCA function
bits_as_cols = bitstring_df['bitstring'].apply(lambda x: pd.Series(list(x)))
# first and last columns will be empty and all columns will be strings
tsne_ready = bits_as_cols.iloc[:,1:-1].astype('int32')

# actually do t-SNE
tsne = TSNE(n_components = 2, random_state = 0)
tsne_results = tsne.fit_transform(tsne_ready)
tsne_df = pd.DataFrame(data = tsne_results, columns = ['x', 'y'])
# add in other info for plotting purposes
plotting_df = pd.concat([tsne_df, bitstring_df], axis = 1)

# make a scatterplot

# if coloring by unique biomass reactions, make a colormap
#cdict = {v: k for k, v in enumerate(np.unique(bitstring_df.biomass))}
#cvals = [cdict[c] for c in bitstring_df.biomass]

# if coloring by input metabolites, make a colormap
# start by getting the column of input metabolites, splitting it into a list of
# lists, then pasting the sublists together so that there's one string to use
# for making the colormap
in_groups = ['-'.join(sorted(ins.split('-'))) for ins in bitstring_df.inputs]
cdict = {v: k for k, v in enumerate(np.unique(in_groups))}
cvals = [cdict[c] for c in in_groups]

# if sizing by number of occurrences, multiply values by 10
#plotting_df.occurrences = [10*x for x in plotting_df.occurrences]

# make the figure large
plt.figure(figsize = (20,22))

# make the text legible
matplotlib.rcParams.update({
    'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
    'axes.labelsize': 18
})

plt.scatter(
    plotting_df.x, plotting_df.y,
#    s = plotting_df.occurrences, # size according to number of observations
#    c = plotting_df.rxn_count, # color according to number of reactions
    c = cvals, cmap = 'nipy_spectral',  # color by biomass reaction
    s = 100 # make points big
)
plt.title('20 Different Biomass Reactions Pruned on 1000 2-Metabolite Environments Each')
plt.show()

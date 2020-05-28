# figure_6bc.py
# make a UMAP visualizing the results of multiple_env_min_prune.py

import sys
import pandas as pd
import umap
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

# process file containing reaction inclusion bitstrings so we can do UMAP
filename = sys.argv[1]
bitstring_df = pd.read_csv(filename)
# want each reaction bit in its own column and only pass those columns to umap
bits_as_cols = bitstring_df['bitstring'].apply(lambda x: pd.Series(list(x)))
# first and last columns will be empty and all columns will be strings
umap_ready = bits_as_cols.iloc[:,1:-1].astype('int32')

# do UMAP
reducer = umap.UMAP()
umap_results = reducer.fit_transform(umap_ready)
umap_df = pd.DataFrame(data = umap_results, columns = ['x', 'y'])
# add in other info for plotting purposes
plotting_df = pd.concat([umap_df, bitstring_df], axis = 1)

# make a scatterplot

# make a colormap for biomass reactions
bm_cdict = {v: k for k, v in enumerate(np.unique(bitstring_df.biomass))}
bm_cvals = [bm_cdict[c] for c in bitstring_df.biomass]

# make a colormap for input metabolites
# start by getting the column of input metabolites, splitting it into a list of
# lists, then pasting the sublists together so that there's one string to use
# for making the colormap
in_groups = ['-'.join(sorted(ins.split('-'))) for ins in bitstring_df.inputs]
in_cdict = {v: k for k, v in enumerate(np.unique(in_groups))}
in_cvals = [in_cdict[c] for c in in_groups]

# make the figure large
plt.figure(figsize = (8,7))

# make the text legible
matplotlib.rcParams.update({
    'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
    'axes.labelsize': 18
})

# do one scatterplot colored by biomass reactions
plt.scatter(
    plotting_df.x, plotting_df.y,
    c = bm_cvals, cmap = 'nipy_spectral',
    s = 10
)
#plt.title('10 Different Biomass Reactions Pruned on 1000 2-Metabolite Environments Each')
plt.xlabel('UMAP_1')
plt.ylabel('UMAP_2')

# save the figure
plt.savefig('data/figure_6b.png', dpi = 600)

# and one colored by input metabolites
plt.figure(2)
plt.figure(figsize = (8,7))
plt.scatter(
    plotting_df.x, plotting_df.y,
    c = in_cvals, cmap = 'nipy_spectral',
    s = 10
)
#plt.title('10 Different Biomass Reactions Pruned on 1000 2-Metabolite Environments Each')
plt.xlabel('UMAP_1')
plt.ylabel('UMAP_2')

# save the figure
plt.savefig('data/figure_6c.png', dpi = 600)

# pruning_comparison_pca.py
# make some plots to compare the reaction inclusion bitstrings of several
# networks

import sys
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# get reaction inclusion bitstrings
filename = sys.argv[1]
bitstring_df = pd.read_csv(filename) #, sep = '\t')
# want each reaction bit in its own column and need to not have the rxn_count
# column in the dataframe we pass to the PCA function
bits_as_cols = bitstring_df['bitstring'].apply(lambda x: pd.Series(list(x)))
# first and last columns will be empty and all columns will be strings
pca_ready = bits_as_cols.iloc[:,1:-1].astype('int32')

# do PCA with scikit-learn
pca = PCA(n_components = 2)
pcs = pca.fit_transform(pca_ready)
pc_df = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2'])
pca_results = pd.concat([pc_df, bitstring_df], axis = 1)

# now make a scatterplot

# if sizing by number of occurrences, multiply all values by 10 so the points
# are visible
#pca_results.occurrences = [10*x for x in pca_results.occurrences]

# if coloring by unique biomass reactions, make a colormap
#cdict = {v: k for k, v in enumerate(np.unique(pca_results.biomass))}
#cvals = [cdict[c] for c in pca_results.biomass]

# if coloring by input metabolites, make a colormap
# start by getting the column of input metabolites, splitting it into a list of
# lists, then pasting the sublists together so that there's one string to use
# for making the colormap
in_groups = ['-'.join(sorted(ins.split('-'))) for ins in pca_results.inputs]
cdict = {v: k for k, v in enumerate(np.unique(in_groups))}
cvals = [cdict[c] for c in in_groups]

# make the figure large
plt.figure(figsize = (20,22))

# make the text legible
matplotlib.rcParams.update({
    'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
    'axes.labelsize': 18
})

# actually make the plot
plt.scatter(
    pca_results.PC1, pca_results.PC2,
#    s = pca_results.occurrences, # size according to number of observations
#    c = pca_results.rxn_count, # color according to number of reactions
    c = cvals, cmap = 'gist_ncar',  # color by biomass reaction
    s = 100 # make points big
)
plt.title('20 Biomass Reactions Pruned on 100 2-Metabolite Environments Each')
plt.show()

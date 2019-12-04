# pca_plot.py
# make some plots to compare the results of multiple rounds of random pruning
# of a network

import sys
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# get reaction inclusion bitstrings from output of pruning.py
filename = sys.argv[1]
bitstring_df = pd.read_csv(filename)
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
#pca_results.occurences = [10*x for x in pca_results.occurences]
plt.figure(figsize = (20,22))
plt.scatter(
    pca_results.PC1, pca_results.PC2,
#    s = pca_results.occurences, # size according to number of observations
#    c = pca_results.rxn_count # color according to number of reactions
    c = pca_results.biomass # color according to biomass reaction
)
#plt.colorbar()
plt.show()

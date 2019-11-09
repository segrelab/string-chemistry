# pruning_plots.py
# make some plots to compare the results of multiple rounds of random pruning
# of a network

import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# get reaction inclusion bitstrings from output of pruning.py
bitstring_df = pd.read_csv('bitstrings.csv')
# want each reaction bit in its own column and need to not have the rxn_count
# column in the dataframe we pass to the PCA function
bits_as_cols = bitstring_df['bitstring'].str.split('', expand = True)
# first and last columns will be empty and all columns will be strings
pca_ready = bits_as_cols.iloc[:,1:-1].astype('int32')
# do PCA with scikit-learn
pca = PCA(n_components = 2)
pcs = pca.fit_transform(pca_ready)
pc_df = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2'])
pca_results = pd.concat([
    pc_df, 
    bitstring_df[['bitstring', 'occurences']]],
    axis = 1
)
# now make a scatterplot
plt.scatter(pca_results.PC1, pca_results.PC2, s = pca_results.occurences)
plt.show()

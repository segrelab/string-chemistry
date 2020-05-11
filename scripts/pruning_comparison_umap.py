# figure_6bc.py
# make a UMAP to compare some reaction bitstrings

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

# set figure size
plt.figure(figsize = (5,6))

# make the text legible
matplotlib.rcParams.update({
    'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
    'axes.labelsize': 18
})

plt.scatter(plotting_df.x, plotting_df.y, s = 10)
#plt.title('10 Different Biomass Reactions Pruned on 1000 2-Metabolite Environments Each')
plt.xlabel('UMAP_1')
plt.ylabel('UMAP_2')

plt.show()

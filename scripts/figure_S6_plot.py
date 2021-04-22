# figure_S6_plot.py
'''
Make several UMAPs visualizing the results of figure_S6_data.py
'''

import pandas as pd
import umap
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

# just need name of file with output
def do_umap(filename):
    # get the reaction-inclusion vector out of the input file and make it into
    # a bunch of 1/0 columns instead of one column of strings of 1s and 0s
    data = pd.read_csv(filename)
    # want each reaction bit in its own column and only pass those columns to umap
    rxn_incl_cols = data['rxn_incl'].apply(lambda x: pd.Series(list(x)))
    # first and last columns will be empty and all columns will be strings
    umap_ready = rxn_incl_cols.iloc[:,1:-1].astype('int32')
    # do UMAP
    reducer = umap.UMAP()
    umap_results = reducer.fit_transform(umap_ready)
    umap_df = pd.DataFrame(data = umap_results, columns = ['x', 'y'])
    # add in other info for plotting purposes
    plotting_df = pd.concat([umap_df, data], axis = 1)
    return(plotting_df)

# make a scatterplot of the UMAP results where each point is colored according
# to the biomass reaction used by the associated network
def biomass_plot(axes, data, label):
    # make a colormap
    bm_cdict = {v: k for k, v in enumerate(np.unique(data.biomass))}
    bm_cvals = [bm_cdict[c] for c in data.biomass]
    # make the scatterplot
    axes.scatter(
        data.x, data.y,
        c = bm_cvals, cmap = 'nipy_spectral',
        # points need to be as small as possible so they don't overlap
        s = 1
    )
    axes.set_xlabel('UMAP_1')
    axes.set_ylabel('UMAP_2')
    # add the panel label
    axes.text(
        # position the letter above the plot on the left side
        -0.11, 1.1, label, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = axes.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    # don't need to return anything because matplotlib is a bit weird

# make a scatterplot of the UMAP results where each point is colored according
# to the nutrients available to the associated network
def nutrient_plot(axes, data, label):
    # make a colormap
    # start by getting the column of input metabolites, splitting it into a list of
    # lists, then pasting the sublists together so that there's one string to use
    # for making the colormap
    in_groups = ['-'.join(sorted(ins.split('-'))) for ins in data.env]
    in_cdict = {v: k for k, v in enumerate(np.unique(in_groups))}
    in_cvals = [in_cdict[c] for c in in_groups]
    # make the scatterplot
    axes.scatter(
        data.x, data.y,
        c = in_cvals, cmap = 'nipy_spectral',
        # points need to be as small as possible so they don't overlap
        s = 1
    )
    axes.set_xlabel('UMAP_1')
    axes.set_ylabel('UMAP_2')
    # add the panel label
    axes.text(
        # position the letter above the plot on the left side
        -0.11, 1.1, label, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = axes.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    # return nothing because matplotlib

# make a scatterplot of UMAP results where each point is colored according to
# the growth rate of the associated network
def growth_plot(figure, axes, data, label):
    plot = axes.scatter(
        data.x, data.y, 
        c = data.growth, cmap = 'Blues', 
        # points need to be as small as possible so they don't overlap
        s = 1
    )
    figure.colorbar(plot, ax = axes)
    axes.set_xlabel('UMAP_1')
    axes.set_ylabel('UMAP_2')
    # add the panel label
    axes.text(
        # position the letter above the plot on the left side
        -0.11, 1.1, label, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = axes.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    # return nothing because matplotlib

# get UMAP results for networks with and without export reactions
print('Doing UMAP')
umap_out = do_umap('data/figure_S6_data.csv')

print('Plotting UMAP results')
# set up the subplots in a row
fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (12,3))

# panels B-D use the data from the networks with export reactions
biomass_plot(ax[0], umap_out, 'a')
nutrient_plot(ax[1], umap_out, 'b')
growth_plot(fig, ax[2], umap_out, 'c')

# tight_layout just fixes all sorts of problems with subplots overlapping
plt.tight_layout()
plt.savefig('data/figure_S6.png', dpi = 600)
#plt.show()

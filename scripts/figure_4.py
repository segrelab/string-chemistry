# figure_4.py
# make several UMAPs visualizing the results of multiple_env_min_prune.py

import pandas as pd
import umap
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

# need name of file with multiple_env_prune output, name for the plot colored
# by biomass reaction, name for the plot colored by environment and name for
# plot colored by growth rates
def do_umap(filename, bm_plot_name, env_plot_name, growth_plot_name):
    # get the reaction-inclusion vector out of the input file and make it into
    # a bunch of 1/0 columns instead of one column of strings of 1s and 0s
    data = pd.read_csv(filename)
    # want each reaction bit in its own column and only pass those columns to umap
    rxn_incl_cols = data['rxn_incl'].apply(lambda x: pd.Series(list(x)))
    # first and last columns will be empty and all columns will be strings
    umap_ready = rxn_incl_cols.iloc[:,1:-1].astype('int32')

    # do UMAP
    print('Doing UMAP')
    reducer = umap.UMAP()
    umap_results = reducer.fit_transform(umap_ready)
    umap_df = pd.DataFrame(data = umap_results, columns = ['x', 'y'])
    # add in other info for plotting purposes
    plotting_df = pd.concat([umap_df, data], axis = 1)

    # make a colormap for biomass reactions
    print('Plotting')
    bm_cdict = {v: k for k, v in enumerate(np.unique(data.biomass))}
    bm_cvals = [bm_cdict[c] for c in data.biomass]

    # make a colormap for input metabolites
    # start by getting the column of input metabolites, splitting it into a list of
    # lists, then pasting the sublists together so that there's one string to use
    # for making the colormap
    in_groups = ['-'.join(sorted(ins.split('-'))) for ins in data.env]
    in_cdict = {v: k for k, v in enumerate(np.unique(in_groups))}
    in_cvals = [in_cdict[c] for c in in_groups]

    # make the text legible
    matplotlib.rcParams.update({
        'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
        'axes.labelsize': 18
    })

    # do one scatterplot colored by biomass reactions
    # make the figure large
    plt.figure(figsize = (8,7))
    plt.scatter(
        plotting_df.x, plotting_df.y,
        c = bm_cvals, cmap = 'nipy_spectral',
        s = 10
    )
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.savefig(f'data/{bm_plot_name}.png', dpi = 600)

    # do one scatterplot colored by environments
    plt.figure(2)
    plt.figure(figsize = (8,7))
    plt.scatter(
        plotting_df.x, plotting_df.y,
        c = in_cvals, cmap = 'nipy_spectral',
        s = 10
    )
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.savefig(f'data/{env_plot_name}.png', dpi = 600)

    # do one scatterplot colored by growth
    plt.figure(2)
    plt.figure(figsize = (8,7))
    fig, ax = plt.subplots()
    plot = ax.scatter(
        plotting_df.x, 
        plotting_df.y, 
        c = plotting_df.growth, 
        cmap = 'Blues', 
        s = 10
    )
    fig.colorbar(plot, ax = ax)
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.savefig(f'data/{growth_plot_name}.png', dpi = 600)

    # return nothing; the plots are the only required output

# panels B and C are with export reactions, panels D and E are without
do_umap(
    'data/multiple_env_min_prune_ab_5_2ins_1000envs_5outs_10orgs_yesexp.csv',
    'figure_4B',
    'figure_4C',
    'figure_4_export_growth'
)
do_umap(
    'data/multiple_env_min_prune_ab_5_2ins_1000envs_5outs_10orgs_noexp.csv',
    'figure_4D',
    'figure_4E',
    'figure_4_noexport_growth'
)

# comparing_nets_plots.py
'''
Make some plots to compare string chemistry networks to the E. coli and yeast
metabolic networks using the data generated by comparing_nets_data.py
'''

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

def plot_deg_dist(axes, data, color, label):
    '''
    Plot a degree distribution on a log-log scatterplot given a dataframe with
    a "degree" column and a "freq" column
    '''
    axes.scatter(
        data = data, 
        x = 'degree', 
        y = 'freq', 
        c = color, 
        label = label
    )
    return(axes)

def plot_flux_dist(axes, data, color, label):
    '''
    Plot a flux distribution on a log-log scatterplot given a dataframe with a 
    "flux" column and a "freq" column
    '''
    axes.scatter(data = data, x = 'flux', y = 'freq', c = color)
    axes.set_xlabel('Flux')
    axes.set_ylabel('Frequency of Reactions With Flux')
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.text(
        # position the letter above the plot on the left side
        -0.11, 1.1, label, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = axes.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    return(axes)

# read in files produced by comparing_nets_data.py
scn_deg_dists = pd.read_csv('data/scn_deg_dists.csv')
scn_flux_dists = pd.read_csv('data/scn_flux_dists.csv')
ecoli_deg_dist = pd.read_csv('data/ecoli_deg_dist.csv')
ecoli_flux_dist = pd.read_csv('data/ecoli_flux_dist.csv')
yeast_deg_dist = pd.read_csv('data/yeast_deg_dist.csv')
yeast_flux_dist = pd.read_csv('data/yeast_flux_dist.csv')

# make one plot for all of the degree distributions and three plots for the
# flux distributions in a row below
fig = plt.figure(figsize = (9,8))
gs = GridSpec(2,2, figure = fig)
deg_ax = fig.add_subplot(gs[0,0])
scn_flux_ax = fig.add_subplot(gs[0,1])
ecoli_flux_ax = fig.add_subplot(gs[1,0])
yeast_flux_ax = fig.add_subplot(gs[1,1])

# give each degree distribution a different color since they're all going on
# the same axes
# need to plot both points and error bars for the string chemistry networks
deg_ax.scatter(
    data = scn_deg_dists, 
    x = 'degree', 
    y = 'mean_freq', 
    c = 'blue',
    label = 'String Chemistries'
)
deg_ax.errorbar(
    data = scn_deg_dists,
    x = 'degree',
    y = 'mean_freq',
    yerr = 'std_freq',
    linestyle = 'None',
    c = 'blue',
    label = ''
)
deg_ax = plot_deg_dist(deg_ax, ecoli_deg_dist, 'red', 'E. coli')
deg_ax = plot_deg_dist(deg_ax, yeast_deg_dist, 'green', 'S. cerevisiae')
# set panel and axis labels and make them log-scaled
deg_ax.set_xlabel('Degree')
deg_ax.set_ylabel('Frequency of Metabolites With Degree')
deg_ax.set_xscale('log')
deg_ax.set_yscale('log')
deg_ax.text(
    # position the letter above the plot on the left side
    -0.11, 1.1, 'a', va = 'top', ha = 'right',
    # no idea what this does; got it from StackOverflow
    transform = deg_ax.transAxes,
    # make the letter large and bold
    fontsize = 14, fontweight = 'bold' 
)
# make a legend so it's clear which colors represent which thing
deg_ax.legend()

# now do the flux distributions all on different axes
scn_flux_ax.scatter(
    data = scn_flux_dists,
    x = 'flux',
    y = 'mean_freq',
    c = 'blue'
)
scn_flux_ax.errorbar(
    data = scn_flux_dists,
    x = 'flux',
    y = 'mean_freq',
    yerr = 'std_freq',
    c = 'blue',
    linestyle = 'None'
)
scn_flux_ax.set_xlabel('Flux')
scn_flux_ax.set_ylabel('Frequency of Reactions With Flux')
scn_flux_ax.set_xscale('log')
scn_flux_ax.set_yscale('log')
scn_flux_ax.text(
    # position the letter above the plot on the left side
    -0.11, 1.1, 'b', va = 'top', ha = 'right',
    # no idea what this does; got it from StackOverflow
    transform = scn_flux_ax.transAxes,
    # make the letter large and bold
    fontsize = 14, fontweight = 'bold' 
)
plot_flux_dist(ecoli_flux_ax, ecoli_flux_dist, 'red', 'c')
plot_flux_dist(yeast_flux_ax, yeast_flux_dist, 'green', 'd')

# make sure subplots don't overlap with each other
plt.tight_layout()
#plt.savefig('data/deg_flux_dists.png', dpi = 600)
plt.show()

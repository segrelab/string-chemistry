# figure_2_plot.py
'''
Read in data in data/figure_2_data.csv and make two lineplots showing the sizes
of various string chemistry networks in comparison to the E. coli and yeast
metabolic networks, then read in data in data/comparing_nets_data.csv and make
two scatterplots showing the degree and flux distributions in several string
chemistry networks and the E. coli and yeast networks
'''

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

def make_size_plot(ax, data, y_var, y_lab, panel):
    '''
    Plot either metabolite or reaction counts of string chemistry networks as
    a lineplot with a log-scaled y axis
    '''
    ax = sns.lineplot(
        data = data, x = 'max_pol', y = y_var, 
        hue = 'monos', palette = 'viridis_r',
        ax = ax
    )
    # rename the legend from 'monos' to 'Types of Monomers'
    ax.legend(title = 'Types of Monomers')
    # make y-axis log-scaled and set axis labels
    ax.set(yscale = 'log', xlabel = 'Maximum String Length', ylabel = y_lab)
    # give this panel a label
    ax.text(
        # position the letter above the plot on the left side
        -0.08, 1.05, panel, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = ax.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    return(ax)

def plot_scn(axes, data, x_col):
    '''
    Plot the mean frequencies of either degrees or fluxes on a log-log
    scatterplot with errorbars on the points representing the standard 
    deviations of those frequencies
    '''
    axes.scatter(
        data = data, 
        x = x_col, 
        y = 'mean_freq', 
        c = 'blue',
        label = 'String Chemistries'
    )
    axes.errorbar(
        data = data,
        x = x_col,
        y = 'mean_freq',
        yerr = 'std_freq',
        linestyle = 'None',
        c = 'blue',
        label = ''
    )
    return(axes)

def plot_dist(axes, data, x_col, color, label):
    '''
    Plot a distribution of either degrees or fluxes on a log-log scatterplot
    '''
    axes.scatter(
        data = data, 
        x = x_col, 
        y = 'freq', 
        c = color, 
        label = label
    )
    return(axes)

def fix_axes(axes, x_lab, y_lab, panel):
    '''
    Make a set of axes log-scaled with appropriate labels
    '''
    axes.set_xlabel(x_lab)
    axes.set_ylabel(f'Frequency of {y_lab} With {x_lab}')
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.text(
        # position the letter above the plot on the left side
        -0.08, 1.05, panel, va = 'top', ha = 'right',
        # no idea what this does; got it from StackOverflow
        transform = axes.transAxes,
        # make the letter large and bold
        fontsize = 14, fontweight = 'bold' 
    )
    # make a legend so it's clear which colors represent which thing
    axes.legend()
    return(axes)

# read in data
# read data
data = pd.read_csv(
    'data/figure_2_data.csv',
    names = ['monos', 'max_pol', 'met_count', 'rxn_count']
)
scn_deg_dists = pd.read_csv('data/scn_deg_dists.csv')
scn_flux_dists = pd.read_csv('data/scn_flux_dists.csv')
ecoli_deg_dist = pd.read_csv('data/ecoli_deg_dist.csv')
ecoli_flux_dist = pd.read_csv('data/ecoli_flux_dist.csv')
yeast_deg_dist = pd.read_csv('data/yeast_deg_dist.csv')
yeast_flux_dist = pd.read_csv('data/yeast_flux_dist.csv')

# set up a two by two grid of subplots
fig = plt.figure(figsize = (10,10))
gs = GridSpec(2,2, figure = fig)
met_ax = fig.add_subplot(gs[0,0])
rxn_ax = fig.add_subplot(gs[0,1])
deg_ax = fig.add_subplot(gs[1,0])
flux_ax = fig.add_subplot(gs[1,1])

# plot metabolite and reaction counts as lineplots with log-scaled y-axes
met_ax = make_size_plot(met_ax, data, 'met_count', 'Number of Metabolites', 'a')
rxn_ax = make_size_plot(rxn_ax, data, 'rxn_count', 'Number of Reactions', 'b')

# give each degree distribution a different color since they're all going on
# the same axes
deg_ax = plot_scn(deg_ax, scn_deg_dists, 'degree')
deg_ax = plot_dist(deg_ax, ecoli_deg_dist, 'degree', 'red', 'E. coli')
deg_ax = plot_dist(
    deg_ax, yeast_deg_dist, 'degree', 'green', 'S. cerevisiae'
)
deg_ax = fix_axes(deg_ax, 'Degree', 'Metabolites', 'c')

# now do the same thing with the flux distributions
flux_ax = plot_scn(flux_ax, scn_flux_dists, 'flux')
flux_ax = plot_dist(flux_ax, ecoli_flux_dist, 'flux', 'red', 'E. coli')
flux_ax = plot_dist(flux_ax, yeast_flux_dist, 'flux', 'green', 'S. cerevisiae')
flux_ax = fix_axes(flux_ax, 'Normalized Flux', 'Reactions', 'd')
# make sure subplots don't overlap with each other
plt.tight_layout()
plt.show()

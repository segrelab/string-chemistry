# comparing_nets.py
'''
Make some string chemistry networks of various sizes and load some real 
metabolic networks and make some plots that compare various characteristics
'''

import string_chem_net as scn
import cobra
import re
import networkx as nx
from networkx.algorithms import bipartite
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def get_edgelist_from_cobra(model):
    '''
    Given a COBRApy model, make an edgelist of the graph where all of the
    metabolites are joined by an edge if they are a product-reactant pair
    '''
    edgelist = list()
    for r in model.reactions:
        # split the reaction into products and reactants
        (front, middle, end) = re.split('(\-\->|<=>|<\-\-)', r.reaction)
        # split on + to get individual metabolites and strip extra spaces and
        # drop stoichiometric coefficients if present
        reacs = [m.strip(' ').split(' ')[-1] for m in front.split('+')]
        prods = [m.strip(' ').split(' ')[-1] for m in end.split('+')]
        # add edges to the reaction from all reactants
        for m in reacs:
            # skip empty strings; they come from exchange reactions
            if m == '':
                pass
            else:
                edgelist.append((m, r.id))
        # do the same for products
        for m in prods:
            if m == '':
                pass
            else:
                edgelist.append((r.id, m))
    return(edgelist)

def plot_deg_dist(G, ax, title):
    '''
    Plot the degree distribution of a graph of a metabolic network but only
    plot degrees of metabolite nodes
    '''
    # first, use networkx magic to identify all the metabolite nodes
    node_sets = bipartite.sets(G)
    # the metabolite nodes will be in the smaller set since there are always
    # fewer metabolites than reactions
    met_nodes = list()
    if len(node_sets[0]) < len(node_sets[1]):
        met_nodes = node_sets[0]
    else:
        met_nodes = node_sets[1]
    # get degrees of only the metabolite nodes
    degs = [G.degree(n) for n in G.nodes() if n in met_nodes]
    # find out what proportion of nodes have each degree; this takes several
    # steps either because pandas is annoying or I'm not good at using it
    deg_df = pd.DataFrame(pd.Series(degs).value_counts(), columns = ['count'])
    deg_df.index.name = 'degree'
    deg_df = deg_df.reset_index()
    deg_df['prop'] = deg_df['count'] / deg_df['count'].sum()
    # make a scatterplot with a log-scaled y-axis
    ax = sns.scatterplot(data = deg_df, x = 'degree', y = 'prop', ax = ax)
    ax.set(
        xlabel = 'Degree',
        ylabel = 'Proportion of Nodes',
        yscale = 'log',
        title = title
    )
    return(ax)

# generate some string chemistry networks
sc_4_5_scn = scn.CreateNetwork('abcd', 5)
sc_5_4_scn = scn.CreateNetwork('abcde', 4)
sc_3_7_scn = scn.CreateNetwork('abc', 7)

# get edgelists for these networks
sc_4_5_edges = scn.make_edgelist(sc_4_5_scn.rxn_list)
sc_5_4_edges = scn.make_edgelist(sc_5_4_scn.rxn_list)
sc_3_7_edges = scn.make_edgelist(sc_3_7_scn.rxn_list)

# read in the real metabolic networks using COBRApy
ecoli = cobra.io.read_sbml_model('data/iJO1366.xml')
yeast = cobra.io.read_sbml_model('data/yeastGEM.xml')

# make edgelists from those models
ecoli_edges = get_edgelist_from_cobra(ecoli)
yeast_edges = get_edgelist_from_cobra(yeast)

# make all of these into networkx graphs
sc_4_5_net = nx.Graph(sc_4_5_edges)
sc_5_4_net = nx.Graph(sc_5_4_edges)
sc_3_7_net = nx.Graph(sc_3_7_edges)
ecoli_net = nx.Graph(ecoli_edges)
yeast_net = nx.Graph(yeast_edges)

# make degree distribution plots for each network and arrange them side by side
(fig, axs) = plt.subplots(1,4)
sc_4_5_ax = plot_deg_dist(sc_4_5_net, axs[0], 'String Chem (4,5)')
sc_5_4_ax = plot_deg_dist(sc_5_4_net, axs[1], 'String Chem (5,4)')
sc_3_7_ax = plot_deg_dist(sc_3_7_net, axs[2], 'String Chem (3,7)')
ecoli_ax = plot_deg_dist(ecoli_net, axs[3], 'E. coli')
#yeast_ax = plot_deg_dist(yeast_net, axs[4])
plt.show()

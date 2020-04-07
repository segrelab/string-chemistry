# figure_1,py
# makes the three example string chemistry graphs shown in figure 1 of the 
# paper

import pygraphviz as gv
import string_chem_net as scn

figure = gv.AGraph(splines = 'true')

# least complex graph will have one kind of monomer and a maximum length of two
least = figure.add_subgraph()
least_net = scn.CreateNetwork('a', 2)
# start with nodes
for met in least_net.met_list:
    least.add_node(met, shape = 'box')
for rxn in least_net.rxn_list:
    least.add_node(rxn, shape = 'oval')
# then do edges
for edge in scn.make_edgelist(least_net.rxn_list):
    least.add_edge(edge)

# next we'll add one kind of monomer
mid = figure.add_subgraph()
mid_net = scn.CreateNetwork('bc', 2)
for met in mid_net.met_list:
    mid.add_node(met, shape = 'box')
for rxn in mid_net.rxn_list:
    mid.add_node(rxn, shape = 'oval')
for edge in scn.make_edgelist(mid_net.rxn_list):
    mid.add_edge(edge)

# finally we'll extend the maximum polymer length by one
most = figure.add_subgraph()
most_net = scn.CreateNetwork('de', 3)
for met in most_net.met_list:
    most.add_node(met, shape = 'box')
for rxn in most_net.rxn_list:
    most.add_node(rxn, shape = 'oval')
for edge in scn.make_edgelist(most_net.rxn_list):
    most.add_edge(edge)

# draw the figure
figure.draw('data/figure_1.png', prog = 'dot')

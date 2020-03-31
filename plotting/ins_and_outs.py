# ins_and_outs.py
# draw a relatively big network with input and output reactions to show what
# that looks like

from scripts import string_chem_net as scn
import random
import pygraphviz as gv

# make a non-trivially large network
net = scn.CreateNetwork('ab', 4)

# make a directed graph and tell it to avoid covering nodes with other nodes
# and edges when drawing it
graph = gv.AGraph(directed = True, splines = 'true')

# start with node- metabolites as squares and reactions as ovals
for met in net.met_list:
    graph.add_node(met, shape = 'box')

for rxn in net.rxn_list:
    graph.add_node(rxn, shape = 'oval')

# all edges can be the same for now
edgelist = scn.make_edgelist(net.rxn_list)
for edge in edgelist:
    graph.add_edge(edge)

# make some exchange reactions and make their edges a different color
ins = random.sample(net.met_list, 2)
outs = random.sample(net.met_list, 5)

# one input reaction per food source
for food in ins:
    graph.add_edge([f'-> {food}', food], color = 'red')

# one biomass reaction with all biomass precursors
for out in outs:
    graph.add_edge([out, ' + '.join(outs) + ' ->'], color = 'red')

# draw the graph
graph.draw('../data/ins_and_outs.png', prog = 'fdp')

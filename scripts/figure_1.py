# figure_1,py
# makes the three example string chemistry graphs shown in figure 1 of the 
# paper

import string_chem_net as scn
import pygraphviz as gv
import string

# tried making all three be subgraphs of the same graph, but the resulting
# image came out very wide and oddly low resolution, so we're making 3 separate
# images

# make the three networks we want to visualize, then loop over that list
least_net = scn.CreateNetwork('a', 2)
mid_net = scn.CreateNetwork('ab', 2)
most_net = scn.CreateNetwork('ab', 3)

nets = [least_net, mid_net, most_net]

for net in nets:
    figure = gv.AGraph(
        size = '11,8', # set size of output image
        dpi = '600', # set resolution of output image
        splines = 'true' # make sure nodes aren't overlapping with each other
    )
    # start with nodes
    for met in net.met_list:
        figure.add_node(met, shape = 'box', color = 'blue')
    for rxn in net.rxn_list:
        figure.add_node(rxn, shape = 'oval', color = 'red')
    # then do edges
    for edge in scn.make_edgelist(net.rxn_list):
        figure.add_edge(edge)

    # draw the figure and label it with the appropriate panel letter
    panel = string.ascii_lowercase[nets.index(net)]
    figure.draw(f'data/figure_1_{panel}.png', prog = 'dot')

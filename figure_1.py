# figure_1,py
# makes the three example string chemistry graphs shown in figure 1 of the 
# paper

import graphviz as gv
import string_chem_net as scn

figure = gv.Digraph('Figure 1')

# least complex graph will have one kind of monomer and a maximum length of two
with figure.subgraph() as least:
    SCN = scn.CreateNetwork('a', 2)
    for node in SCN.met_list:
        least.node(node)
    edges = scn.make_edgelist(SCN.rxn_list)
    for edge in edges:
        least.edge(*edge)
# next we'll add one kind of monomer
with figure.subgraph() as mid:
    SCN = scn.CreateNetwork('bc', 2)
    for node in SCN.met_list:
        mid.node(node)
    edges = scn.make_edgelist(SCN.rxn_list)
    for edge in edges:
        mid.edge(*edge)
# finally we'll extend the maximum polymer length by one
with figure.subgraph() as most:
    SCN = scn.CreateNetwork('de', 3)
    for node in SCN.met_list:
        most.node(node)
    edges = scn.make_edgelist(SCN.rxn_list)
    for edge in edges:
        most.edge(*edge)

figure.view()

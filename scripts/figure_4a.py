# figure_4.py
# make a figure to show the results of pruning a network

import sys
import string_chem_net as scn
import random
import pandas as pd
import pygraphviz as gv
import numpy as np
import cobra

# count number of 1s in a bitstring
def count_bitstring(bitstring):
    count = 0
    for bit in [int(bit) for bit in list(bitstring)]:
        if bit == 1:
            count += 1
    return(count)

# visualize effects of pruning
def viz_pruned_net(pruned_net, full_net, full_graph, i):
    # make a copy with differently-colored edges indicating where it was pruned
    # (use copy so full_graph isn't modified and we can do this many times)
    pruned_graph = full_graph.copy()
    # now change the colors of the reaction nodes/edges that were pruned
    for rxn in full_net.reactions:
        if rxn not in pruned_net.reactions:
            # change color of that reaction's node
            rxn_node = pruned_graph.get_node(rxn.id)
            rxn_node.attr['color'] = 'grey'
            # change color of all attached edges
            for met in rxn.metabolites:
                dropped_edge = pruned_graph.get_edge(met.id, rxn.id)
                dropped_edge.attr['color'] = 'grey'

    # now change the colors of the metabolites that are now dropped
    for met in full_net.metabolites:
        if met not in pruned_net.metabolites:
            met_node = pruned_graph.get_node(met.id)
            met_node.attr['color'] = 'grey'

    # draw the graph
    pruned_graph.draw(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_{i}.png', prog = 'fdp'
    )
    # return nothing; that file is the only necessary output

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps, export) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Arguments: monomers, max polymer length, number of food sources, ' +
        'number of biomass precursors, number of times to randomly prune, ' +
        'whether or not to allow export of all metabolites (y/n)'
    )

print('Creating universal string chemistry network.')
SCN = scn.CreateNetwork(monos, int(max_pol))
if export == 'y':
    cobra_model = scn.make_cobra_model(
        SCN.met_list, SCN.rxn_list, allow_export = True
    )
elif export == 'n':
    cobra_model = scn.make_cobra_model(
        SCN.met_list, SCN.rxn_list, allow_export = False
    )
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
scn.choose_inputs(int(ins), cobra_model, bm_rxn)
cobra_model.objective = bm_rxn

# make sure that there's at least one feasible solution before trying to prune
solution = cobra_model.optimize()
while solution.status == 'infeasible' or (solution.fluxes == 0).all():
    # remove existing biomass and input reactions
    cobra_model.remove_reactions([bm_rxn])
    in_rxns = [rxn for rxn in pruned_model.boundary if rxn.id.startswith('->')]
    cobra_model.remove_reactions(in_rxns)
    # choose new ones and see if those yield a solvable network
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    scn.choose_inputs(int(ins), cobra_model, bm_rxn)
    cobra_model.objective = bm_rxn
    solution = cobra_model.optimize()

print('Pruning network.')
min_flux_pruned = scn.min_flux_prune(cobra_model, bm_rxn)
min_flux_count = len(min_flux_pruned.reactions)
min_flux_bitstring = scn.make_bitstring(cobra_model, min_flux_pruned)

# will hold number of reactions in each network after reps runs of random_prune
random_pruned_counts = list()
# will hold bitstrings of all unique reactions and the count of times each one
# came up
random_pruned_dict = dict()
# will hold all the unique networks found by random_prune after reps runs
random_pruned_nets = list()
for i in range(1, int(reps)):
    if i % 100 == 0:
        print(f'On random prune {i}.')
    pruned_net = scn.random_prune(cobra_model, bm_rxn)
    # in order to know whether we've seen this model before, we can't just 
    # compare models, since no two models are ever 'equal', so we'll compare
    # reaction presence bitstrings. 
    # We also want to keep track of how many times we see each model, so we 
    # will make a dict with the bitstrings as keys
    # sort is in-place
    bitstring = scn.make_bitstring(cobra_model, pruned_net)
    if bitstring not in random_pruned_dict.keys():
        # make sure all reaction lists are sorted so that all isomorphic
        # networks have the same reaction list
        random_pruned_dict[bitstring] = 1
        random_pruned_nets.append(pruned_net)
    else:
        # if we already found this network once, then increment the counter
        # in random_pruned_nets by 1
        random_pruned_dict[bitstring] += 1
    # so that we can see the distribution of network sizes, record the length
    # of the reaction list each time, regardless of whether or not we've seen
    # this network before
    random_pruned_counts.append(len(pruned_net.reactions))

print('Preparing output.')
# make a dataframe of all the bitstrings we've generated
# turn the dictionary into a list of lists before coercion
bitstring_df = pd.DataFrame(list(map(list, random_pruned_dict.items())))
# add the bitstring from the min flux prune
bitstring_df = bitstring_df.append(
    pd.Series([min_flux_bitstring, 1]), ignore_index = True
)
# name columns
bitstring_df.columns = ['bitstring', 'occurrences']
# add in a column for the number of reactions in each network 
bitstring_df['rxn_count'] = list(map(count_bitstring, bitstring_df.bitstring))
# write output
bitstring_df.to_csv(
    f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_{reps}reps.csv'
)

# make a graphviz object for the initial network
full_graph = gv.AGraph(
    size = '5,5', 
    dpi = '400', 
    splines = 'true'
)

# visualize reaction fluxes with edge thickness (penwidth)
solution = cobra_model.optimize()

# distinguish metabolite and reaction nodes by shape
for met in cobra_model.metabolites:
    full_graph.add_node(met.id, shape = 'box')
for rxn in cobra_model.reactions:
    full_graph.add_node(rxn.id, shape = 'oval')
    # distinguish exchange fluxes with colored edges
    if rxn == bm_rxn or rxn in cobra_model.boundary:
        # red if they have flux or grey if they don't
        for met in rxn.metabolites:
            if solution.fluxes.loc[rxn.id] == 0:
                full_graph.add_edge(
                    [met.id, rxn.id],
                    color = 'grey'
                )
            else:
                full_graph.add_edge(
                    [met.id, rxn.id],
                    color = 'red'
                )
    else:
        for met in rxn.metabolites:
            # make reactions with no flux have grey edges
            if solution.fluxes.loc[rxn.id] == 0:
                full_graph.add_edge(
                    [met.id, rxn.id],
                    color = 'grey'
                )
            # make reactions with flux have bolded edges
            else:
                full_graph.add_edge(
                    [met.id, rxn.id],
                    penwidth = 2
                )

if export == 'y':
    full_graph.draw(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_exp_full.png',
        prog = 'fdp'
    )
elif export == 'n':
    full_graph.draw(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_noexp_full.png',
        prog = 'fdp'
    )

# visualize pruned networks
# have a counter so there can be unique filenames
i = 0
viz_pruned_net(min_flux_pruned, cobra_model, full_graph, i)
for random_pruned in random_pruned_nets:
    i += 1
    viz_pruned_net(random_pruned, cobra_model, full_graph, i)

# save the S matrix for the full network
if export == 'y':
    np.savetxt(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_exp_full_S.csv',
        cobra.util.create_stoichiometric_matrix(cobra_model), 
        delimiter = ','
    )
elif export == 'n':
    np.savetxt(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_noexp_full_S.csv',
        cobra.util.create_stoichiometric_matrix(cobra_model), 
        delimiter = ','
    )

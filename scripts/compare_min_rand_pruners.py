# compare_min_rand_pruners.py
'''
A somewhat-disorganized collection of things that prune a specified universal
string chemistry network with the minimum flux pruner and the random pruner
(pruned many times with the random pruner) and hypothetically visualize the
differences in pruned networks (not thoroughly debugged; results may vary)
'''

import sys
import string_chem_net as scn
import random
import pandas as pd
import numpy as np
import cobra

# count number of 1s in a reaction-inclusion vector
def count_included_rxns(rxn_incl):
    count = 0
    for bit in [int(rxn) for rxn in list(rxn_incl)]:
        if bit == 1:
            count += 1
    return(count)

# randomly prune a universal network several times and return all of the 
# pruned networks, their reaction inclusion vectors, and their reaction counts
def do_many_rand_prunes(full_model, bm_rxn, reps):
    # will hold number of reactions in each network 
    pruned_rxn_counts = list()
    # will hold reaction inclusion vectors of all unique networks and the 
    # number of times each one came up
    pruned_counts_dict = dict()
    # will hold all the unique networks found by random_prune after reps runs
    pruned_nets = list()
    for i in range(1, reps+1):
        if i % 100 == 0:
            print(f'On random prune {i}.')
        pruned_net = scn.random_prune(full_model, bm_rxn)
        # in order to know whether we've seen this model before, we can't just 
        # compare models, since no two models are ever 'equal', so we'll compare
        # reaction presence bitstrings. 
        # We also want to keep track of how many times we see each model, so we 
        # will make a dict with the bitstrings as keys
        # sort is in-place
        rxn_incl = scn.make_rxn_incl(full_model, pruned_net)
        if rxn_incl not in pruned_counts_dict.keys():
            # make sure all reaction lists are sorted so that all isomorphic
            # networks have the same reaction list
            pruned_counts_dict[rxn_incl] = 1
            pruned_nets.append(pruned_net)
        else:
            # if we already found this network once, then increment the counter
            # in random_pruned_nets by 1
            pruned_counts_dict[rxn_incl] += 1
        # so that we can see the distribution of network sizes, record the length
        # of the reaction list each time, regardless of whether or not we've seen
        # this network before
        pruned_rxn_counts.append(len(pruned_net.reactions))
    return(pruned_rxn_counts, pruned_counts_dict, pruned_nets)

def viz_pruned_nets(full_model, full_graph, min_pruned, rand_pruned_nets, export):
    # call min-flux network 'min' and the rest 'rand1', 'rand2' etc
    scn.viz_pruned_net(min_pruned, full_model, full_graph, export, 'min')
    i = 0
    for rand_pruned in rand_pruned_nets:
        i += 1
        scn.viz_pruned_net(rand_pruned, full_model, full_graph, export, f'rand{i}')
    # return nothing; viz_pruned_net will generate files as output

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Arguments: monomers, max polymer length, number of food sources, ' +
        'number of biomass precursors, number of times to randomly prune'
    )

print('Creating universal string chemistry networks.')
SCN = scn.CreateNetwork(monos, int(max_pol))
export_model = scn.make_cobra_model(
    SCN.met_list, SCN.rxn_list, allow_export = True
)
no_export_model = scn.make_cobra_model(
    SCN.met_list, SCN.rxn_list, allow_export = False
)
# give both networks the same biomass reaction as the objective
exp_bm_rxn = scn.choose_bm_mets(int(outs), export_model)
# have to make a copy of the reaction object or shit gets weird
no_exp_bm_rxn = exp_bm_rxn.copy()
no_export_model.add_reaction(no_exp_bm_rxn)
export_model.objective = exp_bm_rxn
no_export_model.objective = no_exp_bm_rxn
# make import reactions on no_export model first so we can just get them from
# the model's boundary for the export model
scn.choose_inputs(int(ins), no_export_model, no_exp_bm_rxn)
export_model.add_reactions([rxn.copy() for rxn in no_export_model.boundary])

# make sure that there's at least one feasible solution for both networks 
# before trying to prune either
export_solution = export_model.optimize()
no_export_solution = no_export_model.optimize()
while export_solution.status == 'infeasible' or \
    (export_solution.fluxes == 0).all() or \
    no_export_solution.status == 'infeasible' or \
    (no_export_solution.fluxes == 0).all():
    # remove biomass reaction from both networks
    export_model.remove_reactions([exp_bm_rxn])
    no_export_model.remove_reactions([no_exp_bm_rxn])
    # remove input reactions from both networks
    no_export_model.remove_reactions([rxn for rxn in no_export_model.boundary])
    export_model.remove_reactions(
        [rxn for rxn in export_model.boundary if rxn.id.startswith('->')]
    )
    # choose new input and biomass reactions
    exp_bm_rxn = scn.choose_bm_mets(int(outs), export_model)
    no_exp_bm_rxn = exp_bm_rxn.copy()
    no_export_model.add_reaction(no_exp_bm_rxn)
    export_model.objective = exp_bm_rxn
    no_export_model.objective = no_exp_bm_rxn
    scn.choose_inputs(int(ins), no_export_model, no_exp_bm_rxn)
    export_model.add_reactions(
        [rxn.copy() for rxn in no_export_model.boundary]
    )
    # see if this combination works for both networks
    export_solution = export_model.optimize()
    no_export_solution = no_export_model.optimize()

# do min-flux pruning and get reaction-inclusion vectors for both networks
print('Using minimum-flux pruner.')
min_pruned_export = scn.min_flux_prune(export_model, exp_bm_rxn)
min_pruned_no_export = scn.min_flux_prune(no_export_model, no_exp_bm_rxn)
min_export_count = len(min_pruned_export.reactions)
min_no_export_count = len(min_pruned_no_export.reactions)
min_export_rxn_incl = scn.make_rxn_incl(
    export_model, 
    min_pruned_export
)
min_no_export_rxn_incl = scn.make_rxn_incl(
    no_export_model, 
    min_pruned_no_export
)

# randomly prune each network as many times as specified
print('Randomly pruning with export reactions.')
(
    rand_export_pruned_rxn_counts,
    rand_export_prune_counts,
    rand_export_pruned_nets
) = do_many_rand_prunes(export_model, exp_bm_rxn, int(reps))

print('Randomly pruning without export reactions.')
(
    rand_no_export_pruned_rxn_counts,
    rand_no_export_prune_counts,
    rand_no_export_pruned_nets
) = do_many_rand_prunes(no_export_model, no_exp_bm_rxn, int(reps))

print('Preparing text output.')
# make dataframes of all the reaction inclusion vectors from random pruning
export_df = pd.DataFrame(
    list(map(list, rand_export_prune_counts.items()))
)
no_export_df = pd.DataFrame(
     list(map(list, rand_no_export_prune_counts.items()))
)

# add the inclusion vectors from the min flux prunes
export_df = export_df.append(
    pd.Series([min_export_rxn_incl, 1]), ignore_index = True
)
no_export_df = no_export_df.append(
    pd.Series([min_no_export_rxn_incl, 1]), ignore_index = True
)

# name existing columns
export_df.columns = ['rxn_incl', 'occurrences']
no_export_df.columns = ['rxn_incl', 'occurrences']
# add column indicating whether or not export was allowed for those networks
export_df['export'] = ['export'] * (len(rand_export_prune_counts) + 1)
no_export_df['export'] = ['no export'] * (len(rand_no_export_prune_counts) + 1)
# add column indicating which pruning algorithm was used
export_df['pruner'] = ['random'] * len(rand_export_prune_counts) + ['min flux']
no_export_df['pruner'] = ['random'] * len(rand_no_export_prune_counts) + ['min flux']
# add column of reaction counts
export_df['rxn_count'] = list(map(count_included_rxns, export_df.rxn_incl))
no_export_df['rxn_count'] = list(map(count_included_rxns, no_export_df.rxn_incl))

# concatenate dataframes and write to file
all_data = export_df.append(no_export_df)
all_data.to_csv(
    f'data/pruning_results_{monos}_{max_pol}_{ins}ins_{outs}outs_{reps}reps.csv',
    index = False
)

# write stoichiometric matrices to files
np.savetxt(
    f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_exp_full_S.csv',
    cobra.util.create_stoichiometric_matrix(export_model), 
    delimiter = ','
)
np.savetxt(
    f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_noexp_full_S.csv',
    cobra.util.create_stoichiometric_matrix(no_export_model), 
    delimiter = ','
)

# visualize networks
print('Preparing image output.')
# visualize universal networks and save the graphviz objects
export_graph = viz_universal_net(export_model, exp_bm_rxn, 'exp')
no_export_graph = viz_universal_net(no_export_model, no_exp_bm_rxn, 'noexp')

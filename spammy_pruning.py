# spammy_pruning.py
# makes one network, picks one biomass reaction, and then runs random_prune
# for several different choices of n food source metabolites

import sys
import string_chem_net as scn
import random
import re
import pandas as pd
import scipy

# count number of 1s in a bitstring
def count_bitstring(bitstring):
    count = 0
    for bit in [int(bit) for bit in list(bitstring)]:
        if bit == 1:
            count += 1
    return(count)

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps, big_reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments: monomers, max polymer length, number of food sources, \
number of biomass precursors, number of times to prune each network, number of \
times to reselect food sources.')

# create the reference network and pick some food sources and a biomass rxn
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
scn.choose_inputs(int(ins), cobra_model)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
print(f'Biomass reaction: {bm_rxn.id}')
cobra_model.objective = bm_rxn

i = 0
food_mets = list()
# use a while loop and not a for loop so we can go back on occasion
while i < int(big_reps):
    i = i + 1
    foods_string = ' '.join([met.id for met in cobra_model.boundary])
    print(f'Food source group {i}: {foods_string}')
    # choose some new food sources (remove existing ones)
    cobra_model.remove_reactions(cobra_model.boundary)
    scn.choose_inputs(int(ins), cobra_model, bm_rxn)
    # see if this choice of metabolites can produce the biomass on this network
    solution = cobra_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status == 'infeasible' or bm_rxn_flux < 1e-10:
        print('There were no feasible solutions; reselecting food sources.')
        # redo this iteration of the loop
        i = i - 1
        continue
    # record the metabolites that worked before proceeding
    else:
        food_mets.append([met.id for met in cobra_model.boundary])

    # will hold bitstrings of all unique reactions and the count of times each
    # one came up
    random_pruned_dict = dict()
    # will hold all the unique networks found by random_prune after reps runs
    random_pruned_nets = list()
    for j in range(1, int(reps)):
        if j % 100 == 0:
            print(f'Pruned {j} times')
        pruned_net = scn.random_prune(cobra_model, bm_rxn)
        # in order to know whether we've seen this model before, we can't just 
        # compare models, since no two models are ever 'equal', so we'll compare
        # reaction presence bitstrings. 
        # We also want to keep track of how many times we see each model, so we 
        # will make a dict with the bitstrings as keys
        bitstring = scn.make_bitstring(cobra_model, pruned_net)
        if bitstring not in random_pruned_dict.keys():
            # make sure all reaction lists are sorted so that all isomorphic
            # networks have the same reaction list
            random_pruned_dict[bitstring] = 1
            random_pruned_nets.append(pruned_net)
        else:
            # if we already found this network once, then increment the
            # appropriate counter by 1
            random_pruned_dict[bitstring] += 1

    # make a dataframe of all the bitstrings we've generated
    # turn the dictionary into a list of lists before coercion
    bitstring_df = pd.DataFrame(list(map(list, random_pruned_dict.items())))
    bitstring_df.columns = ['bitstring', 'occurrences']
    # add in a column for the number of reactions in each network 
    bitstring_df['rxn_count'] = list(map(
        count_bitstring, bitstring_df.bitstring
    ))
    bitstring_df.to_csv(
        f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_{i}of{big_reps}bitstrings.csv'
    )
# write the list of food sources that were tried to a file
with open(f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_foods.csv', 'w') as out:
    output = '\n'.join([','.join(sublist) for sublist in food_mets])
    out.write(output)

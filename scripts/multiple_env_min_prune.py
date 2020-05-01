# multiple_env_min_prune.py
# runs the minimum flux pruning algorithm many times on the same universal
# network with the same biomass reaction but many different environments

import sys
import string_chem_net as scn
import pandas as pd
import itertools as it

# get command-line arguments
try:
    (monos, max_pol, ins, outs, envs) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\n' +
        'number of food sources\nnumber of biomass precursors\n' +
        'number of times to reselect food sources')

# create the reference network and pick some food sources and a biomass rxn
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.choose_inputs(int(ins), cobra_model)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
cobra_model.objective = bm_rxn

i = 0
# record all of the food metabolites that were used
food_mets = list()
bitstrings = list()
# counter for how many times it had to reselct the environment to get a
# feasible solution with the full network
j = 0
# use a while loop and not a for loop so we can go back on occasion
while i < int(envs):
    i +=  1 
    # remove existing input reactions
    in_rxns = [rxn for rxn in cobra_model.boundary if rxn.id.startswith('->')]
    cobra_model.remove_reactions(in_rxns)
    # choose new input reactions
    scn.choose_inputs(int(ins), cobra_model, bm_rxn)
    in_rxns = [rxn for rxn in cobra_model.boundary if rxn.id.startswith('->')]
    foods_string = ' '.join([
        # getting the metabolite IDs out of a reaction is annoying
        list(rxn.metabolites.keys())[0].id for rxn in in_rxns
    ])
    # see if this choice of metabolites can produce the biomass on this network
    solution = cobra_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status == 'infeasible' or bm_rxn_flux < 1e-10:
        # redo this iteration of the loop
        i -= 1
        # increment the counter of redos
        j += 1 
        continue
    # record the metabolites that worked and prune the network
    else:
        if i % 100 == 0:
            print(f'On environment {i}')
        # reset the reselection counter
        j = 0
        # get the list of food source metabolites
        food_mets.append('-'.join([
            met.id 
            for rxn in in_rxns 
            for met in rxn.metabolites
        ]))
        pruned_net = scn.min_flux_prune(cobra_model, bm_rxn)
        bitstring = scn.make_bitstring(cobra_model, pruned_net)
        bitstrings.append(bitstring)

# make a dataframe out of the two lists
bitstring_df = pd.DataFrame(list(zip(food_mets, bitstrings)))
bitstring_df.columns = ['inputs','bitstring']
# add a column with the biomass components
bitstring_df['biomass'] = list(it.repeat(
    '-'.join([met.id for met in bm_rxn.metabolites]),
    len(food_mets)
))
bitstring_df.to_csv(
    f'data/multiple_env_min_prune_{monos}_{max_pol}_{ins}ins_{envs}envs_' +
    f'{outs}outs.csv', mode = 'a', index = False
)

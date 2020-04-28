# rxn_met_ratio_data.py
# makes a universal network of a specified size and then does min flux pruning
# on many combinations of biomass reaction and environments and computes the
# ratio of reactions to metabolites in each of the pruned networks

import sys
import string_chem_net as scn
import random
import cobra
import pandas as pd

# get command-line arguments
try:
    (monos, max_pol, ins, outs, env_count, bm_count) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\nnumber of food ' +
        'sources in each environment\nnumber of biomass precursors\n' +
        'number of environments\nnumber of different biomass reactions' 
    )

# create the universal network
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)

# generate all of the environments
envs = list()
i = 0
# count how many consecutive times this fails to find a new environment
fail = 0
while i < int(env_count):
    i += 1
    env = random.sample(cobra_model.metabolites, int(ins))
    # make sure this is a new environment
    if env not in envs:
        envs.append(env)
        # reset fail counter
        fail = 0
    else:
        i -= 1
        fail += 1
    if fail > 1000:
        print('Failed to find a new environment 1000 times in a row.')
        print(f'Proceeding with only {len(envs)} environments')
        break

# don't bother storing any information besides the reaction to metabolite ratio
ratios = list()

# loop over number of different biomass reactions to use
# use a counter so we can redo a trial if a biomass reaction turns out to be
# remarkably difficult to grow on
bm_trial = 0
while bm_trial < int(bm_count):
    bm_trial += 1
    print(f'On biomass reaction {bm_trial}')
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    cobra_model.objective = bm_rxn
    # loop over number of different environments to use
    # use a counter so we can redo runs if we get an infeasible solution
    env_trial = 0
    # just for fun, keep a counter of how many times pruning actually happened
    prune_count = 0
    for env in envs:
        env_trial += 1
        if env_trial % 100 == 0:
            print(f'On environment {env_trial}')
        # add appropriate input reactions and do FBA
        for met in env:
            in_rxn = cobra.Reaction(
                '->' + met.id,
                upper_bound = 100.0, # only allow importing of this metabolite
                lower_bound = 0.0
            )
            in_rxn.add_metabolites({met: 1.0})
            cobra_model.add_reaction(in_rxn)
        solution = cobra_model.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        # don't bother pruning if there is no feasible solution on the full net
        if solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
            # don't add anything to ratios
            pass
        else:
            prune_count += 1
            # do min flux pruning and get reaction to metabolite ratio
            pruned = scn.min_flux_prune(cobra_model, bm_rxn)
            ratio = len(pruned.reactions)/len(pruned.metabolites)
            ratios.append(ratio)
        # remove input reactions in preparation for next round of pruning
        in_rxns = [
            rxn for rxn in cobra_model.boundary if rxn.id.startswith('->')
        ]
        cobra_model.remove_reactions(in_rxns)
    # just for fun, print how many times pruning was successful
    print(f'This biomass reaction could sustain growth in {prune_count} ' +
    f'out of {len(envs)} environments tested.')

output = '\n'.join([str(x) for x in ratios]) + '\n'
with open(f'data/ratios_{monos}_{max_pol}_{ins}ins_{outs}outs_{env_count}envs_' + 
    f'{bm_count}bms.csv', 'w') as out:
    out.write(output)

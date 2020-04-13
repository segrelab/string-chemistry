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
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))

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
    env = 0
    # also have a counter for infeasible solutions so that we make sure we
    # don't spend forever trying to find solutions for a biomass reaction
    # that is just not particularly feasible
    infeas_count = 0
    while env < int(env_count):
        env += 1
        if env % 100 == 0:
            print(f'On environment {env}')
        # get an environment and do FBA
        scn.choose_inputs(int(ins), cobra_model, bm_rxn)
        solution = cobra_model.optimize()
        # don't bother pruning if there is no feasible solution on the full net
        if solution.status == 'infeasible' or (solution.fluxes == 0).all():
            # but also don't count this run
            env -= 1
            infeas_count += 1
        else:
            # do min flux pruning and get reaction to metabolite ratio
            pruned = scn.min_flux_prune(cobra_model, bm_rxn)
            ratio = len(pruned.reactions)/len(pruned.metabolites)
            ratios.append(ratio)
        # remove input reactions in preparation for next round of pruning
        cobra_model.remove_reactions(cobra_model.boundary)
        # if env_count environments have not produced growth with this biomass
        # reaction it's garbage get a new one
        if infeas_count > int(env_count):
            bm_trial -= 1
            break

output = '\n'.join(ratios) + '\n'
with open(f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_{env_count}envs_' + 
    '{bm_count}orgs_ratios.csv') as out:
    out.write(output)

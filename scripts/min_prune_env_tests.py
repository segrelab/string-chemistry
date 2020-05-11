# min_prune_env_tests.py
# does minimum flux pruning on a single universal scale network and then tests
# whether or not that one pruned network can grow in a bunch of randomly
# chosen environments and then repeats the whole process several times

import sys
import string_chem_net as scn
import random
import cobra

# get command-line arguments
try:
    (monos, max_pol, ins, env_count, outs, bm_count) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\nnumber of food ' +
        'sources in each environment\nnumber of environments to test\n' +
        'number of biomass precursors\nNumber of biomass reactions to try\n'
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

# loop over different biomass reactions
# list of lists to store output
growth_lists = list()
i = 0
while i < int(bm_count):
    i += 1
    print(f'On biomass reaction {i}')
    # work with a copy of the model so it remains untouched for the next
    # iteration of the loop
    universal_model = cobra_model.copy()
    # find an environment that supports growth with this environment so we can
    # prune
    bm_rxn = scn.choose_bm_mets(int(outs), universal_model)
    scn.choose_inputs(int(ins), universal_model, bm_rxn)
    universal_model.objective = bm_rxn
    solution = universal_model.optimize()
    # can't just check solution.status because sometimes it's feasible but the
    # flux through the biomass reaction is vanishingly small
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    while solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
        # if the solution isn't feasible, pick a different environment
        in_rxns = [
            # don't want to remove all boundary reactions because that would
            # also remove all of the export reactions
            rxn for rxn in pruned_model.boundary if rxn.id.startswith('->')
        ]
        pruned_model.remove_reactions(in_rxns)
        scn.choose_inputs(int(ins), universal_model, bm_rxn)
        solution = universal_model.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    # now that we know there's at least one environment that supports growth
    # with this biomass reaction, we can prune the universal network
    pruned_model = scn.min_flux_prune(universal_model, bm_rxn)
    # find growth in every environment
    for env in envs:
        # start by removing existing input reactions
        in_rxns = [
            # don't want to remove all boundary reactions because that would
            # also remove all of the export reactions
            rxn for rxn in pruned_model.boundary if rxn.id.startswith('->')
        ]
        pruned_model.remove_reactions(in_rxns)
        # while we have the pruned network with no input reactions, make the
        # reaction-inclusion vector
        bitstring = scn.make_bitstring(universal_model, pruned_model)
        # create new input reactions
        for met in env:
            in_rxn = cobra.Reaction(
                '->' + met.id,
                upper_bound = 1.0, # only allow importing of this metabolite
                lower_bound = 0.0
            )
            in_rxn.add_metabolites({met: 1.0})
            pruned_model.add_reaction(in_rxn)
        # do FBA to find growth in this environment
        solution = pruned_model.optimize()
        # prepare output
        growth = str(solution.fluxes.get(key = bm_rxn.id))
        env_string = ','.join([met.id for met in env])
        growth_lists.append([bm_rxn.id, env_string, growth, bitstring])

# write output
with open(f'data/min_env_test_{monos}_{max_pol}_{ins}ins_{outs}outs_' +
    f'{env_count}envs_{bm_count}bms.tsv', 'w') as out:
    out.write('\n'.join(['\t'.join(line) for line in growth_lists]) + '\n')

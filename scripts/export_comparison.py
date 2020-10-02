# export_comparison.py
# generates a bunch of pruned networks from one universal network where half
# were pruned with export reactions and half were pruned without

import sys
import string_chem_net as scn
import pandas as pd

def do_many_prunes(network, export, reps, ins, outs):
    '''
    Given a string_chem_net object and whether or not to allow export,
    prune the network many times and record the sizes of the pruned networks
    '''
    # create the COBRApy model
    model = scn.make_cobra_model(
        network.met_list,
        network.rxn_list,
        allow_export = export
    )
    # create a DataFrame to hold the reaction and metabolite counts of the
    # pruned networks
    output = pd.DataFrame(columns = ['rxn_count', 'met_count'])
    # prune this network reps times
    for rep in range(reps):
        if rep % 10 == 0:
            print(f'On prune {rep} of {reps}')
        # work with a copy of the model so it remains untouched for the next
        # iteration of the loop
        full_model = model.copy()
        # randomly choose the appropriate number of input and output mets
        bm_rxn = scn.choose_bm_mets(outs, full_model)
        scn.choose_inputs(ins, full_model, bm_rxn)
        full_model.objective = bm_rxn
        # see if there's a feasible solution on the full model
        solution = full_model.optimize()
        # can't just check solution.status because sometimes it's feasible but the
        # flux through the biomass reaction is vanishingly small
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        while solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
            # if the solution isn't feasible, pick a different environment
            in_rxns = [
                # don't want to remove all boundary reactions because that would
                # also remove all of the export reactions
                rxn for rxn in full_model.boundary if rxn.id.startswith('->')
            ]
            full_model.remove_reactions(in_rxns)
            scn.choose_inputs(outs, full_model, bm_rxn)
            solution = full_model.optimize()
            bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        # now that we know there's at least one environment that supports growth
        # with this biomass reaction, we can prune the universal network
        pruned_model = scn.min_flux_prune(full_model, bm_rxn)
        # count reactions and metabolites
        rxn_count = len(pruned_model.reactions)
        # metabolites aren't automatically removed when all of their reactions
        # are removed, so find out how many metabolites are left
        met_count = len([
            m for m in pruned_model.metabolites if len(m.reactions) > 0
        ])
        # add to the output dataframe
        some_output = pd.DataFrame(
            [[rxn_count, met_count]], columns = ['rxn_count', 'met_count']
        )
        output = output.append(some_output, ignore_index = True)
    return(output)

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Specify set of monomers to use, maximum string length, number of ' +
        'environmental nutrients, number of biomass precursors, and number ' +
        'of times to prune'
    )

# make the universal network
SCN = scn.CreateNetwork(monos, int(max_pol))

# prune with and without export reactions
print('Pruning with export reactions')
export_data = do_many_prunes(SCN, True, int(reps), int(ins), int(outs))
print('Pruning without export reactions')
no_export_data = do_many_prunes(SCN, False, int(reps), int(ins), int(outs))

# add columns indicating whether or not there were export reactions, merge the
# dataframes and write the result to a file
export_data['export'] = 'Allowed'
no_export_data['export'] = 'Not Allowed'
all_data = export_data.append(no_export_data, ignore_index = False)
all_data.to_csv(f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_sizes.csv')

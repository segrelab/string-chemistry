# figure_S4_S5_data.py
'''
Prune the (2,5) universal network 100 times with 2 input metabolites and 5
output metabolites, then do it again but with 6 output metabolites and keep
going until you've done every combination of 2-5 input metabolites and 5-10
output metabolites
'''

import sys
import string_chem_net as scn
import itertools as it
import multiprocessing as mp

def prune_model(universal_model, ins, outs):
    '''
    Prune the given universal model with the specified number of input and
    output metabolites
    '''
    # work with a copy of the model so it remains untouched for the next
    # iteration of the loop
    full_model = universal_model.copy()
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
        scn.choose_inputs(ins, full_model, bm_rxn)
        solution = full_model.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    # now that we know there's at least one environment that supports growth
    # with this biomass reaction, we can prune the universal network
    pruned_model = scn.min_flux_prune(full_model, bm_rxn)
    # metabolites aren't automatically removed when all of their reactions
    # are removed, so find out how many metabolites are left
    met_count = len([
        m for m in pruned_model.metabolites if len(m.reactions) > 0
    ])
    # compute the reaction-to-metabolite ratio and % of pruned reactions
    ratio = len(pruned_model.reactions)/met_count
    pruned_count = len(universal_model.reactions)-len(pruned_model.reactions)
    pruned_pct = pruned_count/len(universal_model.reactions)
    output = [ins, outs, ratio, pruned_pct]
    # make everything a string so we can join it later
    return([str(x) for x in output])

# set parameters
monos = 'ab'
max_pol = 5
min_ins = 2 # smallest number of nutrients to use
max_ins = 5 # largest number of nutrients to use
min_outs = 5 # smallest number of biomass precursors to use
max_outs = 10 # largest number of biomass precursors to use
reps = 100 # number of times to prune with each combination

# create the universal network
SCN = scn.CreateNetwork(monos, max_pol)
universal_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list, 
    allow_export = True
)

# make every pair of number of environmental metabolites and number of biomass 
# precursors
conditions = list(it.product(
    range(min_ins, max_ins + 1), # have to add 1 because Python
    range(min_outs, max_outs + 1)
))
# add the universal model to each of the sublists so that we can map this list
# of arguments to the function
args = [(universal_model,) + sublist for sublist in conditions]
# multiply the list by reps to make sure each pair of values is tested reps
# times
args = args * reps

# run this in parallel
pool = mp.Pool(mp.cpu_count())
output_data = pool.starmap(prune_model, args)

# save information as a csv
with open(f'data/figure_S4_S5_data.csv', 'w') as out:
    output = '\n'.join([','.join(row) for row in output_data])
    out.write(output)

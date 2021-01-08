# fva_pruning.py
'''
One of the reviewers suggested that we alter the pruning algorithm to account
for the fact that FBA solutions are frequently not unique, so reactions that
have zero flux in one solution may not always have zero flux, and since our
pruning algorithm immediately deletes any reaction it labels as having zero
flux, this is an important concern. So this is a variant of the pruning 
algorithm that does FVA at every step and only deletes reactions that have zero
flux at every optimal point before identifying the reaction with the lowest
flux and removing that
'''

import sys
import string_chem_net as scn
import cobra
from cobra.flux_analysis import flux_variability_analysis as fva

def fva_prune(cobra_model, bm_rxn):
    '''
    iteratively do FVA, remove all reactions with a max flux of 0, and remove
    the reaction with the smallest nonzero maximum flux until that causes the
    biomass flux to fall to 0
    '''
    # removing reactions happens in-place, so we need to make a copy of the 
    # cobra model before altering it in any way
    cobra_net = cobra_model.copy()
    # do FVA to get maximum fluxes for every reaction
    fva_data = fva(cobra_net, loopless = True)
    # save the biomass flux separately to make sure it doesn't get to zero
    bm_rxn_flux = fva_data.loc[bm_rxn.id]['maximum']
    while True:
        # remove all non-boundary reactions with a maximum flux of 0
        no_flux_rxn_ids = fva_data[fva_data['maximum'] == 0].index
        boundary_rxn_ids = [rxn.id for rxn in cobra_net.boundary]
        ids_to_remove = [
            rxn for rxn in no_flux_rxn_ids if rxn not in boundary_rxn_ids
        ]
        rxns_to_remove = [
            cobra_net.reactions.get_by_id(rxn_id) for rxn_id in ids_to_remove
        ]
        cobra_net.remove_reactions(rxns_to_remove)
        # find remaining reaction with smallest maximum flux and remove it
        flux_bearers = fva_data[fva_data['maximum'] != 0]['maximum']
        min_flux_rxn_id = flux_bearers.abs().idxmin()
        min_flux_rxn = cobra_net.reactions.get_by_id(min_flux_rxn_id)
        # if this reaction is the biomass reaction, we're clearly done pruning
        if min_flux_rxn.id == bm_rxn.id:
            break
        cobra_net.remove_reactions([min_flux_rxn])
        # see if that made the network unsolvable; if so, add the reaction back
        # and exit the while loop
        fva_data = fva(cobra_net, loopless = True)
        bm_rxn_flux = fva_data.loc[bm_rxn.id]['maximum']
        if bm_rxn_flux < 10e-10:
            cobra_net.add_reaction(min_flux_rxn)
            break
    # we kept all of the boundary reactions around until now; drop the ones
    # that have no flux
    # have to recreate fva_data first since we probably just added an essential
    # reaction back to the network after discovering that it was essential
    fva_data = fva(cobra_net, loopless = True)
    cobra_net.remove_reactions(
        fva_data[fva_data['maximum'] == 0].index
    )
    return(cobra_net)

# get command-line arguments
try:
    (monos, max_pol, ins, outs) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Specify number of monomers, max length, number of nutrients, and ' + 
        'number of biomass precursors'
    )

# make the universal model
print('Setting up string chemistry network')
SCN = scn.CreateNetwork(monos, int(max_pol))
full_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list
)

# make the nutrient uptake and biomass reactions
bm_rxn = scn.choose_bm_mets(int(outs), full_model)
full_model.objective = bm_rxn

# see if there's a feasible solution with this biomass reaction and environment
# before trying to prune anything
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
    scn.choose_inputs(int(ins), full_model, bm_rxn)
    solution = full_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)

# now prune with FVA and normally so we can compare results
print('Pruning')
fva_pruned_model = fva_prune(full_model, bm_rxn)
pruned_model = scn.min_flux_prune(full_model, bm_rxn)

# now visualize the full and pruned networks
print('Visualizing networks')
graph = scn.viz_universal_net(full_model, bm_rxn, show_all = False)
# draw the full graph
graph.layout('fdp')
graph.draw('data/fva_full.png')
# change node and edge colors to reflect results of pruning
graph = scn.viz_pruned_net(pruned_model, full_model, graph)
graph.draw('data/fva_pruning_normal_pruned.png')
# revert node and edge colors to originals
for n in graph.nodes():
    # metabolites should be blue
    if n.get_name() in SCN.met_list:
        n.attr['color'] = 'blue'
    # reactions should be red
    elif n.get_name() in SCN.rxn_list:
        n.attr['color'] = 'red'
for e in graph.edges():
    # exchange fluxes should be green and everything else should be black
    if bm_rxn.id in e:
        e.attr['color'] = 'green'
    else:
        e.attr['color'] = 'black'
# now we can grey out the reactions and metabolites not present in this pruned
# network without overlapping weirdly with the previous network
fva_pruned_graph = scn.viz_pruned_net(fva_pruned_model, full_model, graph)
fva_pruned_graph.draw('data/fva_pruning_fva_pruned.png')

# also print out all the reactions and metabolites that should be in each graph
# so we can verify that the visualization functions are working correctly
with open('data/fva_pruning_normal_pruned_rxns.txt', 'w') as out:
    out.write('\n'.join([r.id for r in pruned_model.reactions]) + '\n')
with open('data/fva_pruning_fva_pruned_rxns.txt', 'w') as out:
    out.write('\n'.join([r.id for r in fva_pruned_model.reactions]) + '\n')


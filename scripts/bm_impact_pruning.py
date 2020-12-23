# bm_impact_pruning.py
'''
One of the reviewers recommended that we try a different approach to the
pruning algorithm that weighted the decision to remove individual reactions
by the impact of that reaction's removal on the biomass flux, so this is that
pruning algorithm
'''

import sys
import string_chem_net as scn
from cobra.flux_analysis import single_reaction_deletion as get_kos

def bm_impact_prune(cobra_model, bm_rxn):
    '''
    Prune network by identifying reactions whose removal has minimal impact on
    the biomass flux and iteratively removing them until removing any more
    reactions would eliminate flux through the biomass reaction
    '''
    # removing reactions happens in-place, so we need to make a copy of the 
    # cobra model before altering it in any way
    cobra_net = cobra_model.copy()
    # assign reaction fluxes to everything before starting the loop
    solution = cobra_net.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    while True:
        # remove all non-boundary reactions with no flux
        no_flux_rxn_ids = solution.fluxes[solution.fluxes == 0].index
        boundary_rxn_ids = [rxn.id for rxn in cobra_net.boundary]
        ids_to_remove = [
            rxn for rxn in no_flux_rxn_ids if rxn not in boundary_rxn_ids
        ]
        rxns_to_remove = [
            cobra_net.reactions.get_by_id(rxn_id) for rxn_id in ids_to_remove
        ]
        cobra_net.remove_reactions(rxns_to_remove)
        # get biomass fluxes for all single reaction knockouts and identify the
        # reaction whose deletion has the smallest impact on biomass flux
        kos = get_kos(cobra_net)
        min_flux_rxn_id = list(kos['growth'].idxmax())[0]
        min_flux_rxn = cobra_net.reactions.get_by_id(min_flux_rxn_id)
        # if this reaction is the biomass reaction, we're clearly done pruning
        if min_flux_rxn.id == bm_rxn.id:
            break
        cobra_net.remove_reactions([min_flux_rxn])
        # see if that made the network unsolvable; if so, add the reaction back
        # and exit the while loop
        solution = cobra_net.optimize()
        # sometimes the solution will be feasible but the flux through the
        # biomass reaction will be some absurdly small number and then if you
        # do FBA on the same network again you won't get a feasible solution
        # so can't just check to see if the flux is 0
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
            cobra_net.add_reaction(min_flux_rxn)
            break
    # we kept all of the boundary reactions around until now; drop the ones
    # that have no flux
    # have to recreate the solution object first since we probably just added
    # an essential reaction back to the network after discovering that it was
    # essential
    solution = cobra_net.optimize()
    cobra_net.remove_reactions(solution.fluxes[solution.fluxes == 0].index)
    return(cobra_net)

def min_flux_prune(cobra_model, bm_rxn):
    '''
    Iteratively remove reactions from the network by identifying reactions with
    the smallest flux until removing a reaction causes biomass flux to drop to
    zero
    (copied from string_chem_net instead of imported because I'm adding some
    bits to keep track of how many reactions are in the network at each stage
    of pruning)
    '''
    # removing reactions happens in-place, so we need to make a copy of the 
    # cobra model before altering it in any way
    cobra_net = cobra_model.copy()
    # assign reaction fluxes to everything before starting the loop
    solution = cobra_net.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    while True:
        # remove all non-boundary reactions with no flux
        no_flux_rxn_ids = solution.fluxes[solution.fluxes == 0].index
        boundary_rxn_ids = [rxn.id for rxn in cobra_net.boundary]
        ids_to_remove = [
            rxn for rxn in no_flux_rxn_ids if rxn not in boundary_rxn_ids
        ]
        rxns_to_remove = [
            cobra_net.reactions.get_by_id(rxn_id) for rxn_id in ids_to_remove
        ]
        cobra_net.remove_reactions(rxns_to_remove)
        # find remaining reaction with smallest flux and remove it
        flux_bearers = solution.fluxes[solution.fluxes != 0]
        min_flux_rxn_id = flux_bearers.abs().idxmin()
        min_flux_rxn = cobra_net.reactions.get_by_id(min_flux_rxn_id)
        # if this reaction is the biomass reaction, we're clearly done pruning
        if min_flux_rxn.id == bm_rxn.id:
            break
        cobra_net.remove_reactions([min_flux_rxn])
        # see if that made the network unsolvable; if so, add the reaction back
        # and exit the while loop
        solution = cobra_net.optimize()
        # sometimes the solution will be feasible but the flux through the
        # biomass reaction will be some absurdly small number and then if you
        # do FBA on the same network again you won't get a feasible solution
        # so can't just check to see if the flux is 0
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
            cobra_net.add_reaction(min_flux_rxn)
            break
    # we kept all of the boundary reactions around until now; drop the ones
    # that have no flux
    # have to recreate the solution object first since we probably just added
    # an essential reaction back to the network after discovering that it was
    # essential
    solution = cobra_net.optimize()
    cobra_net.remove_reactions(solution.fluxes[solution.fluxes == 0].index)
    return(cobra_net)

try:
    (monos, max_len, ins, outs) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Specify monomers to use, maximum length of polymer, number of ' + 
        'food sources, and number of biomass precursors'
    )

# setup network of specified size
SCN = scn.CreateNetwork(monos, int(max_len))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
scn.choose_inputs(int(ins), cobra_model, bm_rxn)
cobra_model.objective = bm_rxn
print(f'Reactions in full network: {len(SCN.rxn_list)}')

# prune the network using the new and old pruning algorithms
bm_pruned = bm_impact_prune(cobra_model, bm_rxn)
min_flux_pruned = scn.min_flux_prune(cobra_model, bm_rxn)

# compare biomass fluxes between pruned networks
solution = min_flux_pruned.optimize()
print(f'Biomass flux in normally pruned network: {solution.objective_value}')
solution = bm_pruned.optimize()
print(f'Biomass flux in biomass pruned network: {solution.objective_value}')

# compare reaction lists
bm_pruned_rxn_ids = [rxn.id for rxn in bm_pruned.reactions]
min_flux_rxn_ids = [rxn.id for rxn in min_flux_pruned.reactions]

overlap = sum([1 for rxn in bm_pruned_rxn_ids if rxn in min_flux_rxn_ids])
bm_only = sum([1 for rxn in bm_pruned_rxn_ids if rxn not in min_flux_rxn_ids])
min_only = sum([1 for rxn in min_flux_rxn_ids if rxn not in bm_pruned_rxn_ids])

print(f'Number of reactions in normally pruned network: {len(min_flux_rxn_ids)}')
print(f'Number of reactions in other pruned network: {len(bm_pruned_rxn_ids)}')
print(f'Number of reactions in both: {overlap}')
print(f'Number of reactions in only normally pruned network: {min_only}')
print(f'Number of reactions in only other pruned network: {bm_only}')

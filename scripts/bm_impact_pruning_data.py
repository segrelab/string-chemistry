# bm_impact_pruning_data.py
'''
One of the reviewers recommended that we try a different approach to the
pruning algorithm that weighted the decision to remove individual reactions
by the impact of that reaction's removal on the biomass flux, so this script
contains that algorithm and some extra code to record information about how it
and the minimum flux pruner work so we can visualize the differences in the two
approaches
'''

import sys
import string_chem_net as scn
from cobra.flux_analysis import single_reaction_deletion as get_kos
import pandas as pd

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
    # record which reactions are in the network at each step of the pruning
    # process
    rxn_lists = list()
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
        # get biomass fluxes for all single reaction knockouts
        kos = get_kos(cobra_net, processes = 1)
        # make sure we don't drop a boundary reaction (takes several steps
        # because the index of kos is a bunch of frozensets of reaction ids)
        kos['rxn'] = kos.index
        kos['rxn'] = kos['rxn'].apply(lambda x: list(x)[0])
        kos = kos[~kos['rxn'].isin(boundary_rxn_ids)]
        if kos.empty:
            print('KO dataframe was emtpy after trying to remove exchange rxns')
            print(get_kos(cobra_net))
            print(boundary_rxn_ids)
            sys.exit()
        # now we can identify the reaction with the smallest impact on biomass
        # flux and drop it
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
        # now that we've made sure that reaction can be removed, add the 
        # current list of reactions to rxn_lists
        rxn_lists.append([rxn.id for rxn in cobra_net.reactions])
    # we kept all of the boundary reactions around until now; drop the ones
    # that have no flux
    # have to recreate the solution object first since we probably just added
    # an essential reaction back to the network after discovering that it was
    # essential
    solution = cobra_net.optimize()
    cobra_net.remove_reactions(solution.fluxes[solution.fluxes == 0].index)
    # and add this final reaction list to rxn_lists
    rxn_lists.append([rxn.id for rxn in cobra_net.reactions])
    return(rxn_lists)

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
    # record which reactions are in the network at each stage of pruning
    rxn_lists = list()
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
        # now that we know this network is feasible, record which reactions are
        # in it
        rxn_lists.append([rxn.id for rxn in cobra_net.reactions])
    # we kept all of the boundary reactions around until now; drop the ones
    # that have no flux
    # have to recreate the solution object first since we probably just added
    # an essential reaction back to the network after discovering that it was
    # essential
    solution = cobra_net.optimize()
    cobra_net.remove_reactions(solution.fluxes[solution.fluxes == 0].index)
    # and now all this list to rxn_lists
    rxn_lists.append([rxn.id for rxn in cobra_net.reactions])
    return(rxn_lists)

def compare_nets(universal_network, ins, outs):
    '''
    Given a network without any input reactions or a biomass reaction, add
    those and prune it using both the minimum flux pruner and the biomass
    sensitive pruner and record information about how the two pruning
    processes occurred
    '''
    # copy the universal network so we have an unmodified version around
    net = universal_network.copy()
    # choose new input and biomass metabolites for the universal network
    bm_rxn = scn.choose_bm_mets(outs, net)
    net.objective = bm_rxn
    scn.choose_inputs(ins, net, bm_rxn)
    # prune using both approaches
    bm_pruned_rxns = bm_impact_prune(net, bm_rxn)
    min_pruned_rxns = min_flux_prune(net, bm_rxn)
    # make a dictionary to store information about how this network is pruned
    # by the two algorithms
    info_dict = {
        'step' : list(),
        'type' : list(),
        'rxn_count' : list(),
        'jaccard' : list()
    }
    # counter for pruning step
    i = 0
    # zip the two lists of lists together so we can compare the lists of reactions
    for (min_step, bm_step) in zip(min_pruned_rxns, bm_pruned_rxns):
        # record how many reactions are in both networks and find the Jaccard index
        # for this pair
        min_count = len(min_step)
        bm_count = len(bm_step)
        overlap = sum([1 for rxn in min_step if rxn in bm_step])
        jaccard = overlap / (min_count + bm_count - overlap)
        # add two "rows" to info_dict- one for each model. Will make plotting this
        # info easier later
        info_dict['step'].append(i)
        info_dict['type'].append('min')
        info_dict['rxn_count'].append(min_count)
        info_dict['jaccard'].append(jaccard)
        # both step and jaccard will be the same value
        info_dict['step'].append(i)
        info_dict['type'].append('bm')
        info_dict['rxn_count'].append(bm_count)
        info_dict['jaccard'].append(jaccard)
        # now we can increment the step counter
        i += 1

    # there's no guarantee that both pruning algorithms took the same number of 
    # steps to finish, so if that's the case the zip() above will have ignored the
    # extra lists from the algorithm that took longer, so we need to look at those
    # reaction lists
    if len(bm_pruned_rxns) != len(min_pruned_rxns):
        # figure out which list was the longer one
        longer_list = list()
        longer_type = ''
        shorter_list = list()
        shorter_type = ''
        if len(bm_pruned_rxns) > len(min_pruned_rxns):
            longer_list = bm_pruned_rxns
            longer_type = 'bm'
            shorter_list = min_pruned_rxns
        else:
            longer_list = min_pruned_rxns
            longer_type = 'min'
            shorter_list = bm_pruned_rxns
        # we know that i is the index we left off at, so continue looping through
        # the longer list at i+1 and make a new counter to keep track of the number
        # of extra steps in longer_list
        j = 0
        j += i # do it this way so that i and j are actually independent variables
        for rxns in longer_list[i+1:]:
            longer_count = len(rxns)
            # get number of reactions in last list in shorter_list
            shorter_count = len(shorter_list[-1])
            overlap = sum([1 for rxn in rxns if rxn in shorter_list[-1]])
            jaccard = overlap / (longer_count + shorter_count - overlap)
            # only need to add info for longer_list
            info_dict['step'].append(j)
            info_dict['type'].append(longer_type)
            info_dict['rxn_count'].append(longer_count)
            info_dict['jaccard'].append(jaccard)
            j += 1
    # now that that's taken care of, turn this dict into a pandas dataframe and
    # return it
    info_df = pd.DataFrame(info_dict)
    return(info_df)

if __name__ == '__main__':
    try:
        (monos, max_len, ins, outs, reps) = sys.argv[1:]
    except ValueError:
        sys.exit(
            'Specify monomers to use, maximum length of polymer, number of ' + 
            'food sources, number of biomass precursors, and number of times to ' +
            'prune'
        )

    # setup network of specified size
    SCN = scn.CreateNetwork(monos, int(max_len))
    cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)

    # prepare dataframe to hold info from all rounds of pruning
    all_data = pd.DataFrame({
        'step' : list(),
        'type' : list(),
        'rxn_count' : list(),
        'jaccard' : list()
    })
    # prune reps times and store all information in all_data
    for i in range(int(reps)):
        if (i+1) % 10 == 0:
            print(f'On rep {i+1} of {reps}')
        some_data = compare_nets(cobra_model, int(ins), int(outs))
        # make another column to put i in so we can separate out all the different
        # trajectories
        some_data['trial'] = i+1
        all_data = all_data.append(some_data)
    # write all_data to a csv so we can tweak the plotting details without having
    # to run this script again
    all_data.to_csv('data/bm_impact_pruning_data.csv', index = False)

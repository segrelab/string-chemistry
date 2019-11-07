#!/usr/bin/env python3
# pruning.py
# take a network made by string_chem_net, designate some input and output
# metabolites, do fba to find reaction fluxes, drop all reactions with no flux,
# and then do one of two things:
#
# 1. remove the smallest flux, make sure that doesn't render the solution
# infeasible, and then iterating this process until you can't remove a reaction
#
# 2. remove a reaction at random, make sure that doesn't render the solution
# infeasible, then repeat until you can't remove a reaction
# 
# note that the first method will always give the same result while the second
# is liable to give a range of results in many (but not all) cases

import string_chem_net as scn
import random
import re
import matplotlib.pyplot as plt
import sys

# iteratively remove all reactions with zero flux and then the reaction with
# the smallest flux until you make the network unsolvable
def min_flux_prune(cobra_model):
    # removing reactions happens in-place, so we need to make a copy of the 
    # cobra model before altering it in any way
    cobra_net = cobra_model.copy()
    # assign reaction fluxes to everything before starting the loop
    solution = cobra_net.optimize()
    while True:
        # remove all reactions with no flux
        flux_bearers = solution.fluxes[solution.fluxes != 0]
        min_flux_rxn = flux_bearers.abs().idxmin()
        # find reaction with smallest remaining flux and remove it
        min_flux_rxn = cobra_net.reactions.get_by_id(min_flux_rxn)
        cobra_net.remove_reactions([min_flux_rxn])
        # see if that made the network unsolvable; if so, add the reaction back
        # and exit the while loop
        solution = cobra_net.optimize()
        if solution.status == 'infeasible' or (solution.fluxes == 0).all():
            cobra_net.add_reaction(min_flux_rxn)
            break
    return(cobra_net)

# iteratively choose a reaction with nonzero flux to remove until there are no
# reactions that can be removed without making the network unsolvable
# need to know what the biomass reaction is to make sure it doesn't get removed
def random_prune(cobra_model, bm_rxn):
    # removing reactions happens in-place, so we need to make a copy of the
    # cobra model before altering it in any way
    cobra_net = cobra_model.copy()
    solution = cobra_net.optimize()
    # get list of all reactions with nonzero flux
    flux_bearer_names = solution.fluxes[solution.fluxes != 0].index
    # exclude the biomass reaction, since we know we want to keep that
    flux_bearers = [
        rxn for rxn in cobra_net.reactions 
        if rxn.id in flux_bearer_names and rxn.id != bm_rxn.id
    ]
    # shuffle this list and then do a for loop over it so that we can tell if
    # we tried to remove every single possible reaction and failed (i.e. we are
    # done pruning); if we randomly chose from the list, we wouldn't ever know
    # that we actually tried every single reaction in the list, and we would
    # probably needlessly try the same reaction multiple times
    random.shuffle(flux_bearers)
    # if this ever reaches the length of flux_bearers, we are done pruning
    infeas_count = 0
    while infeas_count < len(flux_bearers):
        for rxn in flux_bearers:
            original = cobra_net.optimize()
            start = (original.fluxes < 10e-10).all()
            # try to remove the reaction from the model
            cobra_net.remove_reactions([rxn])
            # see if you can solve the model
            solution = cobra_net.optimize()
            # turns out all the fluxes sometimes get really small without 
            # actually reaching zero, but might as well be zero, so we are
            # setting a lower bound instead of actually checking for zeros
            if solution.status == 'infeasible' or (solution.fluxes < 10e-10).all():
                # keep track of how many times we've gotten a bad solution
                infeas_count += 1
                # put this reaction back and get the old solution object back
                cobra_net.add_reactions([rxn])
                solution = cobra_net.optimize()
                end = (solution.fluxes < 10e-10).all()
                if start != end:
                    print(rxn)
                    print(original.fluxes[original.fluxes != 0])
                    sys.exit()
            else:
                # recreate flux_bearers and restart the while loop and reset
                # infeas_count
                flux_bearer_names = solution.fluxes[solution.fluxes != 0].index
                flux_bearers = [
                    rxn for rxn in cobra_net.reactions
                    if rxn.id in flux_bearer_names and rxn.id != bm_rxn.id
                ]
                infeas_count = 0
                break
    return(cobra_net)

SCN = scn.CreateNetwork('ab', 5)
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
scn.choose_inputs(cobra_model, 2)
bm_rxn = scn.choose_bm_mets(cobra_model, 5)
cobra_model.objective = bm_rxn

# make sure that there's at least one feasible solution before trying to prune
solution = cobra_model.optimize()
i = 0
while solution.status == 'infeasible' or (solution.fluxes == 0).all():
    # remove existing biomass and input reactions
    cobra_model.remove_reactions([bm_rxn])
    cobra_model.remove_reactions(cobra_model.boundary)
    # choose new ones and see if those yield a solvable network
    scn.choose_inputs(cobra_model, 5)
    bm_rxn = scn.choose_bm_mets(cobra_model, 5)
    cobra_model.objective = bm_rxn
    solution = cobra_model.optimize()

min_flux_pruned = min_flux_prune(cobra_model)
min_flux_count = len(min_flux_pruned.reactions)

# will hold number of reactions in each network after 1000 runs of random_prune
random_pruned_counts = list()
# will hold unique reaction lists from randomly pruned networks because no two
# cobrapy models are ever equal, so we have to compare their lists
# will also be a dict so we can record how many times each list came up
random_pruned_net_rxns = dict()
# will hold all the unique networks found by random_prune after 1000 runs
random_pruned_nets = list()
for i in range(1,1000):
    pruned_net = random_prune(cobra_model, bm_rxn)
    # in order to know whether we've seen this model before, we can't just 
    # compare models, since no two models are ever 'equal', so we'll compare
    # reaction lists. But if we want to keep track of how many times we see
    # each model, we'll want a dict, but we can't use the models as keys and we
    # can't use the reaction lists as keys (since lists aren't hashable) so we
    # need to turn the reaction list into a string to make it a key
    # sort is in-place
    rxn_list = pruned_net.reactions
    rxn_list.sort()
    rxn_list_string = ' '.join([rxn.id for rxn in rxn_list])
    if rxn_list_string not in random_pruned_net_rxns.keys():
        # make sure all reaction lists are sorted so that all isomorphic
        # networks have the same reaction list
        random_pruned_net_rxns[rxn_list_string] = 1
        random_pruned_nets.append(pruned_net)
    else:
        # if we already found this network once, then increment the counter
        # in random_pruned_nets by 1
        random_pruned_net_rxns[rxn_list_string] += 1
    # so that we can see the distribution of network sizes, record the length
    # of the reaction list each time, regardless of whether or not we've seen
    # this network before
    random_pruned_counts.append(len(pruned_net.reactions))

# output section
print(f'Total number of possible reactions: {len(cobra_model.reactions)}')
print(f'Number of reactions left after min flux pruning: {min_flux_count}')
plt.hist(random_pruned_counts)
plt.axvline(min_flux_count)
plt.xlabel('Number of reactions')
plt.ylabel('Number of pruned networks')
plt.title('Results of 1000 Random Prunes')
plt.show()
print('Number of times each network was observed:')
for key in random_pruned_net_rxns.keys():
    rxn_count = len(re.split(' ', key))
    print(f'Network with {rxn_count} reactions: {random_pruned_net_rxns[key]}')

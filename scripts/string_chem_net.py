'''
Given a list of characters to use as monomers and a max polymer length,
create a network of all possible reactions.
'''

import itertools as it
import numpy as np
import random
import cobra
import re

class CreateNetwork():
    # given a set of monomers and a max polymer length, generate a network
    def __init__(self, monos, max_len, no_mirrors = False, make_stoich = False):
        self.met_list = self.make_met_list(monos, max_len, no_mirrors)
        self.met_set = set(self.met_list)
        self.rxn_list = self.make_rxn_list(self.met_list, no_mirrors)
        if make_stoich is True:
            self.S = self.make_stoich_mat(self.rxn_list, self.met_list)

    # given a list of metabolites with mirror-image duplicates, remove all
    # mirror-image duplicates
    # e.g. if a list contains both 'aab' and 'baa', only include one in the
    # output list
    def remove_mirrors(self, with_mirrors):
        no_mirrors = set()
        for item in with_mirrors:
            if item[::-1] not in no_mirrors:
                no_mirrors.add(item)
        return(no_mirrors)

    # given a list of monomers and a maximum polymer length, make a list of all
    # possible metabolites
    def make_met_list(self, monos, max_len, no_mirrors):
        # use itertools.product to get all possible combinations of the
        # characters in monos up to length i, then do that for every length i
        # from 1 to max_len. itertools.product returns tuples, hence the
        # coercion and joining
        met_list = list()
        if no_mirrors is True:
            met_list = [
                ''.join(list(t))
                for i in range(1, max_len+1)
                for t in self.remove_mirrors(it.product(monos, repeat = i))
            ]
        elif no_mirrors is False:
            met_list = [
                ''.join(list(t))
                for i in range(1,max_len+1)
                for t in it.product(monos, repeat = i)
            ]
        return(met_list)

    # split met into two substrings at every possible position
    def make_rxns1(self, met):
        rxn_set = set()
        if len(met) > 1:
            for i in range(1, len(met)):
                start = met[0:i]
                end = met[i:]
                rxn = f'{met}->{start}+{end}'
                # make sure the same reaction with products in the opposite
                # order isn't already in the reactions list
                rev_rxn = f'{met}->{end}+{start}'
                if rev_rxn in rxn_set:
                    pass
                else:
                    rxn_set.add(rxn)
        rxn_list = list(rxn_set)
        return(rxn_list)

    # split met into two substrings at every possible position and make sure
    # that none of the generated metabolites aren't in met_list because we
    # removed mirror-image duplicates from met_list
    def make_rxns2(self, met):
        rxn_set = set()
        if len(met) > 1:
            for i in range(1, len(met)):
                start = met[0:i]
                end = met[i:]
                # make sure that all the metabolites we generate are in met_list
                # which may require reversing them since we only kept one of the
                # two possible orderings of each metabolite
                if start not in self.met_set:
                    start = start[::-1]
                if end not in self.met_set:
                    end = end[::-1]
                rxn = f'{met}->{start}+{end}'
                # make sure the same reaction with products in the opposite
                # order isn't already in the reactions list
                rev_rxn = f'{met}->{end}+{start}'
                if rev_rxn in rxn_set:
                    pass
                else:
                    rxn_set.add(rxn)
        rxn_list = list(rxn_set)
        return(rxn_list)

    # find all decomposition reactions for all metabolites
    # only considering mono/bimolecular breakdowns because irl we know higher
    # molecularity reactions tend not to happen
    def make_rxn_list(self, met_list, no_mirrors):
        rxn_list = list()
        if no_mirrors is True:
            # get a list of lists from map so need to flatten
            rxn_list = [
                x for y in list(map(self.make_rxns2, met_list)) for x in y
            ]
        elif no_mirrors is False:
            # get a list of lists from map so need to flatten
            rxn_list = [
                x for y in list(map(self.make_rxns1, met_list)) for x in y
            ]
        return(rxn_list)

    # make a stoichiometric matrix with one column for each metabolite and one
    # row for each reactions, with the approproiate stoichiometric coefficients
    # at each position within the matrix
    def make_stoich_mat(self, rxn_list, met_list):
        # since these can get very large, we need to use a memory-mapped object
        print(f'Dimensions of matrix: {len(rxn_list)}, {len(met_list)}')
        S = np.memmap('stoich_mat.dat', dtype = np.float32, mode = 'w+',
            shape = (len(rxn_list), len(met_list)))
        # fill the array with 0s
        for i in range(0, len(met_list)):
            S[:,i] = np.zeros(len(rxn_list))
        # for each reaction, change all the appropriate zeros to actual numbers
        for i in range(0, len(rxn_list)):
            # technically we don't need to make all these variables but giving
            # things names makes this significantly more legible imo
            ref_rxn = rxn_list[i]
            bits = ref_rxn.split('->')
            reacs = bits[0].split('+')
            prods = bits[1].split('+')
            # find the indices for each of these substances in met_list so we
            # edit the appropriate column in S
            for met in reacs:
                S[i,met_list.index(met)] = -1
            for met in prods:
                # add 1 instead of assigning 1 in case two products are
                # the same, e.g. 'adad -> ad + ad'
                S[i,met_list.index(met)] += 1.0
        return(S)

# make an edgelist so this network can be easily imported into Cytoscape or
# somthing for visualization
def make_edgelist(rxn_list, rxns_as_nodes = True):
    # start by making a dictionary to look up species involved in a reaction
    # using the string notation of that reaction
    rxn_dict = dict()
    for rxn in rxn_list:
        rxn_dict[rxn] = re.split('(\+|\-\>)', rxn)[0::2]
    edgelist = []
    if rxns_as_nodes is True:
        for rxn in rxn_list:
            # could use less code but this is easier to interpret for humans
            reac = rxn_dict[rxn][0]
            prod1 = rxn_dict[rxn][1]
            prod2 = rxn_dict[rxn][2]
            edgelist.append([rxn, reac])
            edgelist.append([rxn, prod1])
            edgelist.append([rxn, prod2])
    elif rxns_as_nodes is False:
        for rxn in rxn_list:
            # could use less code but this is easier to interpret for humans
            reac = rxn_dict[rxn][0]
            prod1 = rxn_dict[rxn][1]
            prod2 = rxn_dict[rxn][2]
            edgelist.append([reac, prod1])
            edgelist.append([reac, prod2])
    return(edgelist)

# randomly choose reactions to remove given a removal probability
# removing edges is in-place, so we have to make a copy of the original graph so
# we don't modify it as well
def remove_random_rxns(more_rxns, S, prob):
    # randomly pick reactions to remove
    to_remove = []
    for rxn in more_rxns:
        if random.random() < prob:
            to_remove.append(rxn)
        else:
            pass
    # remove those reactions from the reaction list and stoichiometric matrix
    less_rxns = [x for x in more_rxns if x not in to_remove]
    # indices in more_rxns are the row indices of S
    indices = [more_rxns.index(x) for x in less_rxns]
    smaller_S = S[indices,]
    return(less_rxns, smaller_S)

# make a COBRA model using these metbolites and reactants so you can do FBA
def make_cobra_model(met_list, rxn_list):
    model = cobra.Model('string_chem')
    # start with the metabolites
    # we will need to make a dictionary with the COBRA metabolite objects
    # as keys and stoichiometric coefficients as values, so we'll need a way
    # to look up the COBRA metabolite objects using their names, since
    cobra_mets = [cobra.Metabolite(met, compartment = 'c') for met in met_list]
    met_dict = dict()
    for met in met_list:
        met_dict[met] = cobra.Metabolite(met, compartment = 'c')

    # start by making a dictionary to look up species involved in a reaction
    # using the string notation of that reaction
    rxn_dict = dict()
    for rxn in rxn_list:
        rxn_dict[rxn] = re.split('(\+|\-\>)', rxn)[0::2]

    # start by just making the COBRA reaction objects then add metabolites
    cobra_rxns = [
        # make all reactions irreversible by default
        cobra.Reaction(rxn, upper_bound = 100.0, lower_bound = 0)
        for rxn in rxn_list
    ]
    for rxn in cobra_rxns:
        # find all the COBRA metabolites associated with this reaction
        c_mets = [met_dict[met] for met in rxn_dict[rxn.id]]
        # we know the first metabolite in the list is the reactant
        rxn.add_metabolites({c_mets[0] : -1.0})
        # check if reactant splits into two identical products
        if rxn_dict[rxn.id][1] == rxn_dict[rxn.id][2]:
            rxn.add_metabolites({c_mets[1] : 2.0})
        else:
            rxn.add_metabolites({c_mets[1] : 1.0, c_mets[2] : 1.0})
        # now that the reaction is complete, we can just add it to the model
        model.add_reaction(rxn)
    return(model)

# choose n metabolites at random to create exchange reactions that are
# constrained to irreversibly produce that metabolite from nothing
def choose_inputs(n, model, bm_rxn=cobra.Reaction()):
    # make sure we don't choose any metabolites in the biomass reaction but
    # also set an empty reaction as the default biomass reaction just in case
    # we're setting food sources before biomass precursors
    in_mets = random.sample(
        [met for met in model.metabolites if met not in bm_rxn.metabolites],
        n
    )
    for met in in_mets:
        in_rxn = cobra.Reaction(
            '->' + met.id,
            upper_bound = 100.0, # only allow importing of this metabolite
            lower_bound = 0.0
        )
        in_rxn.add_metabolites({met: 1.0})
        model.add_reaction(in_rxn)
    # all of these modifications are happening in-place, so we don't actually
    # need to return anything
    return(None)

# choose n metbaolites at random to create a reaction that consummes all of them
# but produces nothing to simulate biomass production
def choose_bm_mets(n, model):
    # make sure no metabolite is both an input and biomass metabolite; start by
    # getting a list of all the metabolites that are already in boundary
    # reactions
    boundary_mets = [rxn.metabolites for rxn in model.boundary]
    bm_mets = random.sample(
        [met for met in model.metabolites if met not in boundary_mets], n
    )
    bm_rxn = cobra.Reaction(
        '+'.join([met.id for met in bm_mets]) + '->',
        upper_bound = 100.0,
        lower_bound = 0
    )
    for met in bm_mets:
        bm_rxn.add_metabolites({met:-1.0})
    model.add_reaction(bm_rxn)
    # the model object is modified in-place, so there's no need to return it
    return(bm_rxn)

def reverse_rxns(model, n):
    rxns = random.sample(model.reactions, n)
    for rxn in rxns:
        rxn.lower_bound = -100
    # all of these modifications are happening in-place, so we don't actually
    # need to return anything
    return(None)

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
            # sometimes "feasible" solutions have extremely small fluxes
            # through the biomass reaction
            bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
            if solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
                # keep track of how many times we've gotten a bad solution
                infeas_count += 1
                # put this reaction back and get the old solution object back
                cobra_net.add_reactions([rxn])
                solution = cobra_net.optimize()
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

# given a model containing all possible reactions and a pruned model only
# containing some subset of those reactions, make a bitstring indicating which
# reactions are present in the pruned model
def make_bitstring(full_model, pruned_model):
    # make sure the pruned network is actually a subnetwork of the full one
    for rxn in pruned_model.reactions:
        if rxn not in full_model.reactions:
            raise Exception(
                'Could not construct bitstrings because second network was ' +
                'not a subnetwork of the full network.\n' +
                f'Problematic reaction: {rxn.id}')
    all_reactions = full_model.reactions
    # make sure the reactions are in the same order so we can directly compare
    # multiple bitstrings from multiple different networks derived from the
    # same parent network
    all_reactions.sort()
    bits = [1 if rxn in pruned_model.reactions else 0 for rxn in all_reactions]
    bitstring = ''.join([str(bit) for bit in bits])
    return(bitstring)

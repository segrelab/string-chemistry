# comparing_nets_data.py
'''
Make a string chemistry network, prune it a bunch of times, then compare the 
degree distributions, flux distributions, and small-worldness of those networks
to the degree distributions etc. of the E. coli and yeast metabolic networks.
Writes csv files that are read by comparing_nets_plots.py to make plots
'''

import string_chem_net as scn
import cobra
import pandas as pd
import networkx as nx
from networkx.algorithms.cluster import square_clustering
from networkx.algorithms.shortest_paths.generic \
    import average_shortest_path_length
from networkx.algorithms.smallworld import random_reference, lattice_reference
import sys

def prune_once(universal_model, ins, outs, flux_bins):
    '''
    Given a universal string chemistry network, add a random biomass reaction
    and random input reactions, make sure that combination can produce biomass,
    prune the network, and return the degree and flux distributions of the 
    pruned network
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
    # get the degree and flux distributions from the pruned network
    deg_dist = make_deg_dist(pruned_model)
    fluxes = abs(pruned_model.optimize().fluxes)
    # exclude fluxes that are approximately zero
    flux_dist = pd.DataFrame(fluxes[fluxes > 10e-10])
    flux_dist.columns = ['flux']
    # get the various metrics of small-worldness
    sw_things = do_small_world_things(pruned_model)
    return((deg_dist, flux_dist, sw_things))

def make_deg_dist(model):
    '''
    Given a COBRApy model of either a real metabolic network or a pruned string
    chemistry network, find the degree and flux distributions for that network
    '''
    degs = [
        len(m.reactions) for m in model.metabolites
        # skip metabolites with degrees of 0
        if len(m.reactions) != 0
    ]
    # find frequency of each degree
    deg_dist = pd.Series(degs).value_counts(normalize = True)
    deg_dist = pd.DataFrame(deg_dist)
    deg_dist = deg_dist.reset_index()
    deg_dist.columns = ['degree', 'freq']
    return(deg_dist)

def make_scn_flux_dists(fluxes, bins):
    '''
    Given all of the fluxes in all pruned networks, find equally-sized bins
    spanning the whole range of fluxes and assign each flux to a bin
    Then get the frequency of fluxes in each bin in each pruned network
    Then find the mean and standard deviation of those bin frequencies across
    all pruned networks
    '''
    # find bins using fluxes from all trials
    fluxes['flux_bin'] = pd.cut(
        fluxes['flux'], flux_bins
    ).apply(lambda x: x.mid).astype(float)
    # get frequencies of fluxes in each bin grouped by trial
    flux_dists = pd.DataFrame(
        fluxes.groupby('trial')['flux_bin'].value_counts(normalize = True)
    )
    # rename things and bring the index back as columns
    flux_dists.columns = ['freq']
    flux_dists = flux_dists.reset_index()
    # find mean and standard deviations of frequencies within each bin
    flux_dists = flux_dists.groupby('flux_bin').agg(
        {'freq' : ['mean', 'std']}
    ).reset_index()
    flux_dists.columns = ['flux', 'mean_freq', 'std_freq']
    return(flux_dists)

def do_small_world_things(model):
    '''
    Given a COBRApy model, use networkx to find the global connectivity and 
    mean path length and compare them to those of an equal-sized lattice
    '''
    # start by getting an edgelist of the metabolic network represented by the
    # COBRApy model
    edgelist = list()
    for r in model.reactions:
        for m in r.metabolites:
            edgelist.append((m.id, r.id))
    # make a networkx object and find the average clustering coefficient and
    # the average shortest path length
    print('Computing C and L for original network')
    graph = nx.Graph(edgelist)
    ref_C = 0
    for n in graph.nodes:
        ref_C += square_clustering(graph, n)
    ref_C /= len(graph.nodes)
    ref_L = average_shortest_path_length(graph)
    # find the same parameters for an equivalent random graph
    print('Computing C and L for random network')
    rand_graph = random_reference(graph)
    rand_C = 0
    for n in rand_graph.nodes:
        rand_C += square_clustering(rand_graph, n)
    rand_C /= len(rand_graph.nodes)
    rand_L = average_shortest_path_length(rand_graph)
    # do it again but for an equivalent lattice
    print('Computing C and L for lattice')
    lat_graph = lattice_reference(graph)
    lat_C = 0
    for n in lat_graph.nodes:
        lat_C += square_clustering(lat_graph, n)
    lat_C /= len(lat_graph.nodes)
    lat_L = average_shortest_path_length(lat_graph)
    # put all these together in a pandas Series with a specified index to
    # facilitate appending this to a dataframe with these values from other
    # networks
    out = pd.Series(
        [ref_C, ref_L, rand_C, rand_L, lat_C, lat_L],
        index = ['ref_C', 'ref_L', 'rand_C', 'rand_L', 'lat_C', 'lat_L']
    )
    return(out)

def make_binned_flux_dist(model, bins):
    '''
    Given a COBRApy model, optimize production of biomass to get reaction 
    fluxes and return a Pandas DataFrame indicating what proportion of 
    fluxes fall in each of the equally-wide bins
    '''
    fluxes = abs(model.optimize().fluxes)
    # bin the fluxes into the specified number of equally-large bins after
    # dropping all of the very low fluxes that are probably supposed to be 0
    binned_fluxes = pd.cut(
        fluxes[fluxes > 10e-10], bins = bins
    )
    # this will give us the intervals but for plotting purposes we'll use the
    # midpoints of each of those intervals
    binned_fluxes = binned_fluxes.apply(lambda x: x.mid).astype(float)
    # now convert counts in each bin to frequency 
    flux_dist = pd.DataFrame(
        binned_fluxes.value_counts(normalize = True)
    )
    flux_dist = flux_dist.reset_index()
    flux_dist.columns = ['flux', 'freq']
    return(flux_dist)

# number of bins to bin flux distributions into since it's hard to visualize
# flux distributions if you don't bin the fluxes
flux_bins = 20

print('Preparing pruned string chemistry networks')
# create a string chemistry network and prune it reps times on random groups of
# nutrients and biomass precursors
monos = 'ab'
max_pol = 5
ins = 2
outs = 5
reps = 100

SCN = scn.CreateNetwork(monos, max_pol)
universal_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list, 
    allow_export = True
)
# make a dataframe to store the fluxes from each round of pruning
scn_fluxes = pd.DataFrame(columns = ['trial', 'flux'])
# make a dataframe to store the degree distributions for each pruned network
scn_deg_dists = pd.DataFrame(columns = ['trial', 'degree', 'freq'])
# make a dataframe to store the various metrics of small-worldness
scn_sw_things = pd.DataFrame(
    columns = ['ref_C', 'ref_L', 'rand_C', 'rand_L', 'lat_C', 'lat_L']
)
for rep in range(reps + 1):
    if (rep + 1) % 10 == 0:
        print(f'On rep {rep + 1} of {reps}')
    (some_degs, some_fluxes, some_sw_things) = prune_once(
        universal_model, ins, outs, flux_bins
    )
    # add a new column to each dataframe so we can distinguish the results from
    # different trials
    some_degs['trial'] = rep + 1
    scn_deg_dists = scn_deg_dists.append(some_degs)
    some_fluxes['trial'] = rep + 1
    scn_fluxes = scn_fluxes.append(some_fluxes)
    scn_sw_things = scn_sw_things.append(some_sw_things, ignore_index = True)

# get the mean and standard deviation of the frequency of each metabolite
# degree across all 100 pruned networks
scn_deg_dists = scn_deg_dists.groupby('degree').agg(
    {'freq': ['mean', 'std']}
).reset_index()
scn_deg_dists.columns = ['degree', 'mean_freq', 'std_freq']

# do the same for the flux distribution
scn_flux_dists = make_scn_flux_dists(scn_fluxes, flux_bins)

print('Preparing real metabolic networks')
# read in the real metabolic networks using COBRApy
ecoli = cobra.io.read_sbml_model('data/iJO1366.xml')
yeast = cobra.io.read_sbml_model('data/yeastGEM.xml')

# get the degree distributions of all metabolites in the real networks
ecoli_deg_dist = make_deg_dist(ecoli)
yeast_deg_dist = make_deg_dist(yeast)
# now the flux distributions
ecoli_flux_dist = make_binned_flux_dist(ecoli, flux_bins)
yeast_flux_dist = make_binned_flux_dist(yeast, flux_bins)
# now the small-world properties
ecoli_sw_things = do_small_world_things(ecoli)
yeast_sw_things = do_small_world_things(yeast)

# now write these dataframes to csv files to be read by the accompanying
# plotting script so that we can edit that script and rerun it a bunch of times
# without also needing to recreate all of this data
scn_deg_dists.to_csv('data/scn_deg_dists.csv', index = False)
scn_flux_dists.to_csv('data/scn_flux_dists.csv', index = False)
scn_sw_things.to_csv('data/scn_sw_things.csv', index = False)
ecoli_deg_dist.to_csv('data/ecoli_deg_dist.csv', index = False)
ecoli_flux_dist.to_csv('data/ecoli_flux_dist.csv', index = False)
ecolit_sw_things.to_csv('data/ecoli_sw_things.csv', index = False)
yeast_deg_dist.to_csv('data/yeast_deg_dist.csv', index = False)
yeast_flux_dist.to_csv('data/yeast_flux_dist.csv', index = False)
yeast_sw_things.to_csv('data/yeast_sw_things.csv', index = False)

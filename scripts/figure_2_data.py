# figure_2_data.py
'''
counts number of reactions and metabolites in networks with up to 5 kinds of 
monomers and polymers of up to length 10, then prunes the (3,7) string
chemistry network with 100 different randomly-selected pairs of 2-metabolite
environments and 5-metabolite biomass reactions, then compares the degree and
flux distributions of those 100 pruned networks with the degree and flux 
distributions of the E. coli and yeast metabolic networks
'''

import string
import itertools as it
import string_chem_net as scn
import cobra
import pandas as pd
import multiprocessing as mp

def prune_once(universal_model, ins, outs, flux_bins, rep):
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
    # add a column to the degree and flux distribution dataframes to indicate
    # which round of pruning this data came from
    deg_dist['trial'] = rep
    flux_dist['trial'] = rep
    return((deg_dist, flux_dist))

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
    # start by normalizing all fluxes to the maximum flux within each trial
    fluxes['flux'] = fluxes.groupby('trial')['flux'].transform(
        lambda x: x / max(x)
    )
    fluxes['flux'] /= max(fluxes['flux'])
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
    # start by normalizing all fluxes to the maximum flux within each trial
    fluxes['flux'] = fluxes.groupby('trial')['flux'].transform(
        lambda x: x / max(x)
    )
    fluxes['flux'] /= max(fluxes['flux'])
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

def make_binned_flux_dist(model, bins):
    '''
    Given a COBRApy model, optimize production of biomass to get reaction 
    fluxes and return a Pandas DataFrame indicating what proportion of 
    fluxes fall in each of the equally-wide bins
    '''
    fluxes = abs(model.optimize().fluxes)
    # normalize fluxes to the maximum flux
    fluxes /= max(fluxes)
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

# start with counting the number of reactions and metabolites in each string
# chemistry network
# list of lists to hold output; each list will have a monomer count, max 
# polymer length, metabolite count and reaction count
# there will be n x l such lists
list_o_lists = list()
for n in range(1, 6):
    for l in range(1,11):
        # get the first n letters of the alphabet
        chars = string.ascii_lowercase[:int(n)]
        # enumerate possible metabolites and reactions for this n and i
        met_count = 0
        rxn_count = 0
        for i in range(1,int(l)+1):
            new_mets = len(list(it.product(chars, repeat = i)))
            met_count += new_mets
            # tbh I forget how I derived this but I hope it's valid
            rxn_count += new_mets * (i-1)
        row_list = [n, l, met_count, rxn_count]
        list_o_lists.append([str(x) for x in row_list])

out_list = [','.join(row) for row in list_o_lists]
with open('data/figure_2_data.csv', 'w') as out:
    for row in out_list:
        out.write(row + '\n')

# now move onto the degree and flux distributions
# number of bins to bin flux distributions into since it's hard to visualize
# flux distributions if you don't bin the fluxes
flux_bins = 20

# create a string chemistry network and prune it reps times on random groups of
# nutrients and biomass precursors
monos = 'abc'
max_pol = 7
ins = 2
outs = 5
reps = 100
threads = 20

SCN = scn.CreateNetwork(monos, max_pol)
universal_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list, 
    allow_export = True
)
# do the reps rounds of pruning in parallel
pool = mp.Pool(threads)
args = [(universal_model, ins, outs, flux_bins, i+1) for i in range(reps + 1)]
mixed_data = pool.starmap(prune_once, args)
# separate the three types of data from mixed_data
scn_deg_dists = pd.concat([t[0] for t in mixed_data])
scn_fluxes = pd.concat([t[1] for t in mixed_data])

# get the mean and standard deviation of the frequency of each metabolite
# degree across all 100 pruned networks
scn_deg_dists = scn_deg_dists.groupby('degree').agg(
    {'freq': ['mean', 'std']}
).reset_index()
scn_deg_dists.columns = ['degree', 'mean_freq', 'std_freq']
# do the same for the flux distribution
scn_flux_dists = make_scn_flux_dists(scn_fluxes, flux_bins)

# read in the real metabolic networks using COBRApy
ecoli = cobra.io.read_sbml_model('data/iJO1366.xml')
yeast = cobra.io.read_sbml_model('data/yeastGEM.xml')

# get the degree distributions of all metabolites in the real networks
ecoli_deg_dist = make_deg_dist(ecoli)
yeast_deg_dist = make_deg_dist(yeast)
# now the flux distributions
ecoli_flux_dist = make_binned_flux_dist(ecoli, flux_bins)
yeast_flux_dist = make_binned_flux_dist(yeast, flux_bins)

# now write these dataframes to csv files to be read by the accompanying
# plotting script so that we can edit that script and rerun it a bunch of times
# without also needing to recreate all of this data
scn_deg_dists.to_csv('data/scn_deg_dists.csv', index = False)
scn_flux_dists.to_csv('data/scn_flux_dists.csv', index = False)
ecoli_deg_dist.to_csv('data/ecoli_deg_dist.csv', index = False)
ecoli_flux_dist.to_csv('data/ecoli_flux_dist.csv', index = False)
yeast_deg_dist.to_csv('data/yeast_deg_dist.csv', index = False)
yeast_flux_dist.to_csv('data/yeast_flux_dist.csv', index = False)

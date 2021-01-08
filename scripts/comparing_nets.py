# comparing_nets.py
'''
Make some string chemistry networks of various sizes and load some real 
metabolic networks and make some plots that compare various characteristics
'''

import string_chem_net as scn
import cobra
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

def make_flux_dist(fluxes, bins):
    '''
    Given a pandas Series containing fluxes from a metabolic network, return a
    DataFrame containing the frequency of fluxes that fall in each of the
    equally-wide bins (the number of which is the argument bins)
    '''
    binned_fluxes = pd.cut(abs(fluxes), bins = bins)
    # this will be a series of intervals and for plotting purposes we want the
    # mitpoint of each interval (also want it as a float and not an interval
    # object)
    binned_fluxes = binned_fluxes.apply(lambda x: x.mid).astype(float)
    # now get the frequency of fluxes in each of these bins
    flux_df = pd.DataFrame(
        binned_fluxes.value_counts(normalize = True)
    )
    flux_df = flux_df.reset_index()
    flux_df.columns = ['flux', 'freq']
    return(flux_df)

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
scn_flux_dists = pd.DataFrame(columns = ['trial', 'flux', 'freq'])
# make a dataframe to store the degree distributions for each pruned network
scn_deg_dists = pd.DataFrame(columns = ['trial', 'degree', 'freq'])
for rep in range(reps):
    # work with a copy of the model so it remains untouched for the next
    # iteration of the loop
    full_model = universal_model.copy()
    # randomly choose the appropriate number of input and output mets
    bm_rxn = scn.choose_bm_mets(max_pol, full_model)
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
    # get the distribution of fluxes in the pruned network
    solution = pruned_model.optimize()
    some_fluxes = make_flux_dist(solution.fluxes, flux_bins)
    # separate flux distributions from different pruned networks
    some_fluxes['trial'] = rep + 1
    scn_flux_dists = scn_flux_dists.append(some_fluxes)
    # make a dataframe with the degree distribution of this network
    # get the degree (i.e. number of associated reactions) for every metabolite
    degs = [
        len(m.reactions) for m in pruned_model.metabolites
        # skip metabolites with degrees of 0
        if len(m.reactions) != 0
    ]
    # find frequency of each degree
    deg_dist = pd.Series(degs).value_counts(normalize = True)
    some_degs = pd.DataFrame(deg_dist)
    some_degs = some_degs.reset_index()
    some_degs.columns = ['degree', 'freq']
    # separate the frequencies of each pruned network
    some_degs['trial'] = rep + 1
    scn_deg_dists = scn_deg_dists.append(some_degs)

# get the mean and standard deviation of the frequency of each metabolite
# degree across all 100 pruned networks
scn_deg_dists = scn_deg_dists.groupby('degree').agg(
    {'freq': ['mean', 'std']}
).reset_index()
scn_deg_dists.columns = ['degree', 'mean_freq', 'std_freq']

# do the same for the fluxes
scn_flux_dists = scn_flux_dists.groupby('flux').agg(
    {'freq' : ['mean', 'std']}
).reset_index()
scn_flux_dists.columns = ['flux', 'mean_freq', 'std_freq']

print('Preparing real metabolic networks')
# read in the real metabolic networks using COBRApy
ecoli = cobra.io.read_sbml_model('data/iJO1366.xml')
yeast = cobra.io.read_sbml_model('data/yeastGEM.xml')

# get the degree distributions of all metabolites in the real networks
ecoli_degs = pd.Series([
    len(m.reactions) for m in ecoli.metabolites if len(m.reactions) != 0
])
ecoli_deg_dist = pd.DataFrame(ecoli_degs.value_counts(normalize = True))
ecoli_deg_dist = ecoli_deg_dist.reset_index()
ecoli_deg_dist.columns = ['degree', 'freq']
yeast_degs = pd.Series([
    len(m.reactions) for m in yeast.metabolites if len(m.reactions) != 0
])
yeast_deg_dist = pd.DataFrame(yeast_degs.value_counts(normalize = True))
yeast_deg_dist = yeast_deg_dist.reset_index()
yeast_deg_dist.columns = ['degree', 'freq']

# get the flux distributions in the real networks
ecoli_fluxes = ecoli.optimize().fluxes
# ignore all fluxes of approximately 0
ecoli_flux_dist = make_flux_dist(
    ecoli_fluxes[ecoli_fluxes > 10e-10], flux_bins
)
# same for yeast
yeast_fluxes = yeast.optimize().fluxes
yeast_flux_dist = make_flux_dist(
    yeast_fluxes[yeast_fluxes > 10e-10], flux_bins
)

print('Plotting')
# make one plot for all of the degree distributions and three plots for the
# flux distributions in a row below
fig = plt.figure(figsize = (12,6))
gs = GridSpec(2,3, figure = fig)
deg_ax = fig.add_subplot(gs[0,1])
scn_flux_ax = fig.add_subplot(gs[1,0])
ecoli_flux_ax = fig.add_subplot(gs[1,1])
yeast_flux_ax = fig.add_subplot(gs[1,2])

# give each degree distribution a different color since they're all going on
# the same axes
# need to plot both points and error bars for the string chemistry networks
deg_ax.scatter(
    data = scn_deg_dists, 
    x = 'degree', 
    y = 'mean_freq', 
    c = 'blue',
    label = 'String Chemistries'
)
deg_ax.errorbar(
    data = scn_deg_dists,
    x = 'degree',
    y = 'mean_freq',
    yerr = 'std_freq',
    linestyle = 'None',
    c = 'blue'
)
deg_ax.scatter(
    data = ecoli_deg_dist,
    x = 'degree',
    y = 'freq',
    c = 'red',
    label = 'E. coli'
)
deg_ax.scatter(
    data = yeast_deg_dist,
    x = 'degree',
    y = 'freq',
    c = 'green',
    label = 'S. cerevisiae'
)
# set axis labels and make them log-scaled
deg_ax.set_xlabel('Degree')
deg_ax.set_ylabel('Frequency of Metabolites With Degree')
deg_ax.set_xscale('log')
deg_ax.set_yscale('log')
# make a legend so it's clear which colors represent which thing
plt.legend()
# now do the flux distributions all on different axes
scn_flux_ax.scatter(
    data = scn_flux_dists,
    x = 'flux',
    y = 'mean_freq',
    c = 'blue'
)
scn_flux_ax.errorbar(
    data = scn_flux_dists,
    x = 'flux',
    y = 'mean_freq',
    yerr = 'std_freq',
    c = 'blue',
    linestyle = 'None'
)
scn_flux_ax.set_xlabel('Flux')
scn_flux_ax.set_ylabel('Frequency of Reactions With Flux')
scn_flux_ax.set_xscale('log')
scn_flux_ax.set_yscale('log')

ecoli_flux_ax.scatter(
    data = ecoli_flux_dist,
    x = 'flux',
    y = 'freq',
    c = 'red'
)
ecoli_flux_ax.set_xlabel('Flux')
ecoli_flux_ax.set_ylabel('Frequency of Reactions With Flux')
ecoli_flux_ax.set_xscale('log')
ecoli_flux_ax.set_yscale('log')

yeast_flux_ax.scatter(
    data = yeast_flux_dist,
    x = 'flux',
    y = 'freq',
    c = 'green'
)
yeast_flux_ax.set_xlabel('Flux')
yeast_flux_ax.set_ylabel('Frequency of Reactions With Flux')
yeast_flux_ax.set_xscale('log')
yeast_flux_ax.set_yscale('log')
# make sure subplots don't overlap with each other
plt.tight_layout()
#plt.savefig('data/comparing_nets.png', dpi = 600)
plt.show()

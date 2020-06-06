# new_figure_data.py
# we'll figure out which figure this one is after it exists
# just do FBA on a bunch of networks of different sizes to see what percentage
# of reactions have flux on the full network without doing any pruning

import pandas as pd
import string
import string_chem_net as scn

# range of numbers of unique monomers to use
min_monos = 2
max_monos = 5
# range of maximum string lengths to use
min_max_len = 3
max_max_len = 10
# number of different networks of each size to generate
reps = 1000
# number of input metabolites (remains constant)
ins = 2
# number of output metabolites (remains constant)
outs = 5

# make a DataFrame to hold the results
data = pd.DataFrame(
    columns = ['monos', 'max_len', 'nonzero_count']
)

# count how many total reps there will be
total_reps = ((max_monos - min_monos) + 1) * ((max_max_len - min_max_len) + 1)

# loop over number of unique monomers first
for monos_count in range(min_monos, max_monos + 1):
    # get this number of unique characters
    monos = string.ascii_lowercase[:monos_count]
    # then loop over max string lengths
    for max_len in range(min_max_len, max_max_len + 1):
        # make a string chemistry network of this size
        SCN = scn.CreateNetwork(monos, max_len)
        cobra_model = scn.make_cobra_model(
            SCN.met_list,
            SCN.rxn_list,
            allow_export = True
        )
        # now do FBA reps times, with a different set of input and output
        # metabolites each time
        for rep in range(reps+1):
            if rep % 100 == 0:
                print(f'On rep {rep} of {total_reps}: {monos},{max_len}')
            # make a copy to add the new reactions to
            new_model = cobra_model.copy()
            bm_rxn = scn.choose_bm_mets(outs, new_model)
            scn.choose_inputs(ins, new_model, bm_rxn)
            new_model.objective = bm_rxn
            # optimize and count number of reactions with flux
            solution = new_model.optimize()
            nonzero_fluxes = len(solution.fluxes[solution.fluxes != 0])
            # add this number to the dataframe
            new_row = {
                'monos': monos_count, 
                'max_len': max_len, 
                'nonzero_count': nonzero_fluxes
            }
            data.append(new_row, ignore_index = True)

# save the dataframe so we can plot it with R
data.to_csv('data/new_figure_data.csv')

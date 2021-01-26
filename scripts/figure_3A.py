# figure_3A.py
'''
Draw a universal-scale string chemistry network and modify it to represent
which reactions have flux. Also save the corresponding stoichiometric matrix
as a text file for figure_3B.R to use
'''

import sys
import string_chem_net as scn
import random
import pandas as pd
import numpy as np
import cobra

# set parameters
monos = 'ab'
max_pol = 3
ins = 2
outs = 3

# create universal network
print('Creating string chemistry network')
SCN = scn.CreateNetwork(monos, max_pol)
model = scn.make_cobra_model(
    SCN.met_list, SCN.rxn_list, allow_export = True
)
# assign biomass and input reactions
bm_rxn = scn.choose_bm_mets(outs, model)
model.objective = bm_rxn
scn.choose_inputs(ins, model, bm_rxn)

# make sure that these input and output reactions can produce biomass
solution = model.optimize()
while solution.status == 'infeasible' or (solution.fluxes == 0).all():
    # remove biomass and input reactions
    model.remove_reactions([bm_rxn])
    model.remove_reactions(
        [rxn for rxn in model.boundary if rxn.id.startswith('->')]
    )
    # choose new input and biomass reactions
    bm_rxn = scn.choose_bm_mets(outs, model)
    model.objective = exp_bm_rxn
    scn.choose_inputs(ins, model, bm_rxn)
    # see if this combination works
    solution = model.optimize()

# write stoichiometric matrix to file for figure_3B.R
np.savetxt(
    f'data/figure_3B_data.csv',
    cobra.util.create_stoichiometric_matrix(model), 
    delimiter = ','
)

# draw network to visualize which reactions have flux
print('Drawing network')
graph = scn.viz_universal_net(model, bm_rxn, show_all = True)
graph.draw('data/figure_3A.png', prog = 'dot')

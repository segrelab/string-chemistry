# counting_plot.py
# counts number of reactions and metabolites in networks with up to 5 kinds of 
# monomers and polymers of up to length 10
# network_sizes.Rmd makes a plot using this data

import itertools as it
import sys
import string

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
with open('../data/many_counts.csv', 'w') as out:
    for row in out_list:
        out.write(row + '\n')

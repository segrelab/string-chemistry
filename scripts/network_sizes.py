# network_sizes.py
# counts the number of metabolites and reactions that would be in a string
# chemistry network with the specified number of unique monomers and maximum
# string length

import itertools as it
import sys
import string

(n,l) = sys.argv[1:]

# get the first n letters of the alphabet
chars = string.ascii_lowercase[:int(n)]

met_count = 0
rxn_count = 0
for i in range(1,int(l)+1):
    new_mets = len(list(it.product(chars, repeat = i)))
    print(f'{new_mets} metabolites of length {i}')
    met_count += new_mets
    rxn_count += new_mets * (i-1)
print(f'Total Metabolites: {met_count}')
print(f'Reactions: {rxn_count}')

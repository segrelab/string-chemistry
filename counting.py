# counting.py

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
    met_count += new_mets
    rxn_count += new_mets * (i-1)
print(f'Metabolites: {met_count}')
print(f'Reactions: {rxn_count}')

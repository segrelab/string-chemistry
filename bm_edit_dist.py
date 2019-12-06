# bm_edit_dist.py
# given a collection of groups of metabolites (for instance, various sets of 
# biomass precursors), find the average edit distance between each pair of 
# metabolites from each pair of groups of metabolites

import itertools as it
import editdistance as ed
import pandas as pd

# given two lists of strings, find the average edit distance between all pairs
# of strings
def avg_ed(l1, l2):
    dists = list()
    #print('list 1')
    #print(l1)
    #print('list 2')
    #print(l2)
    for met1 in l1:
        # if this metabolite is in the other list, give it a zero for all 
        # combinations
        if met1 in l2:
            dists.extend(list(it.repeat(0, len(l2))))
        # only compute edit distances if the metabolite isn't in the list
        else:
            for met2 in l2:
                dists.append(ed.distance(met1, met2))
    avg_dist = sum(dists)/len(dists)
    return(avg_dist)

# given a list of lists of strings, make a matrix of edit distances between all
# pairs of sublists
def make_ed_mat(met_list_list):
    print(met_list_list)
    dist_mat = list()
    # loop over indices in met_list_list
    for i in range(len(met_list_list)-1):
        # add a new empty list to dist_mat
        dist_mat.append(list())
        # for each metabolite list in met_list_list, including the ith one,
        # compute the average edit distance between it and the ith list
        for met_list in met_list_list:
            dist_mat[i].append(avg_ed(met_list_list[i], met_list))
    return(dist_mat)

# coerce to set to remove some duplicates
bms_with_dups = list(set(
    pd.read_csv('data/ab_5_100x2ins_5outs.tsv', sep = '\t').biomass
))

# might have same metabolites in different orders; need to remove those dups
# split metabolites into lists, alphabetize lists, then turn back into comma-
# delimited strings and use set coercion to remove duplicates
bms_no_dups = list(set(
    [','.join(sorted(met_list.split(','))) for met_list in bms_with_dups]
))
print(bms_no_dups)
#out = make_ed_mat(bms)
#print(out)

# random_net_freq_plot_spammer.py
# calls random_net_freq_plot.py on all specified files

import sys
import subprocess as sp

files = sys.stdin.read().split('\n')
for filename in files[:-1]:
    sp.run(['python', 'random_net_freq_plot.py', filename])

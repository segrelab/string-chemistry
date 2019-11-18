# barchart_spammer.py
# calls matplotlib_bar.py on all specified files

import sys
import subprocess as sp

files = sys.stdin.read().split('\n')
for filename in files[:-1]:
    sp.run(['python', 'matplotlib_bar.py', filename])

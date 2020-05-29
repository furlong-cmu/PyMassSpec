"""proc.py
"""

import os

from pyms.Experiment import load_expr
from pyms.DPA.PairwiseAlignment import PairwiseAlignment, align_with_tree
from pyms.DPA.Alignment import exprl2alignment
#from pyms.Peak.List.IO import store_peaks

# define the input experiments list
exprA_codes = [ "a0806_077", "a0806_078", "a0806_079" ]
exprB_codes = [ "a0806_140", "a0806_141", "a0806_142" ]

# within replicates alignment parameters
Dw = 2.5  # rt modulation [s]
Gw = 0.30 # gap penalty

# do the alignment
print('Aligning expt A')
expr_list = []
expr_dir = "../old demos/61a/output/"
for expr_code in exprA_codes:
    file_name = os.path.join(expr_dir, expr_code + ".expr")
    expr = load_expr(file_name)
    expr_list.append(expr)
F1 = exprl2alignment(expr_list)
T1 = PairwiseAlignment(F1, Dw, Gw)
A1 = align_with_tree(T1, min_peaks=2)

top_ion_list = A1.common_ion()
A1.write_common_ion_csv('output/area2.csv', top_ion_list)

print('Aligning expt B')
expr_list = []
expr_dir = "../old demos/61b/output/"
for expr_code in exprB_codes:
    file_name = os.path.join(expr_dir, expr_code + ".expr")
    expr = load_expr(file_name)
    expr_list.append(expr)
F2 = exprl2alignment(expr_list)
T2 = PairwiseAlignment(F2, Dw, Gw)
A2 = align_with_tree(T2, min_peaks=2)

# between replicates alignment parameters
Db = 10.0 # rt modulation
Gb = 0.30 # gap penalty

top_ion_list = A2.common_ion()
A2.write_common_ion_csv('output/area1.csv', top_ion_list)

print('Aligning input {1,2}')
T9 = PairwiseAlignment([A1,A2], Db, Gb)
A9 = align_with_tree(T9)

top_ion_list = A9.common_ion()
A9.write_common_ion_csv('output/area.csv', top_ion_list)

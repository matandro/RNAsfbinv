#!/usr/bin/env python3
'''
Gather good indexes from sequence / structure
Usage: gather_index.py source_sequence source_structure target_sequence target_structure
'''

import sys
from rnafbinv import shapiro_tree_aligner

if len(sys.argv) < 5:
    print("Usage: gather_index.py <source_sequence> <source_structure> <target_sequence> <target_structure>")
    exit(-1)

source_sequence = sys.argv[1].strip("'").strip('"')
source_structure = sys.argv[2].strip("'").strip('"')
target_sequence = sys.argv[3].strip("'").strip('"')
target_structure = sys.argv[4].strip("'").strip('"')
source_tree = shapiro_tree_aligner.get_tree(source_structure, source_sequence)
target_tree = shapiro_tree_aligner.get_tree(target_structure, target_sequence)
aligned_tree, score = shapiro_tree_aligner.align_trees(source_tree, target_tree)
matched, unmached = shapiro_tree_aligner.get_matching_indexes(aligned_tree)
print(matched)
print(unmached)

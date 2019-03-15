#!/usr/bin/env python3
'''
List motifs from structure for -m generation
'''

import sys
import RNAfbinvCL
import shapiro_tree_aligner
from typing import List
import shapiro_generator

def listMotif(structure:str) -> List[str]:
    if not RNAfbinvCL.is_valid_structure(structure):
        raise ValueError("ERROR: Dot bracket structure most be balanced")

    tree_root = shapiro_tree_aligner.get_tree(structure, 'N'*len(structure))
    tree_stack = [tree_root]
    motifs = []
    index = 0
    while tree_stack:
        top = tree_stack.pop()
        for child in top.children[::-1]:
            tree_stack.append(child)
        if index > 0:
            motifs.append('{}_{}{}'.format(index, top.value.name, top.value.size))
        index += 1
    return motifs


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("ERROR: usage: ListMotifs.py <dot bracket structure>")
        sys.exit(-1)
    structure = RNAfbinvCL.bracket_changer(sys.argv[1].strip().strip("'").strip('"'))
    results = listMotif(structure)
    print("Listing motifs:\nstructure {}:\nshapiro {}:"
          "\nMotif List:".format(structure, shapiro_generator.get_shapiro(structure).shapiro))
    index = 1
    for result in results:
        print("{}) {}".format(index, result))
        index += 1


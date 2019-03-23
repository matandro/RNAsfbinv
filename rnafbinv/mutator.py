#!/usr/bin/env python3
'''
enter description
'''

import random
from rnafbinv import IUPAC, tree_aligner
from typing import Dict, Any
import enum

'''
1) select motif? / select index -> identify motif? (single / double on stem)
2) identify sequence constrains ( / missing constraints)
3) select random possible mutation (addition, removal or modification)
'''


class Action(enum.Enum):
    REPLACE = 1
    ADD = 2
    REMOVE = 3


def perturbate(current_sequence: str, match_tree: tree_aligner.Tree, options: Dict[str, Any]) -> str:
    min_length = len(options.get('target_structure')) - options.get('vlength')
    max_length = len(options.get('target_structure')) + options.get('vlength')
    actions = [Action.REPLACE]
    if len(current_sequence) < max_length:
        actions.append(Action.ADD)
    if len(current_sequence) > min_length:
        actions.append(Action.REMOVE)
    # mutated_sequence = simple_point_mutation(current_sequence, random.choice(actions))
    mutated_sequence = multi_point_mutation(current_sequence, min_length, max_length, random.choice(actions))
    return mutated_sequence


def simple_point_mutation(old_sequence: str, action: Action=Action.REPLACE) -> str:
    index = random.randint(0, len(old_sequence) - 1)
    if action == Action.REPLACE:
        sequence = old_sequence[:index] + random.choice(IUPAC.IUPAC_RNA_BASE.replace(old_sequence[index], '')) + \
                   old_sequence[index + 1:]
    elif action == Action.ADD:
        sequence = old_sequence[:index] + random.choice(IUPAC.IUPAC_RNA_BASE.replace(old_sequence[index], '')) + \
                   old_sequence[index:]
    elif action == Action.REMOVE:
        sequence = old_sequence[:index] + old_sequence[index + 1:]
    return sequence


def multi_point_mutation(old_sequence: str, min_length:int, max_length: int, action: Action=Action.REPLACE,
                         max_size: int=5) -> str:
    def gen_sequence(gen_size: int, old_seq: str= None) -> str:
        res = ''
        # if we replace, select an index to be different for sure
        rand_loc = None
        if old_seq is not None:
            rand_loc = random.randint(0, len(old_seq) - 1)
        # generate new subseq
        for i in range(0, gen_size):
            selection = IUPAC.IUPAC_RNA_BASE
            if i == rand_loc:
                selection = selection.replace(old_seq[rand_loc], '')
            res += random.choice(selection)
        return res
    index = random.randint(0, len(old_sequence) - 1)
    dist = []
    for i in range(1, max_size):
        dist += [i] * (max_size - i + 1)
    if action == Action.REPLACE:
        size = min(random.choice(dist), len(old_sequence) - index)
        new_part = gen_sequence(size, old_sequence[index : index + size])
        sequence = old_sequence[:index] + new_part + old_sequence[index + size:]
    elif action == Action.ADD:
        size = min(random.choice(dist), max_length - len(old_sequence))
        new_part = gen_sequence(size)
        sequence = old_sequence[:index] + new_part + old_sequence[index:]
    elif action == Action.REMOVE:
        size = min(random.choice(dist), len(old_sequence) - index, max(1, len(old_sequence) - min_length))
        sequence = old_sequence[:index] + old_sequence[index + size:]
    return sequence

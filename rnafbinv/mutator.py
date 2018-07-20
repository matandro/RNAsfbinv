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
    mutated_sequence = simple_point_mutation(current_sequence, random.choice(actions))
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

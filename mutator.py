#!/usr/bin/env python3
'''
enter description
'''

import random
import IUPAC

'''
1) select motif? / select index -> identify motif? (single / double on stem)
2) identify sequence constrains ( / missing constraints)
3) select random possible mutation (addition, removal or modification)
'''


def perturbate(current_sequence, match_tree, options):
    mutated_sequence = simple_point_mutation(current_sequence)
    return mutated_sequence


def simple_point_mutation(sequence):
    index = random.randint(0, len(sequence) - 1)
    sequence = sequence[:index] + random.choice(IUPAC.IUPAC_RNA_BASE.replace(sequence[index], '')) + \
               sequence[index + 1:]
    return sequence

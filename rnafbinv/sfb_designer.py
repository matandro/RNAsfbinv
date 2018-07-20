#!/usr/bin/env python3
'''
Main loop for RNAsfbinv.
'''

import logging
import random
import math

from rnafbinv import shapiro_tree_aligner, vienna, tree_aligner, shapiro_generator
from rnafbinv import mutator
from rnafbinv import IUPAC

options = {}
logger = None
vienna_fold = None


def update_options(arg_map):
    global options, logger, vienna_fold
    options.update(arg_map)
    logger = options.get('logger')
    vienna_fold = options.get('RNAfold')
    logger.debug("Updating option map: {}".format(arg_map))


def stop():
    options['stop'] = True


def generate_random_start(length):
    return ''.join(random.choices(IUPAC.IUPAC_RNA_BASE, k=length))


def bp_distance(structure_a, structure_b):
    def make_pair_table(structure):
        result = [0] * len(structure)
        brackets = []
        index = 0
        for c in structure:
            if c == '(':
                brackets.append(index)
            elif c == ')':
                if len(brackets) == 0:
                    raise Exception('Error calculating bp distance: unbalanced brackets')
                base_index = brackets.pop()
                result[index] = index
                result[base_index] = index
            index += 1
        return result

    table_a = make_pair_table(structure_a)
    table_b = make_pair_table(structure_b)
    dist = 0
    for i in range(0, min(len(table_a), len(table_b))):
        if (table_a[i] != table_b[i]):
            if table_a[i] > i:
                dist += 1
            if table_b[i] > i:
                dist += 1
    return dist


def calculate_neutrality(sequence, target_structure):
    accum = 0
    seq_length = len(sequence)
    for i in range(0, seq_length):
        for c in IUPAC.IUPAC_RNA_BASE:
            if options.get('stop') is not None:
                return 0.0
            if sequence[i] != c:
                new_seq = sequence[:i] + c + sequence[i + 1:]
                structure = vienna_fold.fold(new_seq)[options.get('fold')]
                accum += bp_distance(structure, target_structure)
    return 1.0 - (accum / (pow(seq_length, 2) * 3.0))


def score_sequence(sequence, target_tree):
    # Align score tree alignment + sequence alignment
    fold_map = vienna_fold.fold(sequence)
    structure = fold_map[options.get('fold')]
    tree, score = shapiro_tree_aligner.align_trees(shapiro_tree_aligner.get_tree(structure, sequence),
                                                   target_tree)
    # Add energy diff
    target_energy = options.get('target_energy')
    if target_energy is not None and target_energy != -1000:
        score += abs(fold_map['{}_energy'.format(options.get('fold'))] - target_energy)
    # Add mutation robustness diff
    target_neutrality = options.get('target_neutrality')
    if target_neutrality is not None and target_neutrality != -1000:
        score += abs(calculate_neutrality(sequence, structure) - target_neutrality)
    return tree, score


class RnafbinvResult:
    def __init__(self, sequence: str):
        self.sequence = sequence
        fold_map = vienna_fold.fold(sequence)
        self.fold_type = options.get('fold')
        self.energy = fold_map.get("{}_energy".format(options.get('fold')))
        self.structure = fold_map.get(options.get('fold'))
        self.mutational_robustness = calculate_neutrality(self.sequence, self.structure)
        target_tree = shapiro_tree_aligner.get_tree(options['target_structure'], options['target_sequence'])
        self.result_tree = shapiro_tree_aligner.get_tree(self.structure, self.sequence)
        self.align_tree, self.score = shapiro_tree_aligner.align_trees(self.result_tree, target_tree)
        self.tree_edit_distance = tree_aligner.get_align_tree_distance(self.align_tree)
        self.bp_dist = bp_distance(self.structure, options['target_structure'])

    def __str__(self):
        print_data = "Result:\n{}\n{}\nFold energy: {}\nMutational Robustness: {}\nBP distance: {}\n" \
                     "Tree edit distance: {}\nResult tree: {}\nAligned tree ({}): {}" \
            .format(self.sequence, self.structure, self.energy, self.mutational_robustness, self.bp_dist,
                    self.tree_edit_distance, self.result_tree, self.score, self.align_tree)
        return print_data


def generate_res_object(result_seq):
    res_object = RnafbinvResult(result_seq)
    logger.info(str(res_object))
    return res_object


ALPHA = 2.0
INITIAL_TEMP = 37.0


def calc_temp(iteration, max_iteration):
    temp = INITIAL_TEMP / (1 + ALPHA * math.log(1 + iteration))
    ''' linear
    temp = (max_iteration - iteration) * INITIAL_TEMP / max_iteration
    '''
    return temp


def acceptance_probability(old_score, new_score, temperature, k):
    if new_score < old_score:
        return 1.0
    elif temperature == 0:
        return 0.0
    else:
        diff = (new_score - old_score)
        return math.exp(-diff / (k * temperature))


def merge_motifs(target_tree, motifs):
    def get_motif(motif_index):
        for motif in motifs:
            if motif.get('index') == motif_index:
                return motif
        return None

    tree_stack = [target_tree]
    index = 0
    while tree_stack:
        top = tree_stack.pop()
        for child in top.children[::-1]:
            tree_stack.append(child)
        found_motif = get_motif(index)
        if found_motif is not None:
            if top.value.size == found_motif.get('length') and \
                    top.value.name == found_motif.get('name'):
                top.value.preserve = True
            else:
                return False
        index += 1
    return True


def simulated_annealing():
    if len(options) == 0:
        logger.fatal('Options not passed to simulated annealing. must call update options before')
        return None
    # init rng
    rng_seed = options.get('rng')
    if rng_seed is not None:
        random.seed(rng_seed)
    # init loop variables
    no_iterations = options.get('iter')
    no_lookahead = options.get('look_ahead')
    # init initial sequence
    current_sequence = options.get('starting_sequence')
    if current_sequence is None and options.get('random'):
        current_sequence = generate_random_start(len(options['target_structure']))
    elif current_sequence is None:
        current_sequence = vienna.inverse(options['target_structure'],
                                          vienna.inverse_seq_ready(options['target_sequence']))
    final_result = current_sequence
    # setup target tree and get initial sequence score (and max score)
    target_tree = shapiro_tree_aligner.get_tree(options['target_structure'], options['target_sequence'])
    if not merge_motifs(target_tree, options.get('motifs')):
        shapiro_str = shapiro_generator.get_shapiro(options['target_structure']).shapiro
        logging.error('Motif list does not match target structure {}\nTarget Shapiro:{}'.format(options.get('motifs'),
                                                                                                shapiro_str))
        return None
    _, optimal_score = shapiro_tree_aligner.align_trees(target_tree, target_tree)
    match_tree, current_score = score_sequence(current_sequence, target_tree)
    best_score = current_score
    logger.info('Initial sequence ({}): {}\nAlign tree: {}'.format(current_score, current_sequence, match_tree))
    updater = options.get('updater')
    # main loop
    for iter in range(0, no_iterations):
        if options.get('stop') is not None:
            return None
        if best_score == 0:
            break
        progress = False
        for look_ahead in range(0, no_lookahead):
            if options.get('stop') is not None:
                return None
            new_sequence = mutator.perturbate(current_sequence, match_tree, options)
            new_tree, new_score = score_sequence(new_sequence, target_tree)
            temperature = calc_temp(iter, no_iterations)
            probability = acceptance_probability(current_score, new_score, temperature, len(current_sequence))
            logger.debug("iteration {} - TEMP: {} PROBABILITY: {}".format(iter + 1, temperature, probability))
            if random.random() < probability:
                progress = True
                break
            ''' OLD method, decays very fast (new is boltzman probability)
            if new_score < current_score:
                progress = True
                break
            elif random.random() < (2.0 / (iter + 1.0) / no_lookahead):
                progress = True
                break
            '''
        if progress:
            current_sequence = new_sequence
            current_score = new_score
            match_tree = new_tree
        if current_score <= best_score:
            best_score = current_score
            final_result = current_sequence
        logger.debug('Iteration {} current sequence ({}): {}\nAlign tree: {}'.format(iter + 1, current_score,
                                                                                    current_sequence, match_tree))
        if updater is not None:
            updater.update(iter + 1)
    # final print

    return final_result

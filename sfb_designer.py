#!/usr/bin/env python3
'''
Main loop for RNAsfbinv.
'''

import logging
import random

import shapiro_tree_aligner
import mutator
import vienna
import IUPAC


options = {}


def update_options(arg_map):
    options.update(arg_map)


def generate_random_start(length):
    return ''.join(random.choices(IUPAC.IUPAC_RNA_BASE, k=length))


def score_sequence(sequence, target_tree):
    structure = vienna.fold(sequence)[update_options('fold')]
    tree, score = shapiro_tree_aligner.align_trees(shapiro_tree_aligner.get_tree(structure, sequence),
                                                target_tree)
    return tree, score


def simulated_annealing():
    global options
    if len(options) == 0:
        logging.fatal('Options not passed to simulated annealing. must call update options before')
        return
    # init rng
    rng_seed = options.get('rng')
    if rng_seed is not None:
        random.seed(rng_seed)
    # init loop variables
    no_iterations = options.get('iter')
    no_lookahead = options.get('look_ahead')
    # init initial sequence
    current_sequence = options.get('starting_sequence')
    if current_sequence is None and options.get('-r') is not None:
        current_sequence = generate_random_start(len(options['target_structure']))
    final_result = current_sequence = vienna.inverse(options['target_structure'],
                                                     vienna.inverse_seq_ready(options['target_sequence']))
    # setup target tree and get initial sequence score (and max score)
    target_tree = shapiro_tree_aligner.get_tree(options['target_structure'], options['target_sequence'])
    _, optimal_score = shapiro_tree_aligner.align_trees(target_tree, target_tree)
    match_tree, current_score = score_sequence(current_sequence, target_tree)
    best_score = current_score
    # main loop
    for iter in range(0, no_iterations):
        progress = False
        for look_ahead in range(0, no_lookahead):
            new_sequence = mutator.perturbate(current_sequence, match_tree, options)
            new_tree, new_score = score_sequence(new_sequence, target_tree)
            # TODO: currently max is best score, this searches for minimum
            if new_score < current_score:
                progress = True
                break
            elif random.random() < (2.0 / (iter+1.0) / no_lookahead):
                # TODO: change transition probability as a function of the diff \
                # TODO: we have a partial test (current/new_score / optimal_score)
                progress = True
                break
        if progress:
            current_sequence = new_sequence
            current_score = new_score
            match_tree = new_tree
        # TODO: currently max is best score, this searches for minimum
        if current_score < best_score:
            best_score = current_score
            final_result = current_sequence
    # final print
    return final_result

# TODO: Global task list
# TODO: 1) Fix scoring to allow for grater diffs on motif based errors and sequence and add other minor impacts such as
# TODO:    MR and energy
# TODO: 2) Add support for settings via file (non mandatory options)
# TODO: 3) change input via file format to have headers and additional options such as starting sequence and morifs
# TODO: 4) Change mutation code to include issues motif locations and covariation mutations
# TODO: 5) Optimize (1) and (4)
# TODO: 6) Graphical interface
# TODO: 7) Update incaRNAfbinv


#!/usr/bin/env python3
'''
Main file for the RNA fragment based sequence designer.

Contains the argument analysis and log preparation
'''

__author__ = "Matan Drory Retwitzer"
__maintainer__ = "Matan Drory Retwitzer"
__email__ = "matandro@post.bgu.ac.il"

import os
import re
import sys
import logging
import shlex
from typing import List

from rnafbinv import IUPAC
from rnafbinv import vienna
from rnafbinv import sfb_designer

USAGE_TEXT = """Usage: python3 RNAsfbinvCL [Options]
    -h : Shows usage text
    -i <number of iterations> : sets the number of simulated annealing iterations (default is 100)
    --seed <random number generator seed> : a long number that is used by the random number generator
    -t <number of look ahead attempts> : number of look head mutation attempts for each iteration (default is 4)
    -e : designs a circular RNA (default is False)
    -m <motif[,...]> : comma separated list of motifs to preserve
                       motif: <motif No>[M|H|E|I|S|B]<motif No of bases>
                       Use ListMotifs.listMotif(structure) to retrieve a list of legal motifs for a given structure,  
    -s <starting sequence> : the initial sequence for the simulated annealing process
    -r : force starting simulated annealing with a random sequence
    -p <MFE|centroid> : uses RNAfold centroid or MFE folding. (default is MFE)
    --verbose : Additional info message on simulation process
    --debug : Debug information
    -l <log file path> : Logging information will be written to a given file path (rewrites file if exists)
    --length <length diff> : The resulting sequence size is target structure length +- length diff (default it 0)
    
    -f <input file path> : Path of ini file that includes mandatory information. Some options can also be set via file.
                           command line options take precedence.
    \tList of available configurations (* are mandetory and will be requested via command line if not inserted):
    *\tTARGET_STRUCTURE=<target structure>
    *\tTARGET_SEQUENCE=<target sequence>
    \tTARGET_ENERGY=<target energy>
    \tTARGET_MR=<target mutational robustness>
    \tSEED=<random seed>
    \tSTARTING_SEQUENCE=<starting sequence>
    \tITERATION=<number of simulated annealing iterations>
    """


DEF_NO_ITER = 100
DEF_NO_LOOKAHEAD = 4
MOTIF_REGEXP = re.compile(r'(?P<index>[0-9]+)(?P<name>[mMhHeEiIsSbB])(?P<length>[0-9]+)$')


def check_motif(motifs):
    res = []
    for motif in motifs:
        clean = motif.strip()
        if clean != '':
            match = MOTIF_REGEXP.match(clean)
            if match is None:
                return None
            res.append({'index': int(match.group('index')), 'name': match.group('name').upper(),
                        'length': int(match.group('length'))})
    return res


# Handles argument from command line
def generate_arg_map(argv):
    #def check_mandetory():
    #    error = None
    def index(key):
        try:
            i = argv.index(key)
        except ValueError:
            i = None
        return i

    arg_map = {}
    logger = logging.getLogger('RNAsfbinv')
    formatter = logging.Formatter('%(levelname)s:%(asctime)s - %(message)s')
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    arg_map['logger'] = logger
    # -h
    item_index = index('-h')
    if item_index is not None:
        arg_map['error'] = 'Help'
    # -l <log file path>
    item_index = index('-l')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-l requires path for log file'
            return arg_map
        log_path = argv[item_index + 1]
        if os.path.exists(log_path):
            logging.warning('Log file already exists, appending to it')
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # --verbose or --debug (logging level info instead of warning)
    if index('--debug') is not None:
        logger.setLevel(logging.DEBUG)
    elif index('--verbose') is not None:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    # -p <MFE|centroid>
    item_index = index('-p')
    if item_index is not None:
        if item_index + 1 >= len(argv) or argv[item_index + 1] not in ['MFE', 'centroid']:
            arg_map['error'] = '-p requires either MFE or centroid'
            return arg_map
        arg_map['fold'] = argv[item_index + 1]
    else:
        arg_map['fold'] = 'MFE'
    # -i <number of iterations>
    item_index = index('-i')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-i requires number of iterations'
            return arg_map
        try:
            iterations = int(argv[item_index + 1])
        except ValueError:
            arg_map['error'] = '-i number of iteration should be an integer'
            return arg_map
        arg_map['iter'] = iterations
    else:
        arg_map['iter'] = DEF_NO_ITER
    # --seed <RNG seed, long>
    item_index = index('--seed')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '--seed requires RNG seed'
            return arg_map
        try:
            seed = int(argv[item_index + 1])
        except ValueError:
            arg_map['error'] = '--seed RNG seed should be a long integer'
            return arg_map
        arg_map['rng'] = seed
    # -t <number of look ahead attempts>
    item_index = index('-t')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-t requires number of look ahead attempts'
            return arg_map
        try:
            look_ahead = int(argv[item_index + 1])
        except ValueError:
            arg_map['error'] = '-t number of look ahead attempts should be an integer'
            return arg_map
        arg_map['look_ahead'] = look_ahead
    else:
        arg_map['look_ahead'] = DEF_NO_LOOKAHEAD
    # -e sets the RNA as circular
    item_index = index('-e')
    if item_index is not None:
        arg_map['circular'] = True
    else:
        arg_map['circular'] = False
    # -m <motif to preserve[,...]>
    item_index = index('-m')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-m requires comma separated list of motif to preserve'
            return arg_map
        motifs = argv[item_index + 1].split(',')
        motifs = check_motif(motifs)
        if motifs is None:
            arg_map['error'] = 'motif format should be: <motif No>[M|H|E|I|S|B]<motif No of bases>,...'
            return arg_map
        arg_map['motifs'] = motifs
    else:
        arg_map['motifs'] = []
    # -s <starting sequence>
    item_index = index('-s')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-s requires starting sequence'
            return arg_map
        # TODO: check if sequence is legal. maybe it should only include AGCU/T. check same length?
        arg_map['starting_sequence'] = argv[item_index + 1].upper().replace('T', 'U')
    # -r sets random sequence start
    item_index = index('-r')
    if item_index is not None:
        if arg_map.get('starting_sequence') is not None:
            arg_map['error'] = '-s <starting sequence> and -r do not work together'
            return arg_map
        arg_map['random'] = True
    else:
        arg_map['random'] = False
    # --length <length>
    item_index = index('--length')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '--length require maximum length offset'
        try:
            length = int(argv[item_index + 1])
        except ValueError:
            arg_map['error'] = '-t maximum length offset should be an integer'
            return arg_map
        arg_map['vlength'] = length
    else:
        arg_map['vlength'] = 0
    # -o <output log file> TODO: replace
    #item_index = index('-o')
    #if item_index is not None:
    #    if item_index + 1 >= len(argv):
    #        arg_map['error'] = '-o requires path to output log file'
    #        return arg_map
    #    error = setup_log_file(argv[item_index + 1])
    #    if error is not None:
    #        arg_map['error'] = error
    #        return arg_map
    # -f <input file path>- Keep last

    item_index = index('-f')
    if item_index is not None:
        if item_index + 1 >= len(argv):
            arg_map['error'] = '-f requires path to input file'
            return arg_map
        temp_map = arg_map
        arg_map = read_input_file(argv[item_index + 1])
        arg_map.update(temp_map)
    if arg_map.get('error') is None:
        arg_map = read_mandatory_params(arg_map, item_index is not None)
        return arg_map
    return arg_map


# check for valid basic (only '( 'and ')' for brackets) dot bracket structure
def is_valid_structure(structure):
    bracket_count = 0
    for c in structure:
        if c == '(':
            bracket_count += 1
        elif c == ')':
            bracket_count -= 1
            if bracket_count < 0:
                return False
        elif c != '.':
            return False
    return bracket_count == 0


# changes any type of bracket that isn't round brackets into one
def bracket_changer(change_structure):
    changed_structure = ''
    for c in change_structure:
        if c in '{[<(':
            changed_structure += '('
        elif c in '}]>)':
            changed_structure += ')'
        elif c == '.':
            changed_structure += c
        elif not c.isspace():
            raise ValueError('Only except brackets and . (ignore whitespaces) [{}]'.format(c))
    return changed_structure


# reads an input file with mandatory parameters. the file should contain 4 lines:
# dot bracket structure
# IUPAC sequence
# desired folding energy (in kcal/mol). -1000 means no target
# desired neutrality. -1000 means no target
def read_input_file(input_file_path):
    input_map = {}

    def generate_file_map():
        res_map = {}
        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                clean_line = line.strip()
                # ignore comments and empty lines
                if clean_line[0] != '#' and clean_line != ';' and clean_line != '':
                    if '=' not in clean_line:
                        input_map['error'] = 'illegal input file format. Supports commends (starting with #) ' \
                                             'or Key=Value lines only'
                        return None
                    key, value = clean_line.split('=', 1)
                    res_map[key.strip().upper()] = value.strip()
        return res_map

    if not os.path.isfile(input_file_path):
        input_map['error'] = "file {} does not exists".format(input_file_path)
    else:
        config_map = generate_file_map()
        if config_map is not None:
            file_errors = []
            structure = config_map.get('TARGET_STRUCTURE')
            if structure is not None:
                structure = bracket_changer(structure)
                if not is_valid_structure(structure):
                    file_errors.append("target structure can hold '.' '(\\{\\[\\<' and ')\\}\\]\\>'" \
                                   " only and must be balanced")
                input_map['target_structure'] = structure
            sequence = config_map.get('TARGET_SEQUENCE')
            if sequence is not None:
                if not IUPAC.is_valid_sequence(sequence):
                    file_errors.append('target sequence must hold legal IUPAC letters')
                if not len(structure) == len(sequence):
                    file_errors.append('sequence and structure must be at the same length')
                input_map['target_sequence'] = sequence.upper().replace('T', 'U')
            target_energy = config_map.get('TARGET_ENERGY')
            if target_energy is not None:
                try:
                    target_energy = float(target_energy)
                    input_map['target_energy'] = target_energy
                except ValueError:
                    file_errors.append('target energy must be a floating number')
            target_mr = config_map.get('TARGET_MR')
            if target_mr is not None:
                try:
                    target_mr = float(target_mr)
                    input_map['target_neutrality'] = target_mr
                except ValueError:
                    file_errors.append('target neutrality must be a floating number')
            seed = config_map.get('SEED')
            if seed is not None:
                try:
                    seed = int(seed)
                    input_map['rng'] = seed
                except ValueError:
                    file_errors.append('RNG seed should be a long integer')
            starting_sequence = config_map.get('STARTING_SEQUENCE')
            if starting_sequence is not None:
                input_map['starting_sequence'] = starting_sequence.upper().replace('T', 'U')
            iteration_no = config_map.get('ITERATION')
            if iteration_no is not None:
                try:
                    iteration_no = int(iteration_no)
                    input_map['iter'] = iteration_no
                except ValueError:
                    file_errors.append('number of iteration should be an integer')

            if len(file_errors) > 0:
                input_map['error'] = 'Input file must have only commend, empty and Key=Value lines.' \
                                     'the following errors were found:\n{}'.format('\n'.join(file_errors))
    return input_map


# prints the usage. if error is not Help, prints an error
def usage(error='Help'):
    if error != 'Help':
        logging.error(error)
    logging.info(USAGE_TEXT)
    return


# setups output file for results only TODO: maybe not via logger
def setup_file(out_path, formatter=logging.Formatter('%(message)s')):
    error = None
    if os.path.exists(out_path):
        logging.warning('Log file already exists, appending to it')
    # TODO: check what error comes out if there are no permission to write to file
    handler = logging.FileHandler(out_path)
    handler.setFormatter(formatter)
    logger = logging.getLogger('file')
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return error


# reads mandatory parameters from command line
def read_mandatory_params(input_map, was_file):
    while input_map.get('target_structure') is None:
        structure = bracket_changer(input('Enter structure in dot bracket notation:\n').strip())
        if is_valid_structure(structure):
            input_map['target_structure'] = structure
        print('structure must be balanced and in dot bracket notation')
    while input_map.get('target_sequence') is None:
        sequence = input('Enter IUPAC format sequence restrictions (must be same length as fold):\n')
        if IUPAC.is_valid_sequence(sequence) and len(input_map['target_structure']) == len(sequence):
            input_map['target_sequence'] = sequence.upper().replace('T', 'U')
        print('sequence must be in IUPAC format and be the same length as structure')
    while input_map.get('target_energy') is None and not was_file:
        str_energy = input('Enter the desired minimum free energy in Kcal/mol (-1000 for none):\n')\
            .strip()
        try:
            input_map['target_energy'] = float(str_energy)
        except ValueError:
            print('target energy must by a floating number')
    while input_map.get('target_neutrality') is None and not was_file:
        str_neutrality = input('Enter the desired neutrality (-1000 for none):\n')\
            .strip()
        try:
            input_map['target_neutrality'] = float(str_neutrality)
            break
        except ValueError:
            print('target neutrality must by a floating number')
    return input_map


def designer(argv: List[str]) -> sfb_designer.RnafbinvResult:
    result = None
    logging.info("Starting RNAsfbinv, arguments: {}".format(argv))
    arg_map = generate_arg_map(argv)
    logging.info("Run argument map:\n{}".format(arg_map))
    error = arg_map.get('error')
    if error is not None:
        usage(error)
    else:
        arg_map.get("logger").debug("Argument map:\n{}".format(arg_map))
        # init RNAfold
        rna_folder = vienna.LiveRNAfold(arg_map.get("logger"))
        rna_folder.start(arg_map.get('circular'))
        arg_map['RNAfold'] = rna_folder
        # run simulated annealing
        arg_map.get("logger").debug("Starting simulated_annealing\nArguments: {}".format(arg_map))
        designed_sequence = sfb_designer.simulated_annealing(arg_map)
        arg_map.get("logger").debug("Finished simulated_annealing\nSequence: {}".format(designed_sequence))
        if designed_sequence is not None:
            logging.info("Finished simulated annealing, resulting sequence: {}".format(designed_sequence))
            result = sfb_designer.generate_res_object(designed_sequence, arg_map)
            print(str(result))
        else:
            logging.error("Failed to design, Exisiting!")
        rna_folder.close()
    return result


def main(command_line_args: str) -> sfb_designer.RnafbinvResult:
    return designer(shlex.split(command_line_args))


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING,
                        format='%(levelname)s:%(asctime)s - %(message)s')
    vienna.set_vienna_path('D:\\Programs\\ViennaRNA')
    designer(sys.argv[1:])

#!/usr/bin/env python3
'''
Main file for the RNA fragment based sequence designer.

Contains the argument analysis and log preparation
'''

__author__ = "Matan Drory Retwitzer"
__maintainer__ = "Matan Drory Retwitzer"
__email__ = "matandro@post.bgu.ac.il"

import os
import sys
import logging

import IUPAC
import vienna
import sfb_designer

USAGE_TEXT = """Usage: python3 RNAsfbinv [Options]
    -h : Shows usage text
    -i <number of iterations> : sets the number of simulated annealing iterations (default is 100)
    --seed <random number generator seed> : a long number that is used by the random number generator
    -t <number of look ahead attempts> : number of look head mutation attempts for each iteration (default is 4)
    -e : designs a circular RNA (default is False)
    -m <motif[,...]> : comma separated list of motifs to preserve,  
    -s <starting sequence> : the initial sequence for the simulated annealing process
    -r : force starting simulated annealing with a random sequence
    -f <input file path> : path of file that includes mandatory information. if not given, the information will be requested by command line.
    -p <MFE|centroid> : uses RNAfold centroid or MFE folding. (default is MFE)
    --verbose : Additional info message on simulation process
    --debug : Debug information
    -l <log file path> : Logging information will be written to a given file path (rewrites file if exists)
    --length <length diff> : The resulting sequence size is target structure length +- length diff (default it 0)"""  # TODO: add example for -m


DEF_NO_ITER = 100
DEF_NO_LOOKAHEAD = 4


# Handles argument from command line
def generate_arg_map(argv):
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
        logger.setLevel(logging.INFO)
    elif index('--verbose') is not None:
        logger.setLevel(logging.DEBUG)
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
            arg_map['error'] = '-p requires RNG seed'
            return arg_map
        try:
            seed = int(argv[item_index + 1])
        except ValueError:
            arg_map['error'] = '-p RNA seed should be a long integer'
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
        # TODO: check format for motif list
        motifs = argv[item_index + 1].split(',')
        arg_map['motifs'] = motifs
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
        arg_map.update(read_input_file(argv[item_index + 1]))
    else:
        arg_map.update(read_mandatory_params())
    if arg_map.get('error') is not None:
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
    return True


# changes any type of bracket that isn't round brackets into one
def bracket_changer(change_structure):
    changed_structure = ''
    for c in change_structure:
        if c in '{[<':
            changed_structure += '('
        elif c in '}]>':
            changed_structure += ')'
        else:
            changed_structure += c
    return changed_structure


# reads an input file with mandatory parameters. the file should contain 4 lines:
# dot bracket structure
# IUPAC sequence
# desired folding energy (in kcal/mol). -1000 means no target
# desired neutrality. -1000 means no target
def read_input_file(input_file_path):
    input_map = {}
    if not os.path.isfile(input_file_path):
        input_map['error'] = "file {} does not exists".format(input_file_path)
    else:
        with open(input_file_path, 'r') as input_file:
            file_errors = []
            structure = bracket_changer(input_file.readline().strip())
            if not is_valid_structure(structure):
                file_errors.append("target structure can hold '.' '(\\{\\[\\<' and ')\\}\\]\\>'" \
                                   " only and must be balanced")
            input_map['target_structure'] = structure
            sequence = input_file.readline().strip()
            if not IUPAC.check_valid_sequence(sequence):
                file_errors.append('target sequence must hold legal IUPAC letters')
            if not len(structure) == len(sequence):
                file_errors.append('sequence and structure must be at the same length')
            input_map['target_sequence'] = sequence.upper().replace('T', 'U')
            try:
                target_energy = float(input_file.readline().strip())
                input_map['target_energy'] = target_energy
            except ValueError:
                file_errors.append('target energy must be a floating number')
            try:
                target_neutrality = float(input_file.readline().strip())
                input_map['target_neutrality'] = target_neutrality
            except ValueError:
                file_errors.append('target neutrality must be a floating number')
            if len(file_errors) > 0:
                input_map['error'] = 'Input file must have 4 lines containing structure, sequence, target energy' \
                                     ' and target neutrality. the following errors were found:\n{}' \
                    .format('\n'.join(file_errors))
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
def read_mandatory_params():
    input_map = {}
    while True:
        structure = bracket_changer(input('Enter structure in dot bracket notation:\n').strip())
        if is_valid_structure(structure):
            input_map['target_structure'] = structure
            break
        print('structure must be balanced and in dot bracket notation')
    while True:
        sequence = input('Enter IUPAC format sequence restrictions (must be same length as fold):\n')
        if IUPAC.check_valid_sequence(sequence) and len(input_map['target_structure']) == len(sequence):
            input_map['target_sequence'] = sequence.upper().replace('T', 'U')
            break
        print('sequence must be in IUPAC format and be the same length as structure')
    while True:
        str_energy = input('Enter the desired minimum free energy in Kcal/mol (-1000 for none):\n')\
            .strip()
        try:
            input_map['target_energy'] = float(str_energy)
            break
        except ValueError:
            print('target energy must by a floating number')
    while True:
        str_neutrality = input('Enter the desired neutrality (-1000 for none):\n')\
            .strip()
        try:
            input_map['target_neutrality'] = float(str_neutrality)
            break
        except ValueError:
            print('target neutrality must by a floating number')
    return input_map


def main(argv):
    logging.info("Starting RNAsfbinv, arguments: {}".format(argv))
    arg_map = generate_arg_map(argv)
    logging.info("Run argument map:\n{}".format(arg_map))
    error = arg_map.get('error')
    if error is not None:
        usage(error)
    else:
        # init RNAfold
        rna_folder = vienna.LiveRNAfold()
        rna_folder.start(arg_map.get('circular'))
        arg_map['RNAfold'] = rna_folder
        # run simulated annealing
        sfb_designer.update_options(arg_map)
        designed_sequence = sfb_designer.simulated_annealing()
        logging.info("Finished simulated annealing, resulting sequence: {}".format(designed_sequence))
        sfb_designer.print_res(designed_sequence)
        rna_folder.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING,
                        format='%(levelname)s:%(asctime)s - %(message)s')
    main(sys.argv[1:])

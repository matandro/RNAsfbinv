#!/usr/bin/env python3
'''
Responsible for running
'''

import os
import re
import sys
import logging
from typing import Dict, List
from subprocess import Popen, PIPE


RNAFOLD_EXE = "RNAfold"
INVERSE_EXE = "RNAinverse"
RES_MATCHER = re.compile(r'(?P<structure>([.()])+) ([({])\s*(?P<energy>[-+]?\d*\.\d+|\d+)\)?.*')
if sys.platform =='win32':
    RNAFOLD_EXE += '.exe'
    INVERSE_EXE += '.exe'


def set_vienna_path(path: str):
    logging.debug('Setting vienna path to ' + path)
    os.environ['VIENNA_PATH'] = path


def output_fold_analyze(output: str) -> Dict[str, str]:
    structure = {}
    try:
        lines = output.split('\n')
        match = RES_MATCHER.match(lines[1].strip())
        if match:
            structure['MFE'] = match.group('structure')
            try:
                structure['MFE_energy'] = float(match.group('energy'))
            except ValueError:
                logging.warning('Could not collect MFE energy: {}'.format(lines[1]))
        else:
            raise Exception(match)
        match = RES_MATCHER.match(lines[3].strip())
        if match:
            structure['centroid'] = match.group('structure')
            try:
                structure['centroid_energy'] = float(match.group('energy'))
            except ValueError:
                logging.warning('Could not collect centroid energy: {}'.format(lines[3]))
        else:
            raise Exception(match)
    except Exception as e:
        logging.warning('could not collect fold data: {}\n{}'.format(str(e), output))
    return structure


END_SEQUENCE = 'ensemble diversity'


class LiveRNAfold:
    def __init__(self):
        self.proc = None

    def __del__(self):
        self.close()

    def _read_until_ready(self) -> List[str]:
        lines = []
        done = False
        line = ''
        while not done:
            c = self.proc.stdout.read(1)
            if c == '\n':
                lines.append(line)
                if END_SEQUENCE in line:
                    done = True
                else:
                    line = ''
            else:
                line += c
        return lines

    def start(self, is_circular: bool=False):
        param_list = [os.path.join(os.getenv('VIENNA_PATH', ""), RNAFOLD_EXE), '-p', '--noPS']
        if is_circular:
            param_list.append('-c')
        logging.debug("Running multi run RNAfold: {}".format(param_list))
        self.proc = Popen(param_list, stdout=PIPE, stdin=PIPE, universal_newlines=True)
        #self._read_until_ready()

    def close(self):
        if self.proc is not None:
            logging.debug("Closing multi run RNAfold")
            self.proc.stdin.write('@\n')
            self.proc.stdin.flush()
            self.proc.stdin.close()
            self.proc.stdout.close()
            self.proc.kill()
            self.proc = None

    def fold(self, sequence: str) -> Dict[str, str]:
        logging.debug("Folding multi run RNAfold: {}".format(sequence))
        write_line = "{}\n".format(sequence)
        self.proc.stdin.write(write_line)
        self.proc.stdin.flush()
        lines = self._read_until_ready()
        structure_map = output_fold_analyze('\n'.join(lines))
        return structure_map


# single use call to RNA fold TODO: add another version that runs on the same client to avoid popen
def fold(sequence: str, is_circular: bool=False, structure_constraints: str = None) -> Dict[str, str]:
    if structure_constraints is not None and len(structure_constraints) != len(sequence):
        return None
    structure_map = None
    try:
        param_list = [os.path.join(os.getenv('VIENNA_PATH', ""), RNAFOLD_EXE), '-p', '--noPS', '-C',
                      '--enforceConstraint']
        if is_circular:
            param_list.append('-c')
        with Popen(param_list, stdout=PIPE, stdin=PIPE, universal_newlines=True) as proc:
            logging.debug("Running RNAfold: {}".format(param_list))
            write_line = "{}\n".format(sequence)
            if structure_constraints is not None:
                write_line += "{}\n".format(structure_constraints)
            else:
                write_line += "{}\n".format('.' * len(sequence))
            write_line += "\n@\n"
            fold_output, errs = proc.communicate(write_line)
            structure_map = output_fold_analyze(fold_output)
    except OSError as e:
        logging.error("Failed to run: '{}'. ERROR: {}".format(param_list, e.errno))
    return structure_map


def output_inverse_analyze(output: str) -> str:
    try:
        lines = output.split('\\n')
        res_sequence = lines[0].split(' ')[0].strip()
    except:
        res_sequence = None
    return res_sequence


def inverse(structure: str, sequence: str=None) -> str:
    res_sequence = None
    try:
        if sequence is None:
            sequence = 'N' * structure
        param_list = [os.path.join(os.getenv('VIENNA_PATH', ""), INVERSE_EXE)]
        with Popen(param_list, stdout=PIPE, stdin=PIPE, universal_newlines=True) as proc:
            logging.debug("Running RNAinverse: {}".format(param_list))
            inverse_output, errs = proc.communicate("{}\n{}\n@\n".format(structure, sequence))
            res_sequence = output_inverse_analyze(inverse_output)
            if res_sequence is not None:
                res_sequence = res_sequence.upper()
            else:
                logging.error("Failed to read output: '{}'. output: {}".format(param_list, inverse_output))
    except OSError as e:
        logging.error("Failed to run: '{}'. ERROR: {}".format(param_list, e.errno))
    return res_sequence


def inverse_seq_ready(sequence: str) -> str:
    res_sequence = ''
    for c in sequence.upper():
        if res_sequence in 'ACGU':
            res_sequence += c.lower()
        else:
            res_sequence += c
    return res_sequence


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_sequence = "AGUAGAUGGCCCGUGGUGUCCCGGAGUGGCUGUAGAGUGAGAUGCAGAUGGAC" \
                    "GACUGAGCCCAUAGGGCCGCUUAUAAUAACCUAUGCCCCCACAUCGUGUAAUU" \
                    "UCAACCCGCAGCACUAUCACAGCCACAGGGUCGAUCA"
    if len(sys.argv) > 1:
        test_sequence = sys.argv[1]

    # TEST RNAfold RUN PER POPEN
    print("Run per popen:\n{}\n".format(fold(test_sequence)))
    test_constraints = '..............................................................((((((........' \
                       '......)))))).......................................................'
    print("Run per popen with constraints:\n{}\n".format(fold(test_sequence, structure_constraints=test_constraints)))

    # TEST RNAfold MULTI RUN PER POPEN
    multi_folder = LiveRNAfold()
    multi_folder.start(False)
    print("Run 1 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence)))
    print("Run 3 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[: int(len(test_sequence) / 2)])))
    print("Run 4 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[int(len(test_sequence) / 2):])))
    multi_folder.close()
    # TEST RNAinverse
    test_structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
    test_sequence = 'NNNNNNNNuNNNNNNNNNNNNNNNNNNNNNNNNuNNNuNNNNNNNNNNNNNNNNNNNNNNyNNNNNNNN'
    print("RNAinverse:\n{}\n".format(inverse(test_structure, test_sequence)))




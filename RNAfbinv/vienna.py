#!/usr/bin/env python3
'''
Responsible for running
'''

import os
import re
import sys
import logging
from subprocess import Popen, PIPE


VIENNA_PATH = "D:\\Programs\\ViennaRNA\\"#""/opt/algorithm/ViennaRNA/bin/"

RNAFOLD_EXE = "RNAfold"
INVERSE_EXE = "RNAinverse"
RES_MATCHER = re.compile(r'(?P<structure>([.()])+) ([({])\s*(?P<energy>[-+]?\d*\.\d+|\d+)\)?.*')
if sys.platform =='win32':
    RNAFOLD_EXE += '.exe'
    INVERSE_EXE += '.exe'


def set_vienna_path(path):
    VIENNA_PATH = path


def output_fold_analyze(output):
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

    def _read_until_ready(self):
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

    def start(self, is_circular=False):
        param_list = [os.path.join(VIENNA_PATH, RNAFOLD_EXE), '-p', '--noPS']
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

    def fold(self, sequence):
        logging.debug("Folding multi run RNAfold: {}".format(sequence))
        self.proc.stdin.write('{}\n'.format(sequence))
        self.proc.stdin.flush()
        lines = self._read_until_ready()
        structure_map = output_fold_analyze('\n'.join(lines))
        return structure_map


# single use call to RNA fold TODO: add another version that runs on the same client to avoid popen
def fold(sequence, is_circular=False):
    structure_map = None
    try:
        param_list = [os.path.join(VIENNA_PATH, RNAFOLD_EXE), '-p', '--noPS']
        if is_circular:
            param_list.append('-c')
        with Popen(param_list, stdout=PIPE, stdin=PIPE, universal_newlines=True) as proc:
            logging.debug("Running RNAfold: {}".format(param_list))
            fold_output, errs = proc.communicate("{}\n@\n".format(sequence))
            structure_map = output_fold_analyze(fold_output)
    except OSError as e:
        logging.error("Failed to run: '{}'. ERROR: {}".format(param_list, e.errno))
    return structure_map


def output_inverse_analyze(output):
    try:
        lines = output.split('\\n')
        res_sequence = lines[0].split(' ')[0].strip()
    except:
        res_sequence = None
    return res_sequence


def inverse(structure, sequence=None):
    res_sequence = None
    try:
        if sequence is None:
            sequence = 'N' * structure
        param_list = [os.path.join(VIENNA_PATH, INVERSE_EXE)]
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


def inverse_seq_ready(sequence):
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
    test_structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
    # TEST RNAfold MULTI RUN PER POPEN
    multi_folder = LiveRNAfold()
    multi_folder.start(False)
    print("Run 1 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence)))
    print("Run 2 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[: int(len(test_sequence) / 2)])))
    print("Run 3 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[int(len(test_sequence) / 2):])))
    multi_folder.close()
    # TEST RNAinverse
    test_sequence = 'NNNNNNNNuNNNNNNNNNNNNNNNNNNNNNNNNuNNNuNNNNNNNNNNNNNNNNNNNNNNyNNNNNNNN'
    print("RNAinverse:\n{}\n".format(inverse(test_structure, test_sequence)))




#!/usr/bin/env python3
'''
Responsible for running
'''

import os
import re
import sys
import logging
from typing import Dict, List
from subprocess import Popen, PIPE, DEVNULL
from threading import Thread
from queue import Queue


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
    def __init__(self, logger=None):
        self.proc = None
        self.read_lines = Queue()
        self.reader_thread = None
        self.logger = logger

    def __del__(self):
        self.close()

    def _read_until_ready(self) -> List[str]:
        lines = []
        done = False
        while not done:
            line = self.read_lines.get()
            if self.logger is not None:
                self.logger.debug("LiveRNAfold, _read_until_ready ({}) [{}]".format(self.read_lines.qsize(), line))
            lines.append(line.strip())
            if END_SEQUENCE in line:
                done = True
        return lines

    def start(self, is_circular: bool=False):
        def _reader_func(stream: PIPE, queue: Queue, logger: logging):
            try:
                while True:
                    line = stream.readline()
                    if len(line) > 0:
                        queue.put(line)
                        if logger is not None:
                            logger.debug("LiveRNAfold, _reader_func ({}) [{}]".format(queue.qsize(), line))
                    else:
                        break
            except Exception as e:
                self.logger.debug("LiveRNAfold, _reader_func THREAD DEAD [{}]".format(e))
                queue.put(END_SEQUENCE)

        param_list = [os.path.join(os.getenv('VIENNA_PATH', ""), RNAFOLD_EXE), '-p', '--noPS']
        if is_circular:
            param_list.append('-c')
        if self.logger is not None:
            self.logger.debug("LiveRNAfold, start: {}".format(param_list))
        self.proc = Popen(param_list, stdout=PIPE, stdin=PIPE, stderr=DEVNULL, universal_newlines=True)
        self.reader_thread = Thread(target=_reader_func, args=(self.proc.stdout, self.read_lines, self.logger))
        self.reader_thread.daemon = True
        self.reader_thread.start()

    def close(self):
        if self.proc is not None:
            if self.logger is not None:
                self.logger.debug("LiveRNAfold, close")
            self.proc.stdin.write('@\n')
            self.proc.stdin.flush()
            self.proc.kill()
            self.proc.stdin.close()
            self.proc.stdout.close()
            self.proc = None
            self.reader_thread = None

    def fold(self, sequence: str) -> Dict[str, str]:
        write_line = "{}\n".format(sequence)
        if self.logger is not None:
            self.logger.debug("LiveRNAfold, writing [\n{}\n]".format(write_line))
        self.proc.stdin.write(write_line)
        self.proc.stdin.flush()
        lines = self._read_until_ready()
        structure_map = output_fold_analyze('\n'.join(lines))
        return structure_map


# single use call to RNA fold
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
            sequence = 'N' * len(structure)
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


def inverse_seq_ready(target_sequence: str, start_sequence: str=None) -> str:
    if start_sequence is None or start_sequence == '':
        start_sequence = 'N'*len(target_sequence)
    elif len(start_sequence) != len(target_sequence):
        raise ValueError("start_sequence must be empty \ None or the length of target_sequence length")
    res_sequence = ''
    for i in range(0, len(target_sequence)):
        c = target_sequence[i]
        if c in 'ACGU':
            res_sequence += c.lower()
        else:
            res_sequence += start_sequence[i]
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
    test_constraints = '.......((...........)).....(((....................' \
                       '................((((.................)))).........' \
                       '...........))).............................'
    print("Run per popen with constraints:\n{}\n".format(fold(test_sequence, structure_constraints=test_constraints)))

    # TEST RNAfold MULTI RUN PER POPEN

    multi_folder = LiveRNAfold()
    multi_folder.start(False)
    print("Run 1 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence)))
    print("Run 2 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[: int(len(test_sequence) / 2)])))
    print("Run 3 multi run per popen:\n{}\n".format(multi_folder.fold(test_sequence[int(len(test_sequence) / 2):])))
    multi_folder.close()
    # TEST RNAinverse
    # test_structure = '((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))'
    # test_sequence = 'NNNNNNNNuNNNNNNNNNNNNNNNNNNNNNNNNuNNNuNNNNNNNNNNNNNNNNNNNNNNyNNNNNNNN'
    test_structure = '..((((((((......(((.......))).....((((......))))...........................(((((.......))))).....(((.......))).......))))))))..'
    test_sequence = 'NNNNNNNYUCNGGGNNNGGNGNNNNUCCNNANCGGNNGUNNAGNNCNNGANNNNNNNNNNNNNNNNNNNNNNNGANNNNNNNNNNNNNNNNNRNCGANRGNNANAGUCYNGAUNNNARANNNNNNNN'
    test_sequence = inverse_seq_ready(test_sequence, test_sequence)
    print("RNAinverse:\n{}\n".format(inverse(test_structure, test_sequence)))

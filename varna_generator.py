#!/usr/bin/env python3
'''
enter description
'''


from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen
import logging
import os


JAVA_PATH = 'java'
VARNA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib', 'VARNAv3-93.jar')


def generate_temp_ct(structure, sequence, title=''):
    comp_map = {}
    closing_stack = []
    for i in range(0, len(sequence)):
        if structure[i] == '(':
            closing_stack.append(i)
        elif structure[i] == ')':
            start_index = closing_stack.pop()
            comp_map[start_index] = i
            comp_map[i] = start_index
    temp_file = NTF(dir='.', delete=False, suffix='.ct', mode='w')
    temp_file.write('{}\t{}\n'.format(len(structure), title))
    for i in range(0, len(sequence)):
        next_item = i + 2
        if next_item > len(sequence):
            next_item = 0
        comp = comp_map.get(i)
        if comp is None:
            comp = 0
        else:
            comp += 1
        line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(i + 1, sequence[i], i, next_item, comp, i)
        temp_file.write(line)
    temp_file.close()
    return temp_file.name


def call_varna(structure, sequence, out_file_path, index_str=None, title=''):
    res = False
    ct_file = None
    try:
        ct_file = generate_temp_ct(structure, sequence, title)
        param_list = [JAVA_PATH, '-cp', VARNA_PATH, 'fr.orsay.lri.varna.applications.VARNAcmd',
                      '-i', ct_file, '-o', out_file_path, '-resolution', '5.0']
        if index_str is not None:
            param_list = param_list + ['-basesStyle1', 'fill=#00FF00,outline=#FF0000', '-applyBasesStyle1on', index_str]
        logging.debug("Generating VARNA image {}".format(param_list))
        with Popen(param_list) as proc:
            proc.wait()
    finally:
        if ct_file is not None and os.path.exists(ct_file):
            os.remove(ct_file)
        f_stat = os.stat(out_file_path)
        if f_stat is not None and f_stat.st_size > 0:
            res = True
    return res


def generate_image(structure, sequence, index_list=None, title=None, output_file_path=None):
    def generate_marked_indexes(ordered_list):
        index_str = None
        if ordered_list is not None:
            index_str = ''
            ordered_list.sort()
            low = None
            last = None
            for item in ordered_list:
                if low is None:
                    low = item
                    last = item
                elif item == last + 1:
                    last = item
                else:
                    index_str += ',' + str(low)
                    if last > low:
                        index_str += '-' + str(last)
                    low = item
                    last = item
            if low is not None:
                index_str += ',' + str(low)
                if last > low:
                    index_str += '-' + str(last)
            if index_str[0] == ',':
                index_str = index_str[1:]
        return index_str
    res = None
    image_name = output_file_path
    if image_name is None:
        image_file = NTF(dir='.', delete=False, suffix='.jpg', mode='w')
        image_file.close()
        image_name = image_file.name
    if call_varna(structure, sequence, image_name, title=title,
                  index_str=generate_marked_indexes(None if index_list is None else sorted(index_list))):
        res = image_name
    return res


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    structure_ = '...(((...(((....))....)))).....((....(((....))))...)...'
    #             1234567890123456789012345678901234567890123456789012345
    #             0000000001111111111222222222233333333334444444444555555
    sequence__ = 'GAUAGGGAGCAGUGACCUUCGAUCCGCGCGCGUAAUAUGCGUCAGCAACGACUAG'
    print(generate_image(structure_, sequence__, output_file_path='TestVARNA.jpg', title='No inedx list'))
    index_list_ = [4,26,5,6,24,43,44,25,41,42]
    print(generate_image(structure_, sequence__, output_file_path='TestVARNAIndex.jpg', title='With inedx list',
                         index_list=index_list_))


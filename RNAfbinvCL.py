#!/usr/bin/env python3
'''
A simple code to activate RNAfbinv 2.0 from command line
'''

from rnafbinv import RNAfbinvCL, vienna
import sys
import os
import configparser
import varna_generator


def read_config():
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_path)
    path_section = config['PATH']
    vienna.set_vienna_path(path_section.get('VIENNA', ''))
    java = path_section.get('JAVA')
    if java is not None and java != '':
        varna_generator.set_java_path(java)
    varna = path_section.get('VARNA')
    if varna is not None and varna != '':
        varna_generator.set_varna_path(varna)


read_config()
result = RNAfbinvCL.main(' '.join(sys.argv[1:]))


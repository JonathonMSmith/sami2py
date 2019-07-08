#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2017, JK
# Full license can be found in License.md
# -----------------------------------------------------------------------------

from __future__ import print_function
import sys
from os import path, mkdir
from setuptools import setup
import subprocess

# generate path for fortran model files
here = path.abspath(path.dirname(__file__))
fortran_path = path.join(here, 'sami2py', 'fortran')
test_data_path = path.join(here, 'sami2py', 'tests', 'test_data')
file_path = path.join(sys.prefix, '.sami2py')

# %% build


if not path.isfile(path.join(fortran_path, 'sami2py.x')):
    try:  # py27 does not have shutil.which()
        cmd = ['gfortran', '-fno-range-check', '-fno-automatic',
               '-ffixed-line-length-none', '-o', 'sami2py.x']
        src = ['nrlmsise00_modified.f', 'grid-1.00.f', 'sami2py-1.00.f',
               'hwm93.f', 'hwm07e_modified.f90', 'apexcord.f90', 'hwm14.f90']
        subprocess.call(cmd + src, cwd=fortran_path)
    except OSError:
        pass

if not path.isfile(path.join(fortran_path, 'sami2py.x')):
    print('\nYou will need to compile the fortran files.  Try\n'
          '$  make -C sami2py/fortran compile\n', file=sys.stderr)

if not path.isdir(file_path):
    mkdir(file_path)
    print(''.join(('Created .sami2py directory in ' + sys.prefix + ' to',
                   'store settings.')))

with open(path.join(file_path, 'fortran_path.txt'), 'w+') as f:
    f.write(fortran_path)
with open(path.join(file_path, 'test_data_path.txt'), 'w+') as f:
    f.write(test_data_path)

setup()

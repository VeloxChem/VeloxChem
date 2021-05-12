#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from pathlib import PurePath
import numpy as np
import sys
import re

from .errorhandler import assert_msg_critical


class InputParser:
    """
    Implements the input parser and provides functions for parsing VeloxChem
    input files into a format, which passes the needed information to the rest
    of the program.

    :param inpname:
        The name of the input file.
    :param outname:
        The name of the output file.

    Instance variables
        - input_dict: The input dictionary.
        - inpname: The name of the input file.
        - outname: The name of the output file.
        - is_basis_set: The flag for parsing a basis set file.
        - basis_set_name: The name of the basis set.
    """

    def __init__(self, inpname, outname=None):
        """
        Initializes the parser and parses the input file.
        """

        self.input_dict = {}

        self.inpname = inpname if isinstance(inpname, str) else str(inpname)
        self.outname = outname

        self.is_basis_set = False
        self.basis_set_name = ''

        self.parse()

    def parse(self):
        """
        Parses the input file.
        """

        self.input_dict = {}

        errmsg = f'InputParser: bad syntax in file {self.inpname}. '
        errmsg += 'You may check for incorrect, incomplete or empty groups.'

        group = None

        with open(str(self.inpname), 'r') as f_inp:

            for line in f_inp:

                # remove comment and extra space
                line = line.strip()
                line = re.sub(r'!.*', '', line)
                line = ' '.join(line.split())

                # skip empty line
                if not line:
                    continue

                # skip first line if reading basis set
                if (not self.is_basis_set) and (line[:10] == '@BASIS_SET'):
                    self.is_basis_set = True
                    self.basis_set_name = line.split()[1]
                    continue

                # begin of group
                if line[0] == '@' and line.lower() != '@end':
                    assert_msg_critical(group is None, errmsg)
                    group = '_'.join(line[1:].strip().lower().split())
                    if self.is_basis_set:
                        self.input_dict[group] = []
                    else:
                        self.input_dict[group] = {}

                # end of group
                elif line.lower() == '@end':
                    assert_msg_critical(group is not None, errmsg)
                    group = None

                # inside a group
                else:

                    # for basis set
                    if self.is_basis_set:
                        if group is not None:
                            self.input_dict[group].append(line)
                        continue

                    # for input file
                    assert_msg_critical(group is not None, errmsg)

                    if ':' in line:
                        key, value = line.split(':')
                        key = '_'.join(key.strip().lower().split())
                        value = value.strip()
                        if value:
                            # single-line input
                            # key and value are in the same line
                            self.input_dict[group][key] = value
                        else:
                            # multi-line input
                            # key is followed by multiple lines
                            self.input_dict[group][key] = []
                    else:
                        # multi-line input
                        # key is followed by multiple lines
                        self.input_dict[group][key].append(line)

        # save basis set name
        if self.is_basis_set:
            self.input_dict['basis_set_name'] = self.basis_set_name
            return

        # check incomplete group
        assert_msg_critical(group is None, errmsg)

        # check empty group
        for group in self.input_dict:
            assert_msg_critical(len(self.input_dict[group]), errmsg)

        # save filename
        if self.outname not in [None, sys.stdout]:
            f_path = PurePath(self.outname)
        else:
            f_path = PurePath(self.inpname)
        self.input_dict['filename'] = str(f_path.with_name(f_path.stem))

        # update checkpoint filename
        abbrev = {
            'scf': 'scf',
            'response': 'rsp',
            'exciton': 'exciton',
            'loprop': 'loprop',
        }

        for job_type in abbrev:
            if job_type not in self.input_dict:
                self.input_dict[job_type] = {}
            if 'checkpoint_file' not in self.input_dict[job_type]:
                checkpoint_file = str(
                    f_path.with_suffix(f'.{abbrev[job_type]}.h5'))
                self.input_dict[job_type]['checkpoint_file'] = checkpoint_file

        # verify options for loprop
        verifyers = {
            'loprop': {
                'checkpoint_file': (lambda v: True, 'Always OK'),
                'localize':
                    (lambda v: v in ['charges'], 'localize: {} illegal value')
            }
        }

        if 'loprop' in self.input_dict:
            for option, value in self.input_dict['loprop'].items():
                verify, msg = verifyers['loprop'][option]
                if not verify(value):
                    raise InputError(msg.format(value))

    def get_dict(self):
        """
        Gets the input dictonary.

        :return:
            A dictionary containing all information form the input file.
        """

        return self.input_dict


class InputError(Exception):
    pass


def parse_seq_fixed(input_seq, flag='float'):
    """
    Parses input sequence that has fixed number of elements, using comma and/or
    white space as delimiter.

    :param input_seq:
        The input sequence.
    :param flag:
        The flag for int/float.

    :return:
        A tuple.
    """

    if isinstance(input_seq, (tuple, list, np.ndarray)):
        return tuple(input_seq)

    elif isinstance(input_seq, str):
        seq_str = input_seq.replace(',', ' ')
        if flag == 'float':
            return tuple([float(x) for x in seq_str.split()])
        elif flag == 'int':
            return tuple([int(x) for x in seq_str.split()])
        else:
            assert_msg_critical(False, f'parse_seq_fixed: invalid flag {flag}')

    else:
        errmsg = 'parse_seq_fixed: invalid type for input_seq'
        assert_msg_critical(False, errmsg)


def parse_seq_range(input_seq):
    """
    Parses input sequence that is specified by range.
    Input example: '0.0 - 0.2525 (0.0025), 0.5 - 1.0 (0.02), 2.0'

    :param input_seq:
        The input sequence.

    :return:
        A tuple.
    """

    if isinstance(input_seq, (tuple, list, np.ndarray)):
        return tuple(input_seq)

    elif isinstance(input_seq, str):
        seq = []

        for w in input_seq.split(','):
            if not w:
                continue

            m = re.search(r'^(.*)-(.*)\((.*)\)$', w)
            if m is None:
                m = re.search(r'^(.*)-(.*)-(.*)$', w)

            if m is not None:
                assert_msg_critical(m is not None,
                                    'parse_seq_range: failed to read seqence')

                start, stop, step = (float(m.group(1)), float(m.group(2)),
                                     float(m.group(3)))
                stop += 0.01 * step

                seq += list(np.arange(start, stop, step))

            else:
                seq.append(float(w))

        return tuple(seq)

    else:
        errmsg = 'parse_seq_range: invalid type for input_seq'
        assert_msg_critical(False, errmsg)


def parse_bool(input_bool):
    """
    Parses input boolean.

    :param input_bool:
        The input.

    :return:
        A boolean.
    """

    if isinstance(input_bool, bool):
        return input_bool

    elif isinstance(input_bool, str):
        bool_val = True if input_bool.lower() in ['yes', 'y'] else False
        return bool_val

    else:
        errmsg = 'parse_bool: invalid type for input_bool'
        assert_msg_critical(False, errmsg)


def parse_str(input_str, flag=None):
    """
    Parses input string.

    :param input_str:
        The input.
    :param flag:
        The flag for upper/lower case.

    :return:
        A string.
    """

    err_str = f'parse_str: Expecting a string but got {input_str}'
    assert_msg_critical(isinstance(input_str, str), err_str)

    if flag is None:
        return input_str
    elif flag == 'upper':
        return input_str.upper()
    elif flag == 'lower':
        return input_str.lower()
    else:
        assert_msg_critical(False, f'parse_str: invalid flag {flag}')


def parse_input(obj, keyword_types, input_dictionary):
    """
    Parses input keywords for object.

    :param obj:
        The object.
    :param keyword_types:
        The data type associated with keyword.
    :param input_dictionary:
        The input dictionary.
    """

    for key, val in input_dictionary.items():
        if key not in keyword_types:
            continue

        # the value of certain keyword can be set to None
        if val is None:
            setattr(obj, key, val)

        elif keyword_types[key] == 'str':
            setattr(obj, key, parse_str(val))

        elif keyword_types[key] == 'str_upper':
            setattr(obj, key, parse_str(val, 'upper'))

        elif keyword_types[key] == 'str_lower':
            setattr(obj, key, parse_str(val, 'lower'))

        elif keyword_types[key] == 'int':
            setattr(obj, key, int(val))

        elif keyword_types[key] == 'float':
            setattr(obj, key, float(val))

        elif keyword_types[key] == 'bool':
            setattr(obj, key, parse_bool(val))

        elif keyword_types[key] == 'seq_fixed_int':
            setattr(obj, key, parse_seq_fixed(val, 'int'))

        elif keyword_types[key] == 'seq_fixed':
            setattr(obj, key, parse_seq_fixed(val))

        elif keyword_types[key] == 'seq_range':
            setattr(obj, key, parse_seq_range(val))

        else:
            err_type = f'parse_input: invalid keyword type for {key}'
            assert_msg_critical(False, err_type)

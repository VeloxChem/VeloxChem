import numpy as np
import sys
import re
import os

from .errorhandler import assert_msg_critical


class InputParser:
    """
    Implements the input parser and provides functions for parsing VeloxChem
    input files into a format, which passes the needed information to the rest
    of the program.

    :param filename:
        The name of the input file.
    :param outname:
        The name of the output file.

    Instance variables
        - input_dict: The input dictionary.
        - success_monitor: The monitor for successful parsing.
        - filename: The name of the input file.
        - is_basis_set: The flag for parsing a basis set file.
        - basis_set_name: The name of the basis set.
    """

    def __init__(self, filename, outname=None):
        """
        Initializes the parser and parses the input file.
        """

        self.input_dict = {}
        self.success_monitor = True

        self.filename = filename
        self.outname = outname
        self.is_basis_set = False
        self.basis_set_name = ''

        self.parse()
        self.consistency_checks()

    # defining main functions

    def parse(self):
        """
        Calls every function needed for the parsing process depending on
        the success of the parsing in different stages.
        """

        errmsg = ''

        try:

            # reading selected file

            self.read_file()

            # checking for syntax correctness of the input file

            self.incomplete_group_check()
            self.empty_group_check()

        except FileNotFoundError:
            errmsg = 'InputParser: cannot open file {}'.format(self.filename)
            self.success_monitor = False

        except SyntaxError:
            errmsg = 'InputParser: bad syntax in file '
            errmsg += self.filename + os.linesep + '     '
            errmsg += 'You may check for incorrect, '
            errmsg += 'incomplete or empty groups.'
            self.success_monitor = False

        assert_msg_critical(self.success_monitor, errmsg)

        if self.success_monitor:

            # manipulation of input string

            self.clear_interspace()

            # processing the data into lists and dictionaries

            self.groupsplit()
            self.convert_dict()

    def get_dict(self):
        """
        Gets the input dictonary.

        :return:
            A dictionary containing all information form the input file.
        """

        return self.input_dict

    def read_file(self):
        """
        Reads and stores content of selected file as a string type. Deletes
        comments (marked by '!') and unnecassary whitespace.
        """

        self.content = ''
        with open(self.filename, 'r') as f_inp:
            for line in f_inp:

                # remove comment and extra white spaces
                line = line.strip()
                line = re.sub(r'!.*', '', line)
                line = ' '.join(line.split())

                # skip first line if reading basis set
                if line[:10] == '@BASIS_SET':
                    self.is_basis_set = True
                    self.basis_set_name = line.split()[1]
                    continue

                # take care of end of group
                if line.lower()[:4] == '@end':
                    line = line.lower()

                # add trailing os.linesep
                if line:
                    self.content += line + os.linesep

    def incomplete_group_check(self):
        """
        Checks for any incomplete groups.
        """

        if re.findall(r'@(?!end)[^@]*@(?!end)|@end[^@]*@end',
                      self.content) != []:
            raise SyntaxError

        last_lines = self.content.split('@end')[-1]
        if re.findall(r'@(?!end)[^@]*', last_lines) != []:
            raise SyntaxError

    def empty_group_check(self):
        """
        Checks for any empty groups.
        """

        if re.findall(r'@\w[\w ]*' + os.linesep + r'\s*@end',
                      self.content) != []:
            raise SyntaxError

    def clear_interspace(self):
        """
        Deletes content that is not within a group.
        """

        self.content = re.sub(r'@end[^@]*@', '@end' + os.linesep + '@',
                              self.content)

    def groupsplit(self):
        """
        Creates a list in which every element is a list itself containing
        every line of a group, while deleting '@' and '@end' tags.
        """

        self.grouplist = re.findall(r'@(?!end)[^@]*@end', self.content)
        for i in range(len(self.grouplist)):
            self.grouplist[i] = self.grouplist[i].strip().replace('@', '')
            self.grouplist[i] = self.grouplist[i].split(os.linesep)[:-1]

    def convert_dict(self):
        """
        Converting the list of lists into a dictionary with groupnames as
        keys and group content as a dictionary itself. The geometry definition
        of the molecule group is stored in a different dictionary. Converting
        the molecular structure into the required format.
        """

        for group in self.grouplist:
            local_dict = {}
            local_list = []

            xyz_flag, constr_flag = False, False

            for entry in group[1:]:
                if ':' in entry:
                    # save input settings in local_dict
                    # except for xyz strings and atomic basis functions
                    key = entry.split(':')[0].strip()
                    key = '_'.join(key.split())
                    if key.lower() == 'xyz':
                        xyz_flag, constr_flag = True, False
                    elif key.lower() == 'constraints':
                        xyz_flag, constr_flag = False, True
                        local_dict['constraints'] = []
                    else:
                        xyz_flag, constr_flag = False, False
                        val = entry.split(':')[1].strip()
                        local_dict[key.lower()] = val
                elif constr_flag:
                    # save constraints
                    local_dict['constraints'].append(entry)
                elif xyz_flag:
                    # save xyz strings
                    local_list.append(entry)
                else:
                    # save atomic basis functions
                    local_list.append(entry)

            group_key = group[0]
            group_key = '_'.join(group_key.split())

            if self.is_basis_set:
                # for basis set, local_list contains atomic basis functions
                self.input_dict[group_key.lower()] = local_list
            else:
                # for input file, local_dict contains input settings
                self.input_dict[group_key.lower()] = local_dict
                if group_key.lower() == 'molecule':
                    # local_list contains xyz strings of the molecule
                    xyzstr = os.linesep.join(local_list)
                    self.input_dict['molecule']['xyzstr'] = xyzstr
                    # also set the default value for units
                    if 'units' not in self.input_dict['molecule']:
                        self.input_dict['molecule']['units'] = 'angs'

        if self.is_basis_set:
            # for basis set, save basis set file name
            self.input_dict['basis_set_name'] = self.basis_set_name

        else:
            # process filename
            if self.outname not in [None, sys.stdout]:
                fname = self.outname
            else:
                fname = self.filename
            if '.' in fname:
                fname = '.'.join(fname.split('.')[:-1])
            self.input_dict['filename'] = fname

            # update checkpoint filename
            fchkp = fname + '.scf.h5'
            frsp = fname + '.rsp.h5'
            fexciton = fname + '.exciton.h5'
            floprop = fname + '.loprop.h5'

            if 'scf' not in self.input_dict:
                self.input_dict['scf'] = {}
            if 'checkpoint_file' not in self.input_dict['scf']:
                self.input_dict['scf']['checkpoint_file'] = fchkp

            if 'response' not in self.input_dict:
                self.input_dict['response'] = {}
            if 'checkpoint_file' not in self.input_dict['response']:
                self.input_dict['response']['checkpoint_file'] = frsp

            if 'exciton' not in self.input_dict:
                self.input_dict['exciton'] = {}
            if 'checkpoint_file' not in self.input_dict['exciton']:
                self.input_dict['exciton']['checkpoint_file'] = fexciton

            if 'loprop' not in self.input_dict:
                self.input_dict['loprop'] = {}
            if 'checkpoint_file' not in self.input_dict['loprop']:
                self.input_dict['loprop']['checkpoint_file'] = floprop

    def consistency_checks(self):
        """
        Checks consistency in the input groups.
        """

        if self.is_basis_set:
            return

        self.verify_options('loprop')

    def verify_options(self, group):
        """
        Detect input errors for selected input group.
        Checks and messages defined by `verifiers` dict.

        :param group:
            The input group.
        """

        verifyers = {
            'loprop': {
                'checkpoint_file': (lambda v: True, 'Always OK'),
                'localize':
                    (lambda v: v in ['charges'], 'localize: {} illegal value')
            }
        }

        if group in self.input_dict:
            for option, value in self.input_dict[group].items():
                verify, msg = verifyers[group][option]
                if not verify(value):
                    raise InputError(msg.format(value))

    @staticmethod
    def parse_frequencies(input_frequencies):
        """
        Parses frequencies input for response solver.
        Input example: '0.0 - 0.2525 (0.0025), 0.5 - 1.0 (0.02), 2.0'

        :param input_frequencies:
            The string of input frequencies.
        """
        if isinstance(input_frequencies, list):
            return input_frequencies

        if isinstance(input_frequencies, np.ndarray):
            return input_frequencies.tolist()

        frequencies = []
        for w in input_frequencies.split(','):
            if '-' in w:
                m = re.search(r'^(.*)-(.*)\((.*)\)$', w)
                if m is None:
                    m = re.search(r'^(.*)-(.*)-(.*)$', w)

                assert_msg_critical(m is not None,
                                    'InputParser: failed to read frequencies')

                frequencies += list(
                    np.arange(
                        float(m.group(1)),
                        float(m.group(2)),
                        float(m.group(3)),
                    ))
            elif w:
                frequencies.append(float(w))
        return frequencies


class InputError(Exception):
    pass

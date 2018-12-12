from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import assert_msg_critical

import re


class InputParser:
    """ Provides functions for parsing VeloxChem input files into a format,
    which passes the needed information to the rest of the program """

    def __init__(self, filename):
        """ Initializes parsing process by setting monitor value to TRUE """

        self.success_monitor = True

        self.filename = filename
        self.is_basis_set = False
        self.basis_set_name = ''

    # defining main functions

    def parse(self):
        """ Calls every function needed for the parsing process depending on
        the success of the parsing in different stages """

        try:

            # reading selected file

            self.read_file()

            # checking for syntax correctness of the input file

            self.incomp_group_check()
            self.empty_group_check()

        except FileNotFoundError:
            print('Selected file was not found!')
            self.success_monitor = False

        except SyntaxError:
            print('Selected input file has bad syntax!')
            print('You may check for incomplete or empty groups.')
            self.success_monitor = False

        if self.success_monitor:

            # manipulation of input string

            self.clear_interspace()

            # processing the data into lists and dictionaries

            self.groupsplit()
            self.convert_dict()

            if self.is_basis_set:
                return self.input_dict

            # converting angstroms to atomic units, if needed

            need_convert_units = True

            if 'molecule' in self.input_dict.keys():
                if 'units' in self.input_dict['molecule'].keys():
                    units = self.input_dict['molecule']['units']
                    if units in ['au', 'bohr', 'bohrs']:
                        need_convert_units = False

            if need_convert_units:
                self.convert_units()

            return self.input_dict

    def parse_success(self):
        """ Performing the parsing process and returning monitor value. """

        self.parse()
        return self.success_monitor

    # defining subordinated functions

    def read_file(self):
        """ Storing content of selected file as a string type. Deleting
        comments (marked by '!'). Deleting unnecassary whitespace. """

        self.content = ''
        with open(self.filename, 'r') as f_inp:
            for line in f_inp:

                # remove comment and extra white spaces
                line = line.strip()
                line = re.sub('!.*', '', line)
                line = ' '.join(line.split())

                # skip first line if reading basis set
                if line[:10] == '@BASIS_SET':
                    self.is_basis_set = True
                    self.basis_set_name = line.split()[1]
                    continue

                # take care of end of group
                if line.lower()[:4] == '@end':
                    line = line.lower()

                # add trailing '\n'
                if line:
                    self.content += line + '\n'

    def incomp_group_check(self):
        """ Checking for any incomplete groups. """

        if re.findall(
                r'@(?!end[\s])[^@]*@(?!end(?![\w]))|@end\s[^@]*@end(?![\w])',
                self.content) != []:
            raise SyntaxError

    def empty_group_check(self):
        """ Checking for any empty groups. """

        if re.findall(r'@[\w ]*\n\s*@end(?![\w])', self.content) != []:
            raise SyntaxError

    def clear_interspace(self):
        """ Deleting content, that's not within a group. """

        self.content = re.sub(r'@end\s[^@]*@', '@end\n@', self.content)

    def groupsplit(self):
        """ Creating a list in which every element is a list itself containing
        every line of a group, while deleting '@' and '@end' tags. """

        self.grouplist = re.findall(r'@(?!end[\s])\s*\w+[^@]*@end(?![\w])',
                                    self.content)
        for i in range(len(self.grouplist)):
            self.grouplist[i] = self.grouplist[i].lstrip('@ \n')
            self.grouplist[i] = self.grouplist[i].rstrip('end\n')
            self.grouplist[i] = self.grouplist[i].rstrip('\n@')
            self.grouplist[i] = self.grouplist[i].split('\n')

    def convert_dict(self):
        """ Converting the list of lists into a dictionary with groupnames as
        keys and group content as a dictionary itself. The geometry definition
        of the molecule group is stored in a different dictionary. Converting
        the molecular structure into the required format. """

        self.input_dict = {}
        for group in self.grouplist:
            local_dict = {}
            local_list = []

            for entry in group[1:]:
                if ':' in entry:
                    key = entry.split(':')[0].strip()
                    key = '_'.join(key.split())
                    val = entry.split(':')[1].strip()
                    if key.lower() != 'xyz':
                        local_dict[key.lower()] = val
                else:
                    local_list.append(entry)

            group_key = group[0]
            group_key = '_'.join(group_key.split())

            if self.is_basis_set:
                self.input_dict[group_key.lower()] = local_list
            else:
                self.atom_list = local_list
                self.input_dict[group_key.lower()] = local_dict

        if self.is_basis_set:
            self.input_dict['basis_set_name'] = self.basis_set_name
            return

        self.input_dict['molecule']['atom_labels'] = []
        self.input_dict['molecule']['x_coords'] = []
        self.input_dict['molecule']['y_coords'] = []
        self.input_dict['molecule']['z_coords'] = []
        for atom in self.atom_list:
            axyz = atom.split()
            self.input_dict['molecule']['atom_labels'].append(axyz[0])
            self.input_dict['molecule']['x_coords'].append(float(axyz[1]))
            self.input_dict['molecule']['y_coords'].append(float(axyz[2]))
            self.input_dict['molecule']['z_coords'].append(float(axyz[3]))

    def convert_units(self):
        """ Converting molecule coordinates from angstroms to atomic units. """

        coords = ['x_coords', 'y_coords', 'z_coords']
        angstroms_in_bohr = 1.0 / bohr_in_angstroms()
        for n in coords:
            for p in range(len(self.input_dict['molecule'][n])):
                self.input_dict['molecule'][n][p] *= angstroms_in_bohr

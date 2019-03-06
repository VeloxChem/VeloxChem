from .errorhandler import assert_msg_critical

import re


class InputParser:
    """ Provides functions for parsing VeloxChem input files into a format,
    which passes the needed information to the rest of the program """

    def __init__(self, filename):
        """ Initializes the parser and parses the input file """

        self.input_dict = {}
        self.success_monitor = True

        self.filename = filename
        self.is_basis_set = False
        self.basis_set_name = ''

        self.parse()

    # defining main functions

    def parse(self):
        """ Calls every function needed for the parsing process depending on
        the success of the parsing in different stages """

        errmsg = ''

        try:

            # reading selected file

            self.read_file()

            # checking for syntax correctness of the input file

            self.incomp_group_check()
            self.empty_group_check()

        except FileNotFoundError:
            errmsg = 'input parser: cannot open file %s' % self.filename
            self.success_monitor = False

        except SyntaxError:
            errmsg = 'input parser: bad syntax in file %s' % self.filename
            errmsg += '\n     You may check for incomplete or empty groups.'
            self.success_monitor = False

        assert_msg_critical(self.success_monitor, errmsg)

        if self.success_monitor:

            # manipulation of input string

            self.clear_interspace()

            # processing the data into lists and dictionaries

            self.groupsplit()
            self.convert_dict()

    def get_dict(self):
        """ Gets the input dictonary """

        return self.input_dict

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

        last_lines = self.content.split('@end')[-1].split('\n')
        for line in last_lines:
            line = line.strip()
            if len(line) > 0 and line[0] == '@':
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
            self.grouplist[i] = self.grouplist[i].strip().replace('@', '')
            self.grouplist[i] = self.grouplist[i].split('\n')[:-1]

    def convert_dict(self):
        """ Converting the list of lists into a dictionary with groupnames as
        keys and group content as a dictionary itself. The geometry definition
        of the molecule group is stored in a different dictionary. Converting
        the molecular structure into the required format. """

        for group in self.grouplist:
            local_dict = {}
            local_list = []

            for entry in group[1:]:
                if ':' in entry:
                    # save input settings in local_dict
                    # except for xyz strings and atomic basis functions
                    key = entry.split(':')[0].strip()
                    key = '_'.join(key.split())
                    val = entry.split(':')[1].strip()
                    if key.lower() != 'xyz':
                        local_dict[key.lower()] = val
                else:
                    # save xyz strings or atomic basis functions in local_list
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
                    xyzstr = '\n'.join(local_list)
                    self.input_dict['molecule']['xyzstr'] = xyzstr
                    # also set the default value for units
                    if 'units' not in self.input_dict['molecule'].keys():
                        self.input_dict['molecule']['units'] = 'angs'

        if self.is_basis_set:
            self.input_dict['basis_set_name'] = self.basis_set_name

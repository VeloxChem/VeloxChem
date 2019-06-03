import unittest
import hashlib
import os

from veloxchem.inputparser import InputParser


class TestBasis(unittest.TestCase):

    def test_basis_sets(self):

        basis_dir = 'basis'
        if not os.path.isdir(basis_dir):
            basis_dir = os.path.join(os.pardir, basis_dir)

        for basis_name in os.listdir(basis_dir):
            if basis_name == 'MIN-CC-PVDZ':
                continue

            if 'RI' in basis_name:
                continue

            basis_file = os.path.join(basis_dir, basis_name)
            if not os.path.isfile(basis_file):
                continue

            # get reference md5

            f_basis = open(basis_file, 'r')
            expected_md5 = f_basis.readlines()[-1].strip()
            f_basis.close()

            # get actual md5

            basis_parser = InputParser(basis_file)
            basis_dict = basis_parser.get_dict()

            basis_string = '@BASIS_SET ' + basis_name + os.linesep + os.linesep

            for key in basis_dict:
                if 'atombasis_' not in key:
                    continue

                content = key.split('_')
                basis_string += '@' + content[0].upper() + ' '
                basis_string += content[1].upper() + os.linesep

                for line in basis_dict[key]:
                    content = line.split()
                    if content[0] in 'SPDFGHIKL':
                        basis_string += content[0] + ' '
                        basis_string += content[1] + '  '
                        basis_string += content[2]
                    else:
                        basis_string += '{:>18s}'.format(content[0])
                        for s in content[1:]:
                            basis_string += ' {:>19s}'.format(s)
                    basis_string += os.linesep

                basis_string += '@END' + os.linesep + os.linesep

            calculated_md5 = hashlib.md5(
                basis_string[:-1].encode('utf-8')).hexdigest()

            # verify md5

            self.assertEqual(expected_md5, calculated_md5)


if __name__ == "__main__":
    unittest.main()

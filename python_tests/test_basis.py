import unittest
import hashlib
import os
from pathlib import Path

from veloxchem.inputparser import InputParser
from veloxchem.molecularbasis import MolecularBasis


class TestBasis(unittest.TestCase):

    def test_basis_sets(self):

        here = Path(__file__).parents[1]
        basis_dir = here / 'basis'

        for basis in basis_dir.iterdir():
            if not basis.is_file():
                continue

            if basis.name == 'MIN-CC-PVDZ':
                continue

            if 'RI' in basis.name:
                continue

            # get reference md5

            with basis.open('r') as f:
                expected_md5 = f.readlines()[-1].strip()

            # get actual md5

            basis_parser = InputParser(str(basis))
            basis_dict = basis_parser.get_dict()

            basis_string = '@BASIS_SET ' + basis.name + os.linesep + os.linesep

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

    def test_avail_basis(self):

        h_basis = [
            '6-311++G(2D,2P)', 'AUG-CC-PVTZ', 'DEF2-SVPD', 'MIN-CC-PVDZ',
            'STO-6G'
        ]
        c_basis = [
            '6-31G', '6-31G(2DF,P)', 'AUG-CC-PVDZ', 'D-AUG-CC-PVQZ',
            'SADLEJ-PVTZ'
        ]
        zn_basis = ['CC-PVTZ', 'MIN-CC-PVDZ', 'STO-3G']

        for bas in h_basis:
            self.assertTrue(bas in MolecularBasis.get_avail_basis('H'))
        for bas in c_basis:
            self.assertTrue(bas in MolecularBasis.get_avail_basis('c'))
        for bas in zn_basis:
            self.assertTrue(bas in MolecularBasis.get_avail_basis('Zn'))


if __name__ == "__main__":
    unittest.main()

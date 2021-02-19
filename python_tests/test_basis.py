import unittest
import hashlib
from pathlib import Path

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

            with basis.open('r') as f:
                basis_lines = f.readlines()

            # verify md5

            expected_md5 = basis_lines[-1].strip()

            calculated_md5 = hashlib.md5(''.join(
                basis_lines[:-1]).encode('utf-8')).hexdigest()

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

from pathlib import Path
import tempfile
import hashlib

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import AtomBasis
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestBasis:

    def build_basis_function(self, expons, coeffs, angl):

        bf = BasisFunction(expons, coeffs, angl)
        bf.normalize()

        return bf

    def test_basis_read(self):

        mol_text = """
            C        0.0000        0.0000       -1.2436
            C       -0.0000        0.0000        1.2436
            H        0.0000        1.7272       -2.3230
            H        0.0000       -1.7272       -2.3230
            H       -0.0000        1.7272        2.3230
            H       -0.0000       -1.7272        2.3230
        """

        basis_text = """@BASIS_SET TESTBASIS
            @ATOMBASIS H
            S 2  1
            3.425250910000e+00  1.543289700000e-01
            6.239137300000e-01  5.353281400000e-01
            @END
            @ATOMBASIS C
            S 1  1
            7.161683700000e+01  1.543289700000e-01
            P 1  1
            2.941249400000e+00  1.559162700000e-01
            @END
        """

        basis_name = 'TESTBASIS'

        if is_mpi_master():

            molecule = Molecule.read_str(mol_text, units='angstrom')

            with tempfile.TemporaryDirectory() as temp_dir:
                fname = Path(temp_dir, basis_name)
                with open(str(fname), 'w') as f_basis:
                    f_basis.write(basis_text)
                basis = MolecularBasis.read(molecule,
                                            basis_name,
                                            temp_dir,
                                            ostream=None)

            ref_basis = MolecularBasis()

            atom_basis = AtomBasis()
            atom_basis.add_basis_function(
                self.build_basis_function(
                    [3.425250910000e+00, 6.239137300000e-01],
                    [1.543289700000e-01, 5.353281400000e-01], 0))
            atom_basis.set_elemental_id(1)
            ref_basis.add_atom_basis(atom_basis)

            atom_basis = AtomBasis()
            atom_basis.add_basis_function(
                self.build_basis_function([7.161683700000e+01],
                                          [1.543289700000e-01], 0))
            atom_basis.add_basis_function(
                self.build_basis_function([2.941249400000e+00],
                                          [1.559162700000e-01], 1))
            atom_basis.set_elemental_id(6)
            ref_basis.add_atom_basis(atom_basis)

            ref_basis.set_label(basis_name)

            assert basis == ref_basis

    def test_basis_md5(self):

        here = Path(__file__).parents[1]
        basis_dir = here / 'basis'

        for basis in basis_dir.iterdir():
            if not is_mpi_master():
                continue

            if not basis.is_file():
                continue
            if 'RI' in basis.name:
                continue

            with basis.open('r') as f:
                basis_lines = f.readlines()

            # verify md5

            expected_md5 = basis_lines[-1].strip()

            calculated_md5 = hashlib.md5(''.join(
                basis_lines[:-1]).encode('utf-8')).hexdigest()

            assert expected_md5 == calculated_md5

    def test_avail_basis(self):

        h_basis = [
            '6-311++G(2D,2P)', 'AUG-CC-PVTZ', 'DEF2-SVPD', 'AO-START-GUESS',
            'STO-6G', 'DEF2-SV(P)', 'PCSEG-0'
        ]
        c_basis = [
            '6-31G', '6-31G(2DF,P)', 'AUG-CC-PVDZ', 'AUG-CC-PVQZ',
            'SADLEJ-PVTZ', 'AUG-CC-PCVTZ', 'PCX-2', 'ANO-S-MB'
        ]
        zn_basis = ['CC-PVTZ', 'AO-START-GUESS', 'STO-3G', 'ANO-L-VQZP']

        if is_mpi_master():

            h_avail_basis = MolecularBasis.get_avail_basis('H')
            c_avail_basis = MolecularBasis.get_avail_basis('C')
            zn_avail_basis = MolecularBasis.get_avail_basis('zn')

            for bas in h_basis:
                assert bas in h_avail_basis
            for bas in c_basis:
                assert bas in c_avail_basis
            for bas in zn_basis:
                assert bas in zn_avail_basis

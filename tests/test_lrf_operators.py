from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.cppsolver import ComplexResponse
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestLrfOperators:

    def run_lrf(self, basis_set_label, freq, operators, components, prefac,
                ref_val):

        molecule = Molecule.read_xyz_string("""1
            xyz
            Ne  0.0 0.0 0.0""")

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        lrf_drv = ComplexResponse()
        lrf_drv.ostream.mute()

        op_a, op_b = operators
        a, b = components

        rsp_settings = {
            'frequencies': [freq],
            'damping': 0.0,
            'a_operator': op_a,
            'b_operator': op_b,
            'a_components': a,
            'b_components': b,
        }
        lrf_drv.update_settings(rsp_settings)
        lrf_results = lrf_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            rsp_func = lrf_results['response_functions'][(a, b, freq)]
            calc_val = prefac * rsp_func
            assert abs(calc_val.real - ref_val) < 1.0e-3
            assert abs(calc_val.imag) < 1.0e-3

    def test_lrf_diplen(self):

        freq = 0.0656
        ref_val = 2.362

        operators = ['electric dipole', 'electric dipole']
        prefac = -1.0

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val)

    def test_lrf_dipvel(self):

        freq = 0.0656
        ref_val = 2.378

        operators = ['linear momentum', 'electric dipole']
        prefac = (-1j * (-1.0) / freq)

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val)

    def test_lrf_dipvel_flipped(self):

        freq = 0.0656
        ref_val = 2.378

        operators = ['electric dipole', 'linear momentum']
        prefac = -1.0 * (-1j * (-1.0) / freq)

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val)

    def test_lrf_quadrupole(self):

        freq = 0.0656

        operators = ['quadrupole', 'quadrupole']

        self.run_lrf('daug-cc-pvdz', freq, operators, ['xx', 'zz'], -1.0, 1.527)
        self.run_lrf('daug-cc-pvdz', freq, operators, ['xy', 'xy'], -1.0, 1.031)
        self.run_lrf('daug-cc-pvdz', freq, operators, ['zz', 'zz'], -1.0, 3.590)

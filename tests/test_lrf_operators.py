from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestLrfOperators:

    def run_lrf(self, basis_set_label, freq, operators, components, prefac,
                ref_val, ri):

        molecule = Molecule.read_xyz_string("""1
            xyz
            Ne  0.0 0.0 0.0""")

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lrf_drv = ComplexResponseSolver()
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

    @pytest.mark.parametrize('ri', [False, True])
    def test_lrf_diplen(self, ri):

        freq = 0.0656
        ref_val = {False: 2.362, True: 2.35906}[ri]

        operators = ['electric dipole', 'electric dipole']
        prefac = -1.0

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_lrf_dipvel(self, ri):

        freq = 0.0656
        ref_val = {False: 2.378, True: 2.37444}[ri]

        operators = ['linear momentum', 'electric dipole']
        prefac = (-1j * (-1.0) / freq)

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_lrf_dipvel_flipped(self, ri):

        freq = 0.0656
        ref_val = {False: 2.378, True: 2.37445}[ri]

        operators = ['electric dipole', 'linear momentum']
        prefac = -1.0 * (-1j * (-1.0) / freq)

        self.run_lrf('daug-cc-pvdz', freq, operators, 'zz', prefac, ref_val, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_lrf_quadrupole(self, ri):

        freq = 0.0656

        operators = ['quadrupole', 'quadrupole']

        self.run_lrf('daug-cc-pvdz', freq, operators, ['xx', 'zz'], -1.0, {
            False: 1.527,
            True: 1.52752
        }[ri], ri)
        self.run_lrf('daug-cc-pvdz', freq, operators, ['xy', 'xy'], -1.0, 1.031,
                     ri)
        self.run_lrf('daug-cc-pvdz', freq, operators, ['zz', 'zz'], -1.0, 3.590,
                     ri)

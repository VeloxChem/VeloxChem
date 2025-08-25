import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.doubleresbeta import DoubleResBetaDriver


@pytest.mark.solvers
class TestTransitionDipoleMomentQR:

    def get_molecule_and_basis(self):

        # this is not a good geometry for water
        # but is only used for testing purposes
        xyz_string = """3

        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        basis_label = 'sto-3g'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        return mol, bas

    def test_tsm_qr(self):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        tsm_drv = DoubleResBetaDriver()
        tsm_drv.ostream.mute()

        tol = 1.0e-5

        tsm_drv.initial_state = 5
        tsm_drv.final_state = 8
        tsm_results = tsm_drv.compute(mol, bas, scf_results)

        if tsm_drv.rank == mpi_master():
            calc_tsm = tsm_results['transition_dipole_moments']
            assert abs(abs(calc_tsm[('x', 5, 8)]) - 0.0) < tol
            assert abs(abs(calc_tsm[('y', 5, 8)]) - 0.0) < tol
            assert abs(abs(calc_tsm[('z', 5, 8)]) - 0.523176) < tol

        tsm_drv.initial_state = 7
        tsm_drv.final_state = 7
        tsm_results = tsm_drv.compute(mol, bas, scf_results)

        if tsm_drv.rank == mpi_master():
            calc_esm = tsm_results['excited_state_dipole_moments']
            assert abs(abs(calc_esm[('x', 7, 7)]) - 0.0) < tol
            assert abs(abs(calc_esm[('y', 7, 7)]) - 0.0) < tol
            assert abs(abs(calc_esm[('z', 7, 7)]) - 0.538214) < tol

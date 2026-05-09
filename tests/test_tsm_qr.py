import pytest
from mpi4py import MPI

import veloxchem.excitedstatemomentdriver as esm_driver_module
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.excitedstatemomentdriver import ExcitedStateMomentDriver
from veloxchem.errorhandler import VeloxChemError


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

    def get_scf_results(self):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        return mol, bas, scf_results

    def test_tsm_qr(self):

        mol, bas, scf_results = self.get_scf_results()

        tsm_drv = ExcitedStateMomentDriver()
        tsm_drv.ostream.mute()

        tol = 1.0e-5

        # tsm_drv._initial_state = 5
        # tsm_drv._final_state = 8
        # tsm_results = tsm_drv.compute(mol, bas, scf_results)

        # if tsm_drv.rank == mpi_master():
        #     calc_tsm = tsm_results['transition_dipole_moment']
        #     assert abs(calc_tsm[0]) < tol
        #     assert abs(calc_tsm[1]) < tol
        #     assert abs(abs(calc_tsm[2]) - 0.523176) < tol

        tsm_drv.state = 7
        tsm_results = tsm_drv.compute(mol, bas, scf_results)

        if tsm_drv.rank == mpi_master():
            calc_esm = tsm_results['excited_state_dipole_moment']
            assert abs(calc_esm[0]) < tol
            assert abs(calc_esm[1]) < tol
            assert abs(calc_esm[2] - 0.053558) < tol

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_tsm_qr_rejects_nonpositive_state_indices(self):

        mol, bas, scf_results = self.get_scf_results()

        tsm_drv = ExcitedStateMomentDriver()
        tsm_drv.ostream.mute()
        tsm_drv.state = 0

        with pytest.raises(VeloxChemError,
                           match='Expecting positive 1-based state index'):
            tsm_drv.compute(mol, bas, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_tsm_qr_rejects_unconverged_complex_response(self, monkeypatch):

        mol, bas, scf_results = self.get_scf_results()

        original_compute = esm_driver_module.ComplexResponseSolver.compute

        def wrapped_compute(self, *args, **kwargs):
            result = original_compute(self, *args, **kwargs)
            self._is_converged = False
            return result

        monkeypatch.setattr(esm_driver_module.ComplexResponseSolver, 'compute',
                            wrapped_compute)

        tsm_drv = ExcitedStateMomentDriver()
        tsm_drv.ostream.mute()
        tsm_drv.state = 1

        with pytest.raises(VeloxChemError,
                           match='Complex response solver did not converge'):
            tsm_drv.compute(mol, bas, scf_results)

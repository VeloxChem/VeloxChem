from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver

import pytest


@pytest.mark.solvers
class TestEFieldRestart:

    @staticmethod
    def get_molecule_and_basis():

        xyz_string = """
        3
        xyz
        O    0.0000000   0.0000000   0.1173000
        H    0.0000000   0.7572000  -0.4692000
        H    0.0000000  -0.7572000  -0.4692000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        return mol, bas

    @staticmethod
    def get_checkpoint_filename(scf_drv, tmp_path, name):

        filename = str(tmp_path / name)
        # to avoid inconsistency across MPI ranks
        return scf_drv.comm.bcast(filename, root=mpi_master())

    @staticmethod
    def run_scf(electric_field, filename):

        mol, bas = TestEFieldRestart.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.electric_field = electric_field
        scf_drv.filename = filename

        scf_results = scf_drv.compute(mol, bas)

        return scf_drv, scf_results

    def test_scf_restart_with_matching_field(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.electric_field = (0.0, 0.0, 0.001)
        scf_drv.filename = self.get_checkpoint_filename(scf_drv, tmp_path,
                                                        'water_matching')

        scf_results_ref = scf_drv.compute(mol, bas)

        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        # a restart is really a restart: the driver accepts the checkpoint
        # on every rank
        assert scf_drv.restart is True

        if scf_drv.rank == mpi_master():
            assert abs(scf_results['scf_energy'] -
                       scf_results_ref['scf_energy']) < 1.0e-6

    def test_scf_restart_without_field_from_field_checkpoint(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.electric_field = (0.0, 0.0, 0.001)
        scf_drv.filename = self.get_checkpoint_filename(scf_drv, tmp_path,
                                                        'water_dropped_field')

        scf_results_ref = scf_drv.compute(mol, bas)

        scf_drv.restart = True
        scf_drv.electric_field = None
        scf_results = scf_drv.compute(mol, bas)

        # SCF restart is permissive by design: the checkpoint orbitals are
        # only an initial guess, so dropping the field is allowed, and the
        # restart is really a restart (the driver accepts the checkpoint)
        assert scf_drv.restart is True

        _, fresh_results = self.run_scf(None, None)

        if scf_drv.rank == mpi_master():
            # the restarted calculation converges to the regular zero-field
            # SCF result ...
            assert abs(scf_results['scf_energy'] -
                       fresh_results['scf_energy']) < 1.0e-6
            # ... and not to the checkpoint's field-perturbed state: the
            # field effect must exceed the equality tolerance, otherwise the
            # test would not discriminate
            assert abs(scf_results['scf_energy'] -
                       scf_results_ref['scf_energy']) > 1.0e-5

    def test_scf_restart_with_field_from_fieldfree_checkpoint(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.filename = self.get_checkpoint_filename(scf_drv, tmp_path,
                                                        'water_added_field')

        scf_results_ref = scf_drv.compute(mol, bas)

        scf_drv.restart = True
        scf_drv.electric_field = (0.0, 0.0, 0.001)
        scf_results = scf_drv.compute(mol, bas)

        # adding a field on restart is likewise allowed: the restart is
        # really a restart (the driver accepts the checkpoint)
        assert scf_drv.restart is True

        _, fresh_results = self.run_scf((0.0, 0.0, 0.001), None)

        if scf_drv.rank == mpi_master():
            # converges to the field-perturbed result, not to the checkpoint's
            # zero-field state (field effect must exceed the tolerance)
            assert abs(scf_results['scf_energy'] -
                       fresh_results['scf_energy']) < 1.0e-6
            assert abs(scf_results['scf_energy'] -
                       scf_results_ref['scf_energy']) > 1.0e-5

    def test_scf_restart_with_changed_field(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.electric_field = (0.0, 0.0, 0.001)
        scf_drv.filename = self.get_checkpoint_filename(scf_drv, tmp_path,
                                                        'water_changed_field')

        scf_results_ref = scf_drv.compute(mol, bas)

        scf_drv.restart = True
        scf_drv.electric_field = (0.0, 0.0, -0.001)
        scf_results = scf_drv.compute(mol, bas)

        # the restart is really a restart: the driver accepts the checkpoint
        assert scf_drv.restart is True

        _, fresh_results = self.run_scf((0.0, 0.0, -0.001), None)

        if scf_drv.rank == mpi_master():
            # converges to the new field's result, not to the checkpoint's
            # old-field state (field effect must exceed the tolerance)
            assert abs(scf_results['scf_energy'] -
                       fresh_results['scf_energy']) < 1.0e-6
            assert abs(scf_results['scf_energy'] -
                       scf_results_ref['scf_energy']) > 1.0e-5

    def test_scf_restart_without_field_from_fieldfree_checkpoint(self,
                                                                 tmp_path):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.filename = self.get_checkpoint_filename(scf_drv, tmp_path,
                                                        'water_fieldfree')

        scf_results_ref = scf_drv.compute(mol, bas)

        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        # a restart is really a restart: the driver accepts the checkpoint
        # on every rank
        assert scf_drv.restart is True

        if scf_drv.rank == mpi_master():
            assert abs(scf_results['scf_energy'] -
                       scf_results_ref['scf_energy']) < 1.0e-6

    def test_response_match_settings_rejects_changed_field(self, tmp_path):

        rsp_drv = LinearResponseEigenSolver()
        rsp_drv.ostream.mute()
        rsp_drv.checkpoint_file = rsp_drv.comm.bcast(
            str(tmp_path / 'water_efield_rsp.h5'), root=mpi_master())

        rsp_drv.electric_field = (0.0, 0.0, 0.001)
        rsp_drv._write_settings_to_checkpoint()

        if rsp_drv.rank == mpi_master():
            # response restart requires an exact settings match
            assert rsp_drv.match_settings(rsp_drv.checkpoint_file) is True

            rsp_drv.electric_field = (0.0, 0.0, -0.001)
            assert rsp_drv.match_settings(rsp_drv.checkpoint_file) is False

            rsp_drv.electric_field = None
            assert rsp_drv.match_settings(rsp_drv.checkpoint_file) is False

            rsp_drv.electric_field = (0.0, 0.0, 0.001)
            assert rsp_drv.match_settings(rsp_drv.checkpoint_file) is True

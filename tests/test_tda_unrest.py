from pathlib import Path
import h5py
from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaeigensolverunrest import TdaUnrestrictedEigenSolver
from veloxchem.inputparser import (unparse_input, read_unparsed_input_from_hdf5)


@pytest.mark.solvers
class TestUnrestrictedTDA:

    def get_water_cation_system(self, basis_label='def2-svp'):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1.0)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        return mol, bas

    def run_unrestricted_scf(self, mol, bas, xcfun_label='hf'):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label

        return scf_drv.compute(mol, bas)

    def run_restarted_unrestricted_scf(self,
                                       mol,
                                       bas,
                                       filename,
                                       xcfun_label='hf'):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.filename = filename
        scf_drv.xcfun = xcfun_label
        scf_drv.compute(mol, bas)
        scf_drv.restart = True

        return scf_drv.compute(mol, bas)

    def configure_tda(self, driver, filename):

        driver.filename = filename
        driver.conv_thresh = 1.0e-5
        driver.max_iter = 60
        driver.nstates = 2
        driver.initial_guess_multiplier = 2
        driver.max_subspace_dim = 8
        driver.non_equilibrium_solv = False
        driver.ri_auxiliary_basis = 'def2-universal-jfit'

    def normalize_setting(self, value):

        if isinstance(value, (list, tuple)):
            return tuple(value)

        return value

    def run_tda(self,
                xcfun_label,
                basis_label,
                ref_exc_enes,
                ref_osc_str,
                tol,
                max_subspace_dim=None):

        xyz_string = """3
        xyz
        O      0.000000   0.000000   0.117790
        H      0.000000   0.755453  -0.471161
        H      0.000000  -0.755453  -0.471161
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 10
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_camb3lyp_svp(self):

        # vlxtag: UKS, Absorption, TDA

        ref_exc_enes = np.array([
            0.09065604, 0.09500459, 0.21881146, 0.47239371, 0.48842790,
            0.51728556, 0.52413401, 0.55026290, 0.55530527, 0.59463209
        ])

        ref_osc_str = np.array([
            0.0013, 0.2065, 0.0000, 0.0325, 0.1953, 0.0015, 0.0000, 0.0256,
            0.0701, 0.1727
        ])

        self.run_tda('cam-b3lyp',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     max_subspace_dim=80)

    def run_tda_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_with_ecp(self):

        # vlxtag: UHF, Absorption, TDA

        ref_exc_enes = np.array(
            [0.08351815, 0.08351815, 0.10137196, 0.10137196, 0.15068886])
        ref_osc_str = np.array([0.0000, 0.0000, 0.0002, 0.0002, 0.0089])

        self.run_tda_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_checkpoint_and_restart(self, tmp_path):

        mol, bas = self.get_water_cation_system()

        filename = str(tmp_path / "water_restart")
        filename = MPI.COMM_WORLD.bcast(filename, root=mpi_master())

        scf_results = self.run_restarted_unrestricted_scf(mol, bas, filename)

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.nstates = 5
        lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.nstates = 10
        lr_results_first = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is True

        lr_drv.restart = False
        lr_drv.nstates = 10
        lr_results_second = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is False

        if lr_drv.rank == mpi_master():
            keys = ['eigenvalues', 'oscillator_strengths']
            tols = [1e-8, 1e-5]
            for key, tol in zip(keys, tols):
                assert np.max(
                    np.abs(lr_results_first[key] -
                           lr_results_second[key])) < tol

    def test_response_settings_roundtrip(self, tmp_path):

        mol, bas = self.get_water_cation_system('sto-3g')
        scf_results = self.run_unrestricted_scf(mol, bas)

        driver = TdaUnrestrictedEigenSolver()
        driver.ostream.mute()

        filename = str(tmp_path / 'tdaunrestrictedeigensolver')
        filename = driver.comm.bcast(filename, root=mpi_master())

        self.configure_tda(driver, filename)

        reference_scf_results = dict(scf_results)
        reference_scf_results['filename'] = filename

        reference_results = driver.compute(mol, bas, reference_scf_results)
        checkpoint_file = Path(driver.comm.bcast(driver.checkpoint_file,
                                                 root=mpi_master()))

        rsp_keywords = {
            key: val[0] for key, val in driver._input_keywords['response'].items()
        }
        method_keywords = {
            key: val[0]
            for key, val in driver._input_keywords['method_settings'].items()
        }

        expected_rsp = unparse_input(driver, rsp_keywords)
        expected_method = unparse_input(driver, method_keywords)

        if driver.rank == mpi_master():
            assert checkpoint_file.is_file()
            checkpoint_rsp_input = read_unparsed_input_from_hdf5(
                str(checkpoint_file), group_name='response_settings')
            checkpoint_method_input = read_unparsed_input_from_hdf5(
                str(checkpoint_file), group_name='method_settings')
        else:
            checkpoint_rsp_input = None
            checkpoint_method_input = None

        checkpoint_rsp_input = driver.comm.bcast(checkpoint_rsp_input,
                                                 root=mpi_master())
        checkpoint_method_input = driver.comm.bcast(checkpoint_method_input,
                                                    root=mpi_master())

        assert checkpoint_rsp_input == expected_rsp
        assert checkpoint_method_input == expected_method

        second_drv = TdaUnrestrictedEigenSolver()
        second_drv.ostream.mute()
        second_drv.update_settings(checkpoint_rsp_input, checkpoint_method_input)

        for key in rsp_keywords:
            assert self.normalize_setting(getattr(second_drv,
                                                 key)) == self.normalize_setting(
                                                     getattr(driver, key))
        for key in method_keywords:
            assert self.normalize_setting(getattr(second_drv,
                                                 key)) == self.normalize_setting(
                                                     getattr(driver, key))

        third_drv = TdaUnrestrictedEigenSolver()
        third_drv.ostream.mute()
        third_drv.restart = False
        copied_filename = driver.comm.bcast(
            str(tmp_path / 'tdaunrestrictedeigensolver_copy'), root=mpi_master())
        copied_checkpoint_file = driver.comm.bcast(
            str(tmp_path / 'tdaunrestrictedeigensolver_copy.h5'),
            root=mpi_master())
        third_drv.filename = copied_filename
        third_drv.checkpoint_file = copied_checkpoint_file

        third_drv.read_settings(str(checkpoint_file))

        assert third_drv.restart is False
        assert third_drv.filename == copied_filename
        assert third_drv.checkpoint_file == copied_checkpoint_file

        copied_scf_results = dict(scf_results)
        copied_scf_results['filename'] = third_drv.filename

        copied_results = third_drv.compute(mol, bas, copied_scf_results)

        for key in rsp_keywords:
            if key in ('restart', 'filename', 'checkpoint_file'):
                continue
            assert self.normalize_setting(getattr(third_drv,
                                                  key)) == self.normalize_setting(
                                                      getattr(driver, key))
        for key in method_keywords:
            assert self.normalize_setting(getattr(third_drv,
                                                  key)) == self.normalize_setting(
                                                      getattr(driver, key))

        if third_drv.rank == mpi_master():
            assert np.max(
                np.abs(copied_results['eigenvalues'] -
                       reference_results['eigenvalues'])) < 1.0e-8

    def test_restart_allows_missing_nstates_dataset(self, tmp_path):

        mol, bas = self.get_water_cation_system('sto-3g')

        filename = str(tmp_path / 'tda_unrestricted_missing_nstates')
        driver = TdaUnrestrictedEigenSolver()
        driver.ostream.mute()
        filename = driver.comm.bcast(filename, root=mpi_master())

        scf_results = self.run_restarted_unrestricted_scf(mol, bas, filename)
        reference_scf_results = dict(scf_results)
        reference_scf_results['filename'] = filename

        self.configure_tda(driver, filename)
        reference_results = driver.compute(mol, bas, reference_scf_results)
        checkpoint_file = driver.comm.bcast(driver.checkpoint_file,
                                            root=mpi_master())

        if driver.rank == mpi_master():
            with h5py.File(checkpoint_file, 'a') as h5f:
                del h5f['nstates']

        driver.comm.barrier()

        restarted_drv = TdaUnrestrictedEigenSolver()
        restarted_drv.ostream.mute()
        self.configure_tda(restarted_drv, filename)
        restarted_drv.checkpoint_file = checkpoint_file
        restarted_drv.restart = True

        restarted_results = restarted_drv.compute(mol, bas, reference_scf_results)

        assert restarted_drv.restart is True

        if restarted_drv.rank == mpi_master():
            assert np.max(
                np.abs(restarted_results['eigenvalues'] -
                       reference_results['eigenvalues'])) < 1.0e-8

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_missing_num_core_orbitals(self):

        mol, bas = self.get_water_cation_system()
        scf_results = self.run_unrestricted_scf(mol, bas, 'cam-b3lyp')

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.core_excitation = True
        lr_drv.nstates = 2

        with pytest.raises(AssertionError,
                           match='num_core_orbitals not set or invalid'):
            lr_drv.compute(mol, bas, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_too_large_num_core_orbitals(self):

        mol, bas = self.get_water_cation_system()
        scf_results = self.run_unrestricted_scf(mol, bas, 'cam-b3lyp')

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 5
        lr_drv.nstates = 2

        with pytest.raises(AssertionError,
                           match='num_core_orbitals too large'):
            lr_drv.compute(mol, bas, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_restricted_subspace(self):

        mol, bas = self.get_water_cation_system()
        scf_results = self.run_unrestricted_scf(mol, bas)

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.restricted_subspace = True
        lr_drv.nstates = 2

        with pytest.raises(AssertionError,
                           match='restricted_subspace not implemented'):
            lr_drv.compute(mol, bas, scf_results)

from copy import deepcopy
from pathlib import Path
import h5py
from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.inputparser import (unparse_input, read_unparsed_input_from_hdf5)
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestTDA:

    def get_water_system(self, basis_label='def2-svp'):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        return mol, bas

    def run_restricted_scf(self, mol, bas, xcfun_label='hf'):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label

        return scf_drv.compute(mol, bas)

    def run_restarted_restricted_scf(self, mol, bas, filename, xcfun_label='hf'):

        scf_drv = ScfRestrictedDriver()
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
                ri_coulomb=False,
                max_subspace_dim=None):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = ri_coulomb
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_svp(self):

        # vlxtag: RHF, Absorption, CIS

        ref_exc_enes = np.array(
            [0.34189725, 0.40717708, 0.43577895, 0.50150019, 0.55485041])

        ref_osc_str = np.array([0.0229, 0.0000, 0.1039, 0.0975, 0.3060])

        self.run_tda('hf',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-6,
                     max_subspace_dim=30)

    def test_b3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28045483, 0.35053288, 0.36504833, 0.43904962, 0.51570229])

        ref_osc_str = np.array([0.0180, 0.0000, 0.0854, 0.0695, 0.3006])

        self.run_tda('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28370223, 0.35557650, 0.36883370, 0.44496286, 0.51683587])

        ref_osc_str = np.array(
            [0.017966, 0.000000, 0.084011, 0.068027, 0.300331])

        self.run_tda('cam-b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_tzvp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28085567, 0.35172399, 0.36474654, 0.43912355, 0.50442105])

        ref_osc_str = np.array(
            [0.034572, 0.000000, 0.107302, 0.059653, 0.249276])

        self.run_tda('cam-b3lyp', 'def2-tzvp', ref_exc_enes, ref_osc_str,
                     1.0e-5)

    def test_ri_blyp_svp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.25974113, 0.33054867, 0.34244014, 0.41789336, 0.50103038])

        ref_osc_str = np.array(
            [0.016643, 0.000000, 0.081274, 0.066394, 0.293944])

        self.run_tda('blyp',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     ri_coulomb=True)

    def run_tda_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """2
        xyz
        Au 0 0 0
        Au 0 0 2.88
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_with_ecp(self):

        # vlxtag: RHF, Absorption, CIS

        ref_exc_enes = np.array([0.11555805, 0.15539632, 0.16343183])
        ref_osc_str = np.array([0.6238, 0.0000, 0.0027])

        self.run_tda_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_checkpoint_and_restart(self, tmp_path):

        mol, bas = self.get_water_system()

        filename = str(tmp_path / "water_restart")
        filename = MPI.COMM_WORLD.bcast(filename, root=mpi_master())

        scf_results = self.run_restarted_restricted_scf(mol, bas, filename)

        lr_drv = TdaEigenSolver()
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

        mol, bas = self.get_water_system('sto-3g')
        scf_results = self.run_restricted_scf(mol, bas)

        driver = TdaEigenSolver()
        driver.ostream.mute()

        filename = str(tmp_path / 'tdaeigensolver')
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

        second_drv = TdaEigenSolver()
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

        third_drv = TdaEigenSolver()
        third_drv.ostream.mute()
        third_drv.restart = False
        copied_filename = driver.comm.bcast(str(tmp_path / 'tdaeigensolver_copy'),
                                            root=mpi_master())
        copied_checkpoint_file = driver.comm.bcast(
            str(tmp_path / 'tdaeigensolver_copy.h5'), root=mpi_master())
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

    def test_update_settings_enables_detach_attach_cubes(self):

        driver = TdaEigenSolver()
        driver.update_settings({'detach_attach_cubes': 'yes'})

        assert driver.detach_attach_cubes is True
        assert driver.detach_attach is True

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize(
        'rsp_dict, message',
        [
            ({'cube_origin': '0.0, 1.0'}, 'cube origin needs 3 numbers'),
            ({'cube_stepsize': '0.1, 0.2'}, 'cube stepsize needs 3 numbers'),
            ({'cube_points': '10, 12'}, 'cube points needs 3 integers'),
        ],
    )
    def test_update_settings_rejects_invalid_cube_settings(self, rsp_dict,
                                                           message):

        driver = TdaEigenSolver()

        with pytest.raises(VeloxChemError, match=message):
            driver.update_settings(rsp_dict)

    def test_driver_deepcopy(self):

        driver = TdaEigenSolver()
        driver.nstates = 4
        driver.core_excitation = True
        driver.num_core_orbitals = 2
        driver.cube_points = [12, 14, 16]

        copied = deepcopy(driver)

        assert isinstance(copied, TdaEigenSolver)
        assert copied is not driver
        assert copied.comm == driver.comm
        assert copied.ostream == driver.ostream
        assert copied.nstates == driver.nstates
        assert copied.core_excitation is True
        assert copied.num_core_orbitals == driver.num_core_orbitals
        assert copied.cube_points == driver.cube_points
        assert copied.cube_points is not driver.cube_points

    def test_restart_allows_missing_nstates_dataset(self, tmp_path):

        mol, bas = self.get_water_system('sto-3g')

        filename = str(tmp_path / 'tda_missing_nstates')
        driver = TdaEigenSolver()
        driver.ostream.mute()
        filename = driver.comm.bcast(filename, root=mpi_master())

        scf_results = self.run_restarted_restricted_scf(mol, bas, filename)
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

        restarted_drv = TdaEigenSolver()
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
                        reason='skip plot for multiple MPI processes')
    def test_plot_valence_spectra(self, monkeypatch):

        plt = pytest.importorskip('matplotlib.pyplot')

        mol, bas = self.get_water_system('sto-3g')
        scf_results = self.run_restricted_scf(mol, bas)

        lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_results = lr_drv.compute(mol, bas, scf_results)

        monkeypatch.setattr(plt, 'show', lambda: None)

        try:
            lr_drv.plot_uv_vis(lr_results,
                               broadening_type='gaussian',
                               broadening_value=0.10)
            lr_drv.plot_ecd(lr_results,
                            broadening_type='lorentzian',
                            broadening_value=0.15)
            lr_drv.plot(lr_results, plot_type='uv')
            lr_drv.plot(lr_results, plot_type='ecd')
            lr_drv.plot(lr_results, plot_type='electronic')
        finally:
            plt.close('all')

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip plot for multiple MPI processes')
    def test_plot_core_spectra(self, monkeypatch):

        plt = pytest.importorskip('matplotlib.pyplot')

        mol, bas = self.get_water_system('def2-svpd')
        scf_results = self.run_restricted_scf(mol, bas, 'b3lyp')

        lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.xcfun = 'b3lyp'
        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 1
        lr_drv.nstates = 3
        lr_results = lr_drv.compute(mol, bas, scf_results)

        monkeypatch.setattr(plt, 'show', lambda: None)

        try:
            lr_drv.plot_xas(lr_results,
                            broadening_type='gaussian',
                            broadening_value=0.10)
            lr_drv.plot_xcd(lr_results,
                            broadening_type='lorentzian',
                            broadening_value=0.15)
            lr_drv.plot(lr_results, plot_type='xas')
            lr_drv.plot(lr_results, plot_type='xcd')
            lr_drv.plot(lr_results, plot_type='electronic')
        finally:
            plt.close('all')

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_plot_rejects_restricted_subspace(self):

        lr_drv = TdaEigenSolver()
        lr_drv.restricted_subspace = True

        with pytest.raises(
                VeloxChemError,
                match='Plotting spectrum for restricted_subspace is not implemented.'
        ):
            lr_drv.plot_uv_vis({})

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_plot_rejects_invalid_modes(self):

        pytest.importorskip('matplotlib.pyplot')

        lr_drv = TdaEigenSolver()

        with pytest.raises(VeloxChemError, match='Invalid plot type'):
            lr_drv.plot({}, plot_type='invalid')

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_openshell_molecule(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match="TdaEigenSolver: not implemented for unrestricted case"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

import numpy as np
from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.distributedarray import DistributedArray


@pytest.mark.solvers
class TestRPABase:

    def get_water_system(self, basis_label='sto-3g'):

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

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_update_settings_and_full_solution_vector(self):

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()

        lr_drv.detach_attach_cubes = True
        lr_drv.detach_attach_charges = True
        lr_drv.update_settings({})
        assert lr_drv.detach_attach is True

        lr_drv.cube_origin = [0.0, 0.0]
        with pytest.raises(AssertionError, match='cube origin needs 3 numbers'):
            lr_drv.update_settings({})
        lr_drv.cube_origin = None

        lr_drv.cube_stepsize = [0.2, 0.2]
        with pytest.raises(AssertionError,
                           match='cube stepsize needs 3 numbers'):
            lr_drv.update_settings({})
        lr_drv.cube_stepsize = None

        lr_drv.cube_points = [80, 80]
        with pytest.raises(AssertionError,
                           match='cube points needs 3 integers'):
            lr_drv.update_settings({})

        solution = DistributedArray(np.array([[1.0, 2.0], [3.0, 4.0]]),
                                    lr_drv.comm,
                                    distribute=False)
        full_vector = lr_drv.get_full_solution_vector(solution)
        assert np.allclose(full_vector, np.array([3.0, 7.0, -1.0, -1.0]))

    def test_spectrum_helpers_valence(self):

        mol, bas = self.get_water_system('sto-3g')
        scf_results = self.run_restricted_scf(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            spectrum_au = lr_drv.get_absorption_spectrum(
                lr_results, [0.3, 0.4], 'au', 0.05, 'au')
            spectrum_ev = lr_drv.get_absorption_spectrum(
                lr_results, [10.0, 11.0], 'ev', 0.10, 'ev')
            ecd_nm = lr_drv.get_ecd_spectrum(lr_results, [200.0, 250.0], 'nm',
                                             0.20, 'ev')

            assert spectrum_au['x_label'] == 'Photon energy [a.u.]'
            assert spectrum_ev['x_label'] == 'Photon energy [eV]'
            assert ecd_nm['x_label'] == 'Wavelength [nm]'
            assert len(spectrum_au['y_data']) == 2
            assert len(ecd_nm['y_data']) == 2

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip plot for multiple MPI processes')
    def test_plot_valence_spectra(self, monkeypatch):

        plt = pytest.importorskip('matplotlib.pyplot')

        mol, bas = self.get_water_system('sto-3g')
        scf_results = self.run_restricted_scf(mol, bas)

        lr_drv = LinearResponseEigenSolver()
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

        lr_drv = LinearResponseEigenSolver()
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
    def test_plot_rejects_invalid_modes(self):

        lr_drv = LinearResponseEigenSolver()
        lr_drv.restricted_subspace = True

        with pytest.raises(AssertionError,
                           match='Plotting spectrum for restricted_subspace ' +
                           'is not implemented.'):
            lr_drv.plot_uv_vis({})

        lr_drv.restricted_subspace = False

        with pytest.raises(AssertionError, match='Invalid plot type'):
            lr_drv.plot({}, plot_type='invalid')

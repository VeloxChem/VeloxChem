from mpi4py import MPI
from pathlib import Path
import h5py
import numpy as np
import pytest

from veloxchem.veloxchemlib import (mpi_master, hartree_in_ev,
                                    hartree_in_inverse_nm)
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.resultsio import read_results
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tpatransitiondriver import TpaTransitionDriver


@pytest.mark.solvers
class TestTpaTransition:

    def run_scf(self, xcfun_label, ri_coulomb=False):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """
        basis_set_label = 'def2-svpd'
        scf_conv_thresh = 1.0e-8

        molecule = Molecule.read_molecule_string(molecule_string)
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = ri_coulomb
        scf_drv.conv_thresh = scf_conv_thresh
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        return scf_results, molecule, basis

    def run_tpatransition(self, xcfun_label, ref_energies, ref_strengths,
                          ri_coulomb=False):

        tpa_nstates = 2
        tpa_conv_thresh = 1.0e-7

        scf_results, molecule, ao_basis = self.run_scf(xcfun_label,
                                                       ri_coulomb=ri_coulomb)

        tpa_drv = TpaTransitionDriver()
        tpa_drv.xcfun = xcfun_label
        tpa_drv.nstates = tpa_nstates
        tpa_drv.conv_thresh = tpa_conv_thresh
        tpa_drv.ostream.mute()
        tpa_results = tpa_drv.compute(molecule, ao_basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            tpa_str = tpa_results['tpa_strengths']['linear']
            assert len(tpa_str) == len(ref_strengths)

            for freq, ref_freq in zip(tpa_results['photon_energies'],
                                      ref_energies):
                assert abs(freq / ref_freq - 1.0) < 1.0e-6

            for val, ref_val in zip(tpa_str, ref_strengths):
                assert abs(val / ref_val - 1.0) < 1.0e-6

    def compute_tpatransition(self, xcfun_label, ri_coulomb=False):

        scf_results, molecule, ao_basis = self.run_scf(xcfun_label,
                                                       ri_coulomb=ri_coulomb)

        tpa_drv = TpaTransitionDriver()
        tpa_drv.xcfun = xcfun_label
        tpa_drv.nstates = 2
        tpa_drv.conv_thresh = 1.0e-7
        tpa_drv.ostream.mute()

        return tpa_drv.compute(molecule, ao_basis, scf_results)

    def run_scf_restart(self, xcfun_label, filename, ri_coulomb=False):

        scf_results, molecule, basis = self.run_scf(xcfun_label,
                                                    ri_coulomb=ri_coulomb)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = ri_coulomb
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.filename = filename
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)
        scf_drv.restart = True
        restarted_results = scf_drv.compute(molecule, basis)

        assert scf_drv.restart

        return restarted_results, molecule, basis

    def test_tpatransition_hf(self):

        ref_energies = [0.1656922537003149, 0.21190963394768278]
        ref_strengths = [4.130284073574981, 15.83816511497292]
        self.run_tpatransition('hf', ref_energies, ref_strengths)

    def test_tpatransition_lda(self):

        ref_energies = [0.13269116242012727, 0.18051942347838879]
        ref_strengths = [12.922686324743287, 55.35603386346646]
        self.run_tpatransition('slda', ref_energies, ref_strengths)

    def test_tpatransition_gga(self):

        ref_energies = [0.13304132289794054, 0.179389130413857]
        ref_strengths = [13.068064328595725, 57.062788455439296]
        self.run_tpatransition('bp86', ref_energies, ref_strengths)

    def test_tpatransition_mgga(self):

        ref_energies = [0.1381414072229231, 0.18265658100098003]
        ref_strengths = [9.579957827267322, 41.94340887233736]
        self.run_tpatransition('tpssh', ref_energies, ref_strengths)

    def test_tpatransition_ri_blyp(self):

        ref_energies = [0.12580806392701718, 0.1710157238983142]
        ref_strengths = [13.672823155898364, 60.526749482300374]
        self.run_tpatransition('blyp',
                               ref_energies,
                               ref_strengths,
                               ri_coulomb=True)

    def test_get_spectrum(self):

        tpa_results = self.compute_tpatransition('hf')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            x_data_au = [0.16, 0.17, 0.18]
            x_data_ev = [x * hartree_in_ev() for x in x_data_au]
            x_data_nm = [1.0 / (hartree_in_inverse_nm() * x) for x in x_data_au]

            spectrum_au = TpaTransitionDriver.get_spectrum(tpa_results,
                                                           x_data_au, 'au',
                                                           0.01, 'au')
            spectrum_ev = TpaTransitionDriver.get_spectrum(
                tpa_results, x_data_ev, 'ev', 0.01 * hartree_in_ev(), 'ev')
            spectrum_nm = TpaTransitionDriver.get_spectrum(tpa_results,
                                                           x_data_nm, 'nm',
                                                           0.01, 'au')

            assert spectrum_au['x_label'] == 'Photon energy [a.u.]'
            assert spectrum_ev['x_label'] == 'Photon energy [eV]'
            assert spectrum_nm['x_label'] == 'Wavelength [nm]'
            assert spectrum_au['y_label'] == 'TPA cross-section [GM]'

            assert spectrum_au['x_data'] == pytest.approx(x_data_au)
            assert spectrum_ev['x_data'] == pytest.approx(x_data_ev)
            assert spectrum_nm['x_data'] == pytest.approx(x_data_nm)

            assert spectrum_au['y_data'] == pytest.approx(spectrum_ev['y_data'],
                                                          rel=1.0e-12,
                                                          abs=1.0e-12)
            assert spectrum_au['y_data'] == pytest.approx(spectrum_nm['y_data'],
                                                          rel=1.0e-12,
                                                          abs=1.0e-12)
            assert all(val >= 0.0 for val in spectrum_au['y_data'])

    def test_get_spectrum_polarization(self):

        tpa_results = {
            'photon_energies': [0.16, 0.18],
            'tpa_strengths': {
                'linear': [1.0, 2.0],
                'circular': [10.0, 20.0],
            },
        }

        x_data = [0.17]
        spectrum_linear = TpaTransitionDriver.get_spectrum(
            tpa_results, x_data, 'au', 0.01, 'au', polarization='linear')
        spectrum_circular = TpaTransitionDriver.get_spectrum(
            tpa_results, x_data, 'au', 0.01, 'au', polarization='circular')

        assert spectrum_circular['y_data'] == pytest.approx(
            [10.0 * spectrum_linear['y_data'][0]])

    @pytest.mark.parametrize('section', ['all', 'summary'])
    @pytest.mark.parametrize('max_states', [1, None])
    def test_print_results(self, tmp_path, section, max_states):

        tpa_results = self.compute_tpatransition('hf')

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        outfile = tmp_path / 'tpatransition_print.out'
        tpa_drv = TpaTransitionDriver(MPI.COMM_WORLD, OutputStream(outfile))
        tpa_drv.print_results(tpa_results, section=section, max_states=max_states)
        tpa_drv.ostream.flush()

        printed = outfile.read_text()
        assert 'TPA Transition Summary' in printed
        assert 'TPA Str. (Linear)' in printed
        assert 'TPA Str. (Circular)' in printed
        assert 'OPA Osc. Str.' in printed
        if section == 'all':
            assert 'Components of TPA Transition Moments' in printed
        elif section == 'summary':
            assert 'Components of TPA Transition Moments' not in printed
        if max_states == 1:
            assert 'additional state(s) omitted' in printed
        elif max_states is None:
            assert 'additional state(s) omitted' not in printed

    def test_plot_spectrum(self):

        pytest.importorskip('matplotlib')
        import matplotlib.pyplot as plt

        tpa_results = self.compute_tpatransition('hf')

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        fig, ax = plt.subplots()
        returned_ax = tpa_drv.plot_spectrum(tpa_results,
                                            x_unit='ev',
                                            broadening_value=0.123984,
                                            broadening_unit='ev',
                                            ax=ax)

        assert returned_ax is ax
        assert ax.get_xlabel() == 'Photon energy [eV]'
        assert ax.get_ylabel() == 'TPA cross-section [GM]'
        assert len(ax.lines) == 1
        assert len(ax.figure.axes) == 2
        plt.close(ax.figure)

    def test_plot_spectrum_circular(self):

        pytest.importorskip('matplotlib')
        import matplotlib.pyplot as plt

        tpa_results = self.compute_tpatransition('hf')

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        fig, ax = plt.subplots()
        returned_ax = tpa_drv.plot_spectrum(tpa_results,
                                            x_unit='ev',
                                            broadening_value=0.123984,
                                            broadening_unit='ev',
                                            polarization='circular',
                                            ax=ax)

        assert returned_ax is ax
        assert len(ax.figure.axes) == 2
        assert ax.figure.axes[1].get_ylabel() == 'TPA strengths (circular) [a.u.]'
        plt.close(ax.figure)

    def test_plot_spectrum_without_ax_returns_none(self):

        pytest.importorskip('matplotlib')

        tpa_results = self.compute_tpatransition('hf')

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()

        assert tpa_drv.plot_spectrum(tpa_results,
                                     x_unit='ev',
                                     broadening_value=0.123984,
                                     broadening_unit='ev') is None

    def test_tpatransition_results_hdf5_roundtrip(self, tmp_path):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        h5file = tmp_path / 'tpatransition_results.h5'
        with h5py.File(h5file, 'w'):
            pass

        tpa_results = {
            'photon_energies': [0.1656922537, 0.2119096339],
            'transition_moments': [
                np.array([[1.0 + 0.0j, 0.1 + 0.0j, 0.2 + 0.0j],
                          [0.1 + 0.0j, 2.0 + 0.0j, 0.3 + 0.0j],
                          [0.2 + 0.0j, 0.3 + 0.0j, 3.0 + 0.0j]]),
                np.array([[4.0 + 0.0j, 0.4 + 0.0j, 0.5 + 0.0j],
                          [0.4 + 0.0j, 5.0 + 0.0j, 0.6 + 0.0j],
                          [0.5 + 0.0j, 0.6 + 0.0j, 6.0 + 0.0j]])
            ],
            'tpa_strengths': {
                'linear': {
                    0: 4.1302840736,
                    1: 15.8381651150,
                },
                'circular': {
                    0: 1.5,
                    1: 2.5,
                },
            },
            'oscillator_strengths': np.array([0.01, 0.02]),
            'elec_trans_dipoles': np.array([[0.11, -0.22, 0.33],
                                            [0.44, -0.55, 0.66]]),
            'excitation_details': [['1a -> 2a (0.90)'], ['1b -> 3b (0.80)']],
            'rsp_type': 'tpa_transition',
        }

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        tpa_drv._write_final_hdf5(str(h5file), tpa_results)

        recovered = read_results(str(h5file), 'rsp')

        assert recovered['photon_energies'] == pytest.approx(
            tpa_results['photon_energies'])
        for rec_tensor, ref_tensor in zip(recovered['transition_moments'],
                                          tpa_results['transition_moments']):
            assert np.allclose(rec_tensor, ref_tensor)
        assert recovered['tpa_strengths'] == tpa_results['tpa_strengths']
        assert np.allclose(recovered['oscillator_strengths'],
                           tpa_results['oscillator_strengths'])
        assert np.allclose(recovered['elec_trans_dipoles'],
                           tpa_results['elec_trans_dipoles'])
        assert recovered['excitation_details'] == tpa_results[
            'excitation_details']

    def test_tpatransition_results_hdf5_roundtrip_without_suffix(self, tmp_path):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        h5stem = tmp_path / 'tpatransition_results'
        h5file = tmp_path / 'tpatransition_results.h5'
        with h5py.File(h5file, 'w'):
            pass

        tpa_results = {
            'photon_energies': [0.1656922537],
            'transition_moments': [
                np.array([[1.0 + 0.0j, 0.1 + 0.0j, 0.2 + 0.0j],
                          [0.1 + 0.0j, 2.0 + 0.0j, 0.3 + 0.0j],
                          [0.2 + 0.0j, 0.3 + 0.0j, 3.0 + 0.0j]])
            ],
            'tpa_strengths': {
                'linear': {
                    0: 4.1302840736,
                },
                'circular': {
                    0: 1.5,
                },
            },
            'oscillator_strengths': np.array([0.01]),
            'elec_trans_dipoles': np.array([[0.11, -0.22, 0.33]]),
            'excitation_details': [['1a -> 2a (0.90)']],
            'rsp_type': 'tpa_transition',
        }

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        tpa_drv._write_final_hdf5(str(h5stem), tpa_results)

        recovered = read_results(str(h5file), 'rsp')

        assert recovered['photon_energies'] == pytest.approx(
            tpa_results['photon_energies'])
        for rec_tensor, ref_tensor in zip(recovered['transition_moments'],
                                          tpa_results['transition_moments']):
            assert np.allclose(rec_tensor, ref_tensor)
        assert recovered['tpa_strengths'] == tpa_results['tpa_strengths']
        assert np.allclose(recovered['oscillator_strengths'],
                           tpa_results['oscillator_strengths'])
        assert np.allclose(recovered['elec_trans_dipoles'],
                           tpa_results['elec_trans_dipoles'])
        assert recovered['excitation_details'] == tpa_results[
            'excitation_details']

    def test_tpatransition_write_final_hdf5_without_filename_is_noop(self):

        tpa_results = {
            'photon_energies': [0.1656922537],
            'transition_moments': {
                0: np.eye(3),
            },
            'tpa_strengths': {
                'linear': {0: 4.1302840736},
                'circular': {0: 1.5},
            },
            'oscillator_strengths': np.array([0.01]),
            'elec_trans_dipoles': np.array([[0.11, -0.22, 0.33]]),
            'excitation_details': [['1a -> 2a (0.90)']],
        }

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        tpa_drv._write_final_hdf5(None, tpa_results)

    def test_update_settings_and_restart(self, tmp_path):

        tpa_drv = TpaTransitionDriver()
        tpa_drv.ostream.mute()
        assert tpa_drv.comm == MPI.COMM_WORLD

        tpa_dict = {
            'nstates': 2,
            'eri_thresh': 1e-13,
            'max_iter': 199,
            'conv_thresh': 1e-5,
            'lindep_thresh': 1e-11,
            'restart': False,
            'checkpoint_file': 'mycheckpoint.h5',
            'timing': True,
            'profiling': True,
            'memory_profiling': True,
            'memory_tracing': True,
        }

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) != val

        tpa_drv.update_settings(tpa_dict)

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) == val

        filename = str(tmp_path / 'tpatrans_restart')
        filename = MPI.COMM_WORLD.bcast(filename, root=mpi_master())

        scf_results, molecule, ao_basis = self.run_scf_restart('hf', filename)

        restart_drv = TpaTransitionDriver()
        restart_drv.ostream.mute()
        restart_drv.filename = filename
        restart_drv.nstates = 2
        restart_drv.conv_thresh = 1.0e-7
        first_results = restart_drv.compute(molecule, ao_basis, scf_results)
        assert restart_drv.checkpoint_file == f'{filename}_rsp.h5'

        restart_drv.restart = True
        restarted_results = restart_drv.compute(molecule, ao_basis, scf_results)
        assert restart_drv.restart

        fresh_drv = TpaTransitionDriver()
        fresh_drv.ostream.mute()
        fresh_drv.filename = filename
        fresh_drv.nstates = 2
        fresh_drv.conv_thresh = 1.0e-7
        fresh_drv.restart = False
        fresh_results = fresh_drv.compute(molecule, ao_basis, scf_results)
        assert not fresh_drv.restart

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            restarted_strengths = list(restarted_results['tpa_strengths']['linear'])
            first_strengths = list(first_results['tpa_strengths']['linear'])
            fresh_strengths = list(fresh_results['tpa_strengths']['linear'])

            assert restarted_strengths == pytest.approx(fresh_strengths,
                                                        abs=1.0e-7)
            assert first_strengths == pytest.approx(fresh_strengths,
                                                    abs=1.0e-7)

            for suffix in [
                    '.h5',
                    '_rsp_tpatrans_rpa.h5',
                    '_rsp_tpatrans_cpp.h5',
                    '_rsp_tpatrans_fock.h5',
            ]:
                assert Path(f'{filename}{suffix}').is_file()

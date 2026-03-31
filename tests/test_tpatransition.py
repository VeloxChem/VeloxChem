from mpi4py import MPI
from pathlib import Path
import pytest

from veloxchem.veloxchemlib import (mpi_master, hartree_in_ev,
                                    hartree_in_inverse_nm)
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
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

    def run_tpatransition(self, xcfun_label, ref_results, ri_coulomb=False):

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
            for (freq, val), (ref_freq, ref_val) in zip(tpa_str.items(),
                                                        ref_results.items()):
                assert abs(freq / ref_freq - 1.0) < 1.0e-6
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

        ref_result = {
            -0.1656922537003149: 4.130284073574981,
            -0.21190963394768278: 15.83816511497292,
        }
        self.run_tpatransition('hf', ref_result)

    def test_tpatransition_lda(self):

        ref_result = {
            -0.13269116242012727: 12.922686324743287,
            -0.18051942347838879: 55.35603386346646,
        }
        self.run_tpatransition('slda', ref_result)

    def test_tpatransition_gga(self):

        ref_result = {
            -0.13304132289794054: 13.068064328595725,
            -0.179389130413857: 57.062788455439296,
        }
        self.run_tpatransition('bp86', ref_result)

    def test_tpatransition_mgga(self):

        ref_result = {
            -0.1381414072229231: 9.579957827267322,
            -0.18265658100098003: 41.94340887233736,
        }
        self.run_tpatransition('tpssh', ref_result)

    def test_tpatransition_ri_blyp(self):

        ref_result = {
            -0.12580806392701718: 13.672823155898364,
            -0.1710157238983142: 60.526749482300374,
        }
        self.run_tpatransition('blyp', ref_result, ri_coulomb=True)

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
            restarted_strengths = list(
                restarted_results['tpa_strengths']['linear'].values())
            first_strengths = list(first_results['tpa_strengths']['linear'].values())
            fresh_strengths = list(fresh_results['tpa_strengths']['linear'].values())

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

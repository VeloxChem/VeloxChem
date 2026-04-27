from mpi4py import MPI
from pathlib import Path
import pytest

from veloxchem.veloxchemlib import (mpi_master, bohr_in_angstrom,
                                    hartree_in_ev, hartree_in_inverse_nm,
                                    fine_structure_constant,
                                    speed_of_light_in_vacuum_in_SI)
from veloxchem.outputstream import OutputStream
from veloxchem.tpadriver import TpaDriver
from veloxchem.tpafulldriver import TpaFullDriver
from veloxchem.tpareddriver import TpaReducedDriver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rsptpa import TPA


@pytest.mark.solvers
class TestTPA:

    def run_scf(self, task):

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        return scf_results

    def run_scf_restart(self, task, filename):

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.filename = filename
        scf_drv.ostream.mute()
        scf_drv.compute(task.molecule, task.ao_basis)
        scf_drv.restart = True
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        assert scf_drv.restart

        return scf_results

    def run_tpa(self, inpfile, tpa_type, w, ref_result):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_results = self.run_scf(task)

        tpa_prop = TPA({
            'damping': task.input_dict['response']['damping'],
            'frequencies': task.input_dict['response']['frequencies'],
            'conv_thresh': '1.0e-8',
            'tpa_type': tpa_type,
        })
        tpa_prop.init_driver(task.mpi_comm, task.ostream)
        tpa_prop.compute(task.molecule, task.ao_basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            tpa_result = tpa_prop.rsp_property

            for key in ref_result:
                assert abs(tpa_result[key][(w, -w, w)].real /
                           ref_result[key].real - 1.0) < 1.0e-6
                assert abs(tpa_result[key][(w, -w, w)].imag /
                           ref_result[key].imag - 1.0) < 1.0e-6

    def test_tpa_full(self):

        # vlxtag: RHF, TPA, CR

        w = 0.05

        ref_result = {
            't4_dict': 11.43071305 + 0.04957732j,
            't3_dict': -42.19841751 - 0.28695214j,
            'NaX3NyNz': -81.62345190 - 0.35812832j,
            'NaA3NxNy': -27.21320341 - 0.03029788j,
            'NaX2Nyz': 270.69041328 + 2.67837597j,
            'NxA2Nyz': 270.83461366 + 0.52758094j,
            'gamma': 401.92066716 + 2.58015589j,
        }

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_tpa.inp')

        self.run_tpa(inpfile, 'full', w, ref_result)

    def test_tpa_reduced(self):

        # vlxtag: RHF, TPA, CR

        w = 0.05

        ref_result = {
            't3_dict': -15.12982062 - 0.19793495j,
            'NaX2Nyz': 96.30910639 + 1.72679037j,
            'NxA2Nyz': 96.36431088 + 0.51886895j,
            'gamma': 177.54359664 + 2.04772438j,
        }

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_tpa.inp')

        self.run_tpa(inpfile, 'reduced', w, ref_result)

    def test_update_settings(self):

        tpa_dict = {
            'frequencies': (0.1, 0.12, 0.16, 0.20),
            'damping': 0.01,
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

        tpa_drv = TpaDriver(MPI.COMM_WORLD, OutputStream(None))

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) != val

        tpa_drv.update_settings(tpa_dict)

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) == val

    def test_get_spectrum(self):

        alpha = fine_structure_constant()
        a0_in_cm = bohr_in_angstrom() * 1.0e-8
        c_in_cm_per_s = speed_of_light_in_vacuum_in_SI() * 100.0
        au2gm = ((4.0 * pytest.importorskip('numpy').pi**2 * alpha *
                  a0_in_cm**5) / c_in_cm_per_s * 1.0e+50)

        rsp_results = {
            'gamma': {
                (0.0, -0.0, 0.0): 1.0 + 2.0j,
                (0.05, -0.05, 0.05): 3.0 + 4.0j,
                (0.10, -0.10, 0.10): 5.0 + 6.0j,
            },
            'frequencies': [0.0, 0.05, 0.10],
        }

        spectrum_au = TpaDriver.get_spectrum(rsp_results, 'au')
        spectrum_ev = TpaDriver.get_spectrum(rsp_results, 'ev')
        spectrum_nm = TpaDriver.get_spectrum(rsp_results, 'nm')

        assert spectrum_au['x_label'] == 'Photon energy [a.u.]'
        assert spectrum_ev['x_label'] == 'Photon energy [eV]'
        assert spectrum_nm['x_label'] == 'Wavelength [nm]'
        assert spectrum_au['y_label'] == 'TPA cross-section [GM]'

        assert spectrum_au['x_data'] == pytest.approx([0.05, 0.10])
        assert spectrum_ev['x_data'] == pytest.approx(
            [0.05 * hartree_in_ev(), 0.10 * hartree_in_ev()])
        assert spectrum_nm['x_data'] == pytest.approx([
            1.0 / (hartree_in_inverse_nm() * 0.05),
            1.0 / (hartree_in_inverse_nm() * 0.10),
        ])

        expected_y_data = [
            4.0j.imag * 0.05**2 * au2gm,
            6.0j.imag * 0.10**2 * au2gm,
        ]

        assert spectrum_au['y_data'] == pytest.approx(expected_y_data)
        assert spectrum_ev['y_data'] == pytest.approx(expected_y_data)
        assert spectrum_nm['y_data'] == pytest.approx(expected_y_data)

    def test_full_driver_defaults_and_restart(self, tmp_path):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_tpa.inp')

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        filename = str(tmp_path / 'tpa_full_restart')
        filename = task.mpi_comm.bcast(filename, root=mpi_master())

        scf_results = self.run_scf_restart(task, filename)

        rsp_settings = {
            'damping': task.input_dict['response']['damping'],
            'frequencies': task.input_dict['response']['frequencies'],
            'conv_thresh': 1.0e-8,
        }

        tpa_drv = TpaFullDriver()
        tpa_drv.ostream.mute()
        tpa_drv.filename = filename
        tpa_drv.update_settings(rsp_settings)
        first_results = tpa_drv.compute(task.molecule, task.ao_basis,
                                        scf_results)
        assert tpa_drv.checkpoint_file == f'{filename}_rsp.h5'

        tpa_drv.restart = True
        restarted_results = tpa_drv.compute(task.molecule, task.ao_basis,
                                            scf_results)
        assert tpa_drv.restart

        fresh_drv = TpaFullDriver()
        fresh_drv.ostream.mute()
        fresh_drv.filename = filename
        fresh_drv.update_settings(rsp_settings)
        fresh_drv.restart = False
        fresh_results = fresh_drv.compute(task.molecule, task.ao_basis,
                                          scf_results)
        assert not fresh_drv.restart

        if task.mpi_rank == mpi_master():
            key = (0.05, -0.05, 0.05)

            for result_key in ['t4_dict', 't3_dict', 'NaX3NyNz', 'NaA3NxNy',
                               'NaX2Nyz', 'NxA2Nyz', 'gamma']:
                assert restarted_results[result_key][key] == pytest.approx(
                    fresh_results[result_key][key], abs=1.0e-7)
                assert first_results[result_key][key] == pytest.approx(
                    fresh_results[result_key][key], abs=1.0e-7)

            for suffix in [
                    '.h5',
                    '_rsp_tpa_1.h5',
                    '_rsp_tpa_2_full.h5',
                    '_rsp_tpa_fock_1_full.h5',
                    '_rsp_tpa_fock_2_full.h5',
            ]:
                assert Path(f'{filename}{suffix}').is_file()

    def test_reduced_driver_defaults_and_restart(self, tmp_path):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_tpa.inp')

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        filename = str(tmp_path / 'tpa_reduced_restart')
        filename = task.mpi_comm.bcast(filename, root=mpi_master())

        scf_results = self.run_scf_restart(task, filename)

        rsp_settings = {
            'damping': task.input_dict['response']['damping'],
            'frequencies': task.input_dict['response']['frequencies'],
            'conv_thresh': 1.0e-8,
        }

        tpa_drv = TpaReducedDriver()
        tpa_drv.ostream.mute()
        tpa_drv.filename = filename
        tpa_drv.update_settings(rsp_settings)
        first_results = tpa_drv.compute(task.molecule, task.ao_basis,
                                        scf_results)
        assert tpa_drv.checkpoint_file == f'{filename}_rsp.h5'

        tpa_drv.restart = True
        restarted_results = tpa_drv.compute(task.molecule, task.ao_basis,
                                            scf_results)
        assert tpa_drv.restart

        fresh_drv = TpaReducedDriver()
        fresh_drv.ostream.mute()
        fresh_drv.filename = filename
        fresh_drv.update_settings(rsp_settings)
        fresh_drv.restart = False
        fresh_results = fresh_drv.compute(task.molecule, task.ao_basis,
                                          scf_results)
        assert not fresh_drv.restart

        if task.mpi_rank == mpi_master():
            key = (0.05, -0.05, 0.05)

            for result_key in ['t3_dict', 'NaX2Nyz', 'NxA2Nyz', 'gamma']:
                assert restarted_results[result_key][key] == pytest.approx(
                    fresh_results[result_key][key], abs=1.0e-7)
                assert first_results[result_key][key] == pytest.approx(
                    fresh_results[result_key][key], abs=1.0e-7)

            for suffix in [
                    '.h5',
                    '_rsp_tpa_1.h5',
                    '_rsp_tpa_2_red.h5',
                    '_rsp_tpa_fock_1_red.h5',
                    '_rsp_tpa_fock_2_red.h5',
            ]:
                assert Path(f'{filename}{suffix}').is_file()

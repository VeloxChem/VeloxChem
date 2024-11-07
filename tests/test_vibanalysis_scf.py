from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestScfVibrationalAnalysisDriver:

    def test_scf_vibrational_analysis_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        vib_settings = {'do_ir': 'yes', 'do_raman': 'yes', 'numerical_hessian':'yes',
                        'numerical_raman': 'yes'}
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if is_mpi_master(task.mpi_comm):

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_intensities = np.array(hf.get('raman'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))
            rel_diff_raman = np.max(
                np.abs(vibanalysis_drv.raman_intensities[0.0] /
                       ref_raman_intensities - 1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3
            assert rel_diff_raman < 1.0e-3

        task.finish()

    def test_scf_resonance_raman(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        vib_settings = {'do_ir': 'no', 'do_resonance_raman': 'yes', 'numerical_hessian':'yes',
                        'numerical_raman': 'yes', 'frequencies': (0.0, 0.4), 'rr_damping': 0.05}
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if is_mpi_master(task.mpi_comm):

            hf = h5py.File(h5file)
            ref_frequencies = np.array(hf.get('frequencies'))
            hf_rr = hf['resonance_raman']
            ref_raman_intensities = np.array([hf_rr.get('0.0'), hf_rr.get('0.4')])
            hf.close()

            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_raman_static = np.max(
                np.abs(vibanalysis_drv.raman_intensities[0.0] /
                       ref_raman_intensities[0] - 1.0))
            rel_diff_raman_dyn = np.max(
                np.abs(vibanalysis_drv.raman_intensities[0.4] /
                       ref_raman_intensities[1] - 1.0))

            assert rel_diff_freq < 1.0e-3
            assert rel_diff_raman_static < 1.0e-3
            assert rel_diff_raman_dyn < 1.0e-3


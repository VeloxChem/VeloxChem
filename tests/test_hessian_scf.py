from pathlib import Path
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestScfHessianDriver:

    @pytest.mark.timeconsuming
    def test_scf_hessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_hessian_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        vib_settings = {
            'do_ir': 'yes',
            'do_raman': 'yes',
            'numerical_hessian': 'yes',
            'numerical_raman': 'yes'
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_activities = np.array(hf.get('raman'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))
            rel_diff_raman = np.max(
                np.abs(vibanalysis_drv.raman_activities[0.0] /
                       ref_raman_activities - 1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3
            assert rel_diff_raman < 1.0e-3

        task.finish()

    @pytest.mark.solvers
    def test_analytical_scf_hessian(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_analytical_hessian_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        method_settings = {}
        cphf_settings = {'conv_thresh': 1e-8}
        hess_settings = {'numerical': 'no'}
        hessian_drv = ScfHessianDriver(scf_drv)
        hessian_drv.update_settings(method_settings,
                                    hess_dict=hess_settings,
                                    cphf_dict=cphf_settings)
        hessian_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()
            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            assert diff_hessian < 1e-6

        hessian_drv.use_subcomms = True
        hessian_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()
            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            assert diff_hessian < 1e-6

    @pytest.mark.solvers
    def test_analytical_pbe_hessian(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_analytical_hessian_pbe.h5')

        task = MpiTask([inpfile, None])

        method_settings = {'xcfun': 'pbe', 'grid_level': 7}
        scf_settings = {}
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        cphf_settings = {'conv_thresh': 1e-8}
        hess_settings = {'numerical': 'no'}

        hessian_drv = ScfHessianDriver(scf_drv)
        hessian_drv.update_settings(method_settings,
                                    hess_dict=hess_settings,
                                    cphf_dict=cphf_settings)
        hessian_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()
            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            assert diff_hessian < 1e-6

        hessian_drv.use_subcomms = True
        hessian_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()
            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            assert diff_hessian < 1e-6

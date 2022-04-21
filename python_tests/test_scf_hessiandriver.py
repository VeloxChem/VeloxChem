from pathlib import Path
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver


class TestScfHessianDriver:

    def test_scfhessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_hessian.inp')
        h5file = str(here / 'inputs' / 'scfhessian.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():
            hessian_settings = {'do_raman':'yes'}
            method_settings = {}
            scf_hessian_drv = ScfHessianDriver(scf_drv)
            scf_hessian_drv.update_settings(method_settings, hessian_settings)
            scf_hessian_drv.compute(task.molecule, task.ao_basis)

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_intensities = np.array(hf.get('raman'))
            hf.close()

            assert np.max(np.abs(scf_hessian_drv.hessian - ref_hessian)) < 1.0e-5

        task.finish()

    @pytest.mark.skipif(True,
       reason='geometric verison not able to run vib. analysis')
    def test_scf_vibrational_analysis(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_hessian.inp')
        h5file = str(here / 'inputs' / 'scfhessian.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():
            hessian_settings = {'do_raman':'yes'}
            method_settings = {}
            scf_hessian_drv = ScfHessianDriver(scf_drv)
            scf_hessian_drv.update_settings(method_settings, hessian_settings)
            scf_hessian_drv.compute(task.molecule, task.ao_basis)

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_intensities = np.array(hf.get('raman'))
            hf.close()

            # This works only with a new enough version of geometric
            # which can do vibrational analysis and can
            # return un-normalized normal modes (needed for IR intensities
            # and Raman activities)
            scf_hessian_drv.vibrational_analysis(task.molecule)
            assert np.max(np.abs(scf_hessian_drv.frequencies 
                                 - ref_frequencies)) < 1.0e-3
            assert np.max(np.abs(scf_hessian_drv.ir_intensities
                                 - ref_ir_intensities)) < 1.0e-4
            assert np.max(np.abs(scf_hessian_drv.raman_intensities
                                 - ref_raman_intensities)) < 1.0e-4
            
        task.finish()

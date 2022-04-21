from pathlib import Path
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.xtbdriver import XTBDriver
from veloxchem.xtbhessiandriver import XTBHessianDriver


class TestXTBHessianDriver:

    @pytest.mark.skipif(not XTBDriver().is_available(),
                        reason='xtb not available')
    def test_xtbhessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_hessian.inp')
        h5file = str(here / 'inputs' / 'water_hessian_xtb.h5')

        task = MpiTask([inpfile, None])

        xtb_method = 'gfn2'

        xtb_drv = XTBDriver(task.mpi_comm)
        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule, task.ostream)

        xtb_hessian_drv = XTBHessianDriver(xtb_drv)
        xtb_hessian_drv.compute(task.molecule)

        if is_mpi_master(task.mpi_comm):
            # This works only with a new enough version of geometric
            # which can do vibrational analysis and can
            # return un-normalized normal modes (needed for IR intensities
            # and Raman activities)
            xtb_hessian_drv.vibrational_analysis(task.molecule)

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            #ref_frequencies = np.array(hf.get('frequencies'))
            #ref_ir_intensities = np.array(hf.get('ir'))
            hf.close()

            assert np.max(
                np.abs(xtb_hessian_drv.hessian - ref_hessian)) < 1.0e-5
            #assert np.max(np.abs(xtb_hessian_drv.frequencies
            #                     - ref_frequencies)) < 1.0e-3
            #assert np.max(np.abs(xtb_hessian_drv.ir_intensities
            #                     - ref_ir_intensities)) < 1.0e-4

        task.finish()

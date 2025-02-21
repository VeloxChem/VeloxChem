from pathlib import Path
import numpy as np
import pytest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.xtbdriver import XtbDriver
from veloxchem.xtbhessiandriver import XtbHessianDriver


class TestXtbHessianDriver:

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='xtb not available')
    def test_xtb_hessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_xtb.inp')
        h5file = str(here / 'data' / 'water_hessian_xtb.h5')

        task = MpiTask([inpfile, None])

        xtb_method = 'gfn2'

        xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
        xtb_drv.ostream.mute()

        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule)

        xtb_hessian_drv = XtbHessianDriver(xtb_drv)
        xtb_hessian_drv.ostream.mute()
        xtb_hessian_drv.compute(task.molecule)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()

            diff_hessian = np.max(np.abs(xtb_hessian_drv.hessian - ref_hessian))

            assert diff_hessian < 1.0e-5

        task.finish()

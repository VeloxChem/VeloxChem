from pathlib import Path

import h5py
import numpy as np
from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.xtbdriver import XtbDriver
from veloxchem.xtbhessiandriver import XtbHessianDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis

from .addons import using_xtb


class TestXtbHessianDriver:

    @using_xtb
    def test_xtb_hessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_hessian_xtb.inp')
        h5file = str(here / 'inputs' / 'water_hessian_xtb.h5')

        task = MpiTask([inpfile, None])

        xtb_method = 'gfn2'

        xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
        xtb_drv.ostream.mute()

        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule)

        xtb_hessian_drv = XtbHessianDriver(xtb_drv)

        vibanalysis_drv = VibrationalAnalysis(xtb_hessian_drv)
        vibanalysis_drv.update_settings(vib_dict = vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule)

        if is_mpi_master(task.mpi_comm):

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3

        task.finish()

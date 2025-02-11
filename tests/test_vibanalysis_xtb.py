from pathlib import Path
import numpy as np
import pytest
import h5py

from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import mpi_master
from veloxchem.xtbdriver import XtbDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestXtbVibrationalAnalysisDriver:

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='xtb not available')
    def test_xtb_vibrational_analysis_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_xtb.inp')
        h5file = str(here / 'data' / 'water_hessian_xtb.h5')

        task = MpiTask([inpfile, None])

        xtb_method = 'gfn2'

        xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
        xtb_drv.ostream.mute()

        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule)

        vib_settings = {'do_ir': 'yes', 'numerical_hessian': 'yes'}
        vibanalysis_drv = VibrationalAnalysis(xtb_drv)
        vibanalysis_drv.update_settings(vib_dict=vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule)

        if task.mpi_rank == mpi_master():

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

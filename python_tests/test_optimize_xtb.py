from pathlib import Path

import numpy as np
import pytest
from veloxchem.mpitask import MpiTask
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.veloxchemlib import bohr_in_angstrom, is_mpi_master
from veloxchem.xtbdriver import XtbDriver
from veloxchem.xtbgradientdriver import XtbGradientDriver

from .addons import using_xtb


@using_xtb
@pytest.mark.filterwarnings(
    'ignore:.*tostring.*tobytes:DeprecationWarning:geometric')
class TestOptimizeXTB:

    def run_opt(self, inpfile, xtb_method, ref_coords):

        task = MpiTask([inpfile, None])

        xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
        xtb_drv.ostream.mute()

        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule)

        grad_drv = XtbGradientDriver(task.mpi_comm, task.ostream)
        opt_drv = OptimizationDriver(grad_drv)
        opt_drv.update_settings({
            'coordsys': 'tric',
            'filename': task.input_dict['filename'],
            'keep_files': 'no',
        })
        opt_mol = opt_drv.compute(task.molecule, xtb_drv)

        if is_mpi_master(task.mpi_comm):
            opt_coords = opt_mol.get_coordinates_in_bohr()
            assert np.max(np.abs(opt_coords - ref_coords)) < 1.0e-6

            inpfile = Path(inpfile)
            optfile = Path(str(inpfile.with_name(inpfile.stem)) + '_optim.xyz')
            logfile = inpfile.with_suffix('.log')
            if optfile.is_file():
                optfile.unlink()
            if logfile.is_file():
                logfile.unlink()

    def test_nh3(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'nh3.inp')

        xtb_method = 'gfn2'

        ref_coords = np.array([
            [-1.981952536628, 1.605526632230, -0.018764716881],
            [-1.950824053259, 2.616502833048, 0.030950778306],
            [-2.483737817114, 1.273354861268, 0.795720974820],
            [-2.523925592999, 1.354865673454, -0.836657036245],
        ]) / bohr_in_angstrom()

        self.run_opt(inpfile, xtb_method, ref_coords)

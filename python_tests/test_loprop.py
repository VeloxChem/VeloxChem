from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.loprop import LoPropDriver


class TestLoProp:

    def run_loprop(self, inpfile, ref_charges, ref_polarizabilities):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        loprop_drv = LoPropDriver(task.mpi_comm, task.ostream)
        loprop_result = loprop_drv.compute(task.molecule, task.ao_basis,
                                           scf_drv.scf_tensors)

        if is_mpi_master(task.mpi_comm):
            charges = loprop_result['localized_charges']
            polarizabilities = loprop_result['localized_polarizabilities']
            assert np.max(np.abs(charges - ref_charges)) < 1.0e-4
            assert np.max(
                np.abs(polarizabilities - ref_polarizabilities)) < 1.0e-4

    def test_loprop_water(self):

        # vlxtag: RHF, LoProp

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'loprop_water.inp')

        ref_charges = np.array([-0.1198, 0.0599, 0.0599])
        ref_polarizabilities = np.array(
            [[3.1116, -0.0000, -0.1131, 0.8247, -0.0000, 3.1116],
             [3.3244, -0.0000, -0.2019, 0.2486, 0.0000, 0.1336],
             [0.1336, -0.0000, -0.2019, 0.2486, -0.0000, 3.3244]])

        self.run_loprop(inpfile, ref_charges, ref_polarizabilities)

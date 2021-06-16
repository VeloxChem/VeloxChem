from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.respchargesdriver import RespChargesDriver


class TestRespCharges:

    def run_resp(self, inpfile, ref_resp_charges, ref_esp_charges):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings({'number_layers': 1},
                                task.input_dict['method_settings'])

        q_resp = chg_drv.compute(task.molecule, task.ao_basis, 'resp')
        q_esp = chg_drv.compute(task.molecule, task.ao_basis, 'esp')

        if is_mpi_master(task.mpi_comm):
            resp_charges = q_resp
            esp_charges = q_esp
            assert np.max(np.abs(resp_charges - ref_resp_charges)) < 1.0e-6
            assert np.max(np.abs(esp_charges - ref_esp_charges)) < 1.0e-6

    def test_resp_methanol(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol.inp')

        ref_resp_charges = np.array(
            [0.041310, 0.021227, 0.041310, 0.041310, -0.454284, 0.309127])
        ref_esp_charges = np.array(
            [0.038183, 0.091353, 0.018171, 0.013240, -0.474887, 0.313940])

        self.run_resp(inpfile, ref_resp_charges, ref_esp_charges)

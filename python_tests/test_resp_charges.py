from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.respchargesdriver import RespChargesDriver


class TestRespCharges:

    def run_resp(self, inpfile, ref_charges, inp_chg_dict, chg_type):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        chg_dict = {'filename': task.input_dict['filename']}
        chg_dict.update(inp_chg_dict)

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(chg_dict, task.input_dict['method_settings'])

        q_fit = chg_drv.compute(task.molecule, task.ao_basis, scf_results,
                                chg_type.lower())

        if is_mpi_master(task.mpi_comm):
            assert np.max(np.abs(q_fit - ref_charges)) < 1.0e-5

            pdb_file = Path(chg_drv.filename).with_suffix('.pdb')
            if pdb_file.is_file():
                pdb_file.unlink()
            scf_h5_file = Path(chg_drv.filename).with_suffix('.scf.h5')
            if scf_h5_file.is_file():
                scf_h5_file.unlink()
            scf_final_h5_file = Path(
                chg_drv.filename).with_suffix('.scf.results.h5')
            if scf_final_h5_file.is_file():
                scf_final_h5_file.unlink()

    def test_resp_methanol(self):

        # vlxtag: RHF, RESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol.inp')

        ref_resp_charges = np.array(
            [0.041310, 0.021227, 0.041310, 0.041310, -0.454284, 0.309127])

        chg_dict = {'number_layers': 1}

        self.run_resp(inpfile, ref_resp_charges, chg_dict, 'resp')

    def test_esp_methanol(self):

        # vlxtag: RHF, ESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol.inp')

        ref_esp_charges = np.array(
            [0.038183, 0.091353, 0.018171, 0.013240, -0.474887, 0.313940])

        chg_dict = {'number_layers': 1}

        self.run_resp(inpfile, ref_esp_charges, chg_dict, 'esp')

    def test_chelpg_esp_methanol(self):

        # vlxtag: RHF, ESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol.inp')

        ref_chelpg_esp_charges = np.array(
            [0.007755, 0.145603, 0.007755, 0.007755, -0.468499, 0.299632])

        chg_dict = {'grid_type': 'chelpg', 'equal_charges': '1=3=4'}

        self.run_resp(inpfile, ref_chelpg_esp_charges, chg_dict, 'esp')

    def test_esp_points_water(self):

        # vlxtag: RHF, ESP_on_points

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_esp.inp')

        ref_esp_charges = np.array([-1.106527, 1.106527])
        esp_fitting_points = [
            '0.0000000    0.0000000   -0.1653507',
            '0.0000000    0.0000000    0.4424329',
        ]
        chg_dict = {'number_layers': 4, 'fitting_points': esp_fitting_points}

        self.run_resp(inpfile, ref_esp_charges, chg_dict, 'esp')

from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.espchargesdriver import EspChargesDriver
from veloxchem.resultsio import read_results


class TestEspCharges:

    def run_esp(self,
                inpfile,
                ref_charges,
                inp_chg_dict,
                custom_mk_radii=None,
                *,
                ri):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.ri = ri
        scf_drv.filename = task.input_dict['filename']
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        chg_dict = {'filename': task.input_dict['filename']}
        chg_dict.update(inp_chg_dict)

        chg_drv = EspChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(chg_dict, task.input_dict['method_settings'])

        if custom_mk_radii is not None:
            chg_drv.custom_mk_radii = custom_mk_radii

        q_fit = chg_drv.compute(task.molecule, task.ao_basis, scf_results)

        if task.mpi_rank == mpi_master():
            assert np.max(np.abs(q_fit - ref_charges)) < 1.0e-5
            esp_h5_results = read_results(f'{chg_drv.filename}.h5', 'esp')
            if chg_drv.fitting_points is not None:
                np.testing.assert_allclose(esp_h5_results['esp_on_points'],
                                           q_fit)
            else:
                np.testing.assert_allclose(esp_h5_results['esp_charges'],
                                           q_fit)

            pdb_file = Path(chg_drv.filename).with_suffix('.pdb')
            pdb_file.unlink(missing_ok=True)

            final_h5_file = Path(chg_drv.filename).with_suffix('.h5')
            final_h5_file.unlink(missing_ok=True)

            scf_h5_file = Path(chg_drv.filename + '_scf.h5')
            scf_h5_file.unlink(missing_ok=True)

    @pytest.mark.parametrize('ri', [False, True])
    def test_esp_methanol(self, ri):

        # vlxtag: RHF, ESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol.inp')

        ref_esp_charges = {
            False: np.array(
                [0.038183, 0.091353, 0.018171, 0.013240, -0.474887, 0.313940]),
            True: np.array([
                0.0381445, 0.0914522, 0.0181317, 0.0132006, -0.4748232,
                0.3138942
            ])
        }[ri]

        chg_dict = {'number_layers': 1}

        self.run_esp(inpfile, ref_esp_charges, chg_dict, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_chelpg_esp_methanol(self, ri):

        # vlxtag: RHF, ESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol.inp')

        ref_chelpg_esp_charges = {
            False: np.array(
                [0.007755, 0.145603, 0.007755, 0.007755, -0.468499, 0.299632]),
            True: np.array([
                0.0077154, 0.1457034, 0.0077154, 0.0077154, -0.4684355,
                0.2995860
            ])
        }[ri]

        chg_dict = {'grid_type': 'chelpg', 'equal_charges': '1=3=4'}

        self.run_esp(inpfile, ref_chelpg_esp_charges, chg_dict, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_esp_points_water(self, ri):

        # vlxtag: RHF, ESP_on_points

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_esp.inp')

        ref_esp_charges = {
            False: np.array([-1.106527, 1.106527]),
            True: np.array([-1.1063509, 1.1063509])
        }[ri]
        esp_fitting_points = [
            '0.0000000    0.0000000   -0.1653507',
            '0.0000000    0.0000000    0.4424329',
        ]
        chg_dict = {'number_layers': 4, 'fitting_points': esp_fitting_points}

        self.run_esp(inpfile, ref_esp_charges, chg_dict, ri=ri)

    def test_esp_writes_h5_without_input_scf_results(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_esp.inp')

        task = MpiTask([inpfile, None])
        chg_dict = dict(task.input_dict.get('esp_charges', {}))
        chg_dict['filename'] = task.input_dict['filename']

        chg_drv = EspChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(chg_dict, task.input_dict['method_settings'])

        q_fit = chg_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            esp_h5_results = read_results(f'{chg_drv.filename}.h5', 'esp')
            np.testing.assert_allclose(esp_h5_results['esp_on_points'], q_fit)

            pdb_file = Path(chg_drv.filename).with_suffix('.pdb')
            pdb_file.unlink(missing_ok=True)

            final_h5_file = Path(chg_drv.filename).with_suffix('.h5')
            final_h5_file.unlink(missing_ok=True)

            scf_h5_file = Path(chg_drv.filename + '_scf.h5')
            scf_h5_file.unlink(missing_ok=True)

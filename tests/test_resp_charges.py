from pathlib import Path
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.respchargesdriver import RespChargesDriver
from veloxchem.resultsio import read_results


class TestRespCharges:

    def run_resp(self,
                 inpfile,
                 ref_charges,
                 inp_chg_dict,
                 chg_type,
                 custom_mk_radii=None):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.filename = task.input_dict['filename']
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        chg_dict = {'filename': task.input_dict['filename']}
        chg_dict.update(inp_chg_dict)

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(chg_dict, task.input_dict['method_settings'])

        if custom_mk_radii is not None:
            chg_drv.custom_mk_radii = custom_mk_radii

        q_fit = chg_drv.compute(task.molecule, task.ao_basis, scf_results,
                                chg_type.lower())

        if task.mpi_rank == mpi_master():
            assert np.max(np.abs(q_fit - ref_charges)) < 1.0e-5
            resp_h5_results = read_results(f'{chg_drv.filename}.h5', 'resp')
            np.testing.assert_allclose(resp_h5_results['resp_charges'], q_fit)

            pdb_file = Path(chg_drv.filename).with_suffix('.pdb')
            pdb_file.unlink(missing_ok=True)

            final_h5_file = Path(chg_drv.filename).with_suffix('.h5')
            final_h5_file.unlink(missing_ok=True)

            scf_h5_file = Path(chg_drv.filename + '_scf.h5')
            scf_h5_file.unlink(missing_ok=True)

    def test_resp_methanol(self):

        # vlxtag: RHF, RESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol.inp')

        ref_resp_charges = np.array(
            [0.041310, 0.021227, 0.041310, 0.041310, -0.454284, 0.309127])

        chg_dict = {'number_layers': 1}

        self.run_resp(inpfile, ref_resp_charges, chg_dict, 'resp')

    def test_resp_methanol_custom_radii(self):

        # vlxtag: RHF, RESP_charges

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol.inp')

        ref_resp_charges = np.array([
            0.04858359, -0.00613574, 0.04858359, 0.04858359, -0.43473053,
            0.29511548
        ])

        chg_dict = {'number_layers': 1}

        self.run_resp(inpfile, ref_resp_charges, chg_dict, 'resp', ['C', 3.0])

    def test_resp_writes_h5_without_input_scf_results(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol.inp')

        task = MpiTask([inpfile, None])
        chg_dict = {'filename': task.input_dict['filename'], 'number_layers': 1}

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(chg_dict, task.input_dict['method_settings'])

        q_fit = chg_drv.compute(task.molecule, task.ao_basis, None, 'resp')

        if task.mpi_rank == mpi_master():
            resp_h5_results = read_results(f'{chg_drv.filename}.h5', 'resp')
            np.testing.assert_allclose(resp_h5_results['resp_charges'], q_fit)

            pdb_file = Path(chg_drv.filename).with_suffix('.pdb')
            pdb_file.unlink(missing_ok=True)

            final_h5_file = Path(chg_drv.filename).with_suffix('.h5')
            final_h5_file.unlink(missing_ok=True)

            scf_h5_file = Path(chg_drv.filename + '_scf.h5')
            scf_h5_file.unlink(missing_ok=True)

    def test_get_origin_and_charge_dipole(self):

        # H2 with bond length 1 bohr and charges [+0.5, -0.5].

        mol_str = 'H  0.0  0.0  0.0\nH  1.0  0.0  0.0'
        molecule = Molecule.read_str(mol_str, units='au')

        charges = np.array([0.5, -0.5])

        chg_drv = RespChargesDriver()
        origin, dipole = chg_drv._get_origin_and_charge_dipole(molecule, charges)

        assert np.allclose(origin, np.array([0.5, 0.0, 0.0]), atol=1.0e-10)
        assert np.allclose(dipole, np.array([-0.5, 0.0, 0.0]), atol=1.0e-10)

        # Consistency check: the dipole should be translation-invariant since
        # the origin moves with the molecule.

        mol_str_shifted = 'H  5.0  3.0  2.0\nH  6.0  3.0  2.0'
        molecule_shifted = Molecule.read_str(mol_str_shifted, units='au')

        origin_shifted, dipole_shifted = (
            chg_drv._get_origin_and_charge_dipole(molecule_shifted, charges))

        assert np.allclose(origin_shifted, np.array([5.5, 3.0, 2.0]),
                           atol=1.0e-10)
        assert np.allclose(dipole, dipole_shifted, atol=1.0e-10)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_origin_and_charge_dipole_wrong_charge_count(self):

        mol_str = 'H  0.0  0.0  0.0\nH  1.0  0.0  0.0'
        molecule = Molecule.read_str(mol_str, units='au')

        chg_drv = RespChargesDriver()

        with pytest.raises(ValueError, match='Expected 2 charges but got 3'):
            chg_drv._get_origin_and_charge_dipole(molecule, np.zeros(3))

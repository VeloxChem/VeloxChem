from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.peforcefieldgenerator import PEForceFieldGenerator


@pytest.mark.solvers
class TestPEForceFieldGenerator:

    def run_loprop(self, inpfile, ref_charges, ref_polarizabilities):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        pe_ff_gen = PEForceFieldGenerator(task.mpi_comm, task.ostream)
        pe_ff_results = pe_ff_gen.compute(task.molecule, task.ao_basis,
                                          scf_results)

        if task.mpi_rank == mpi_master():
            charges = pe_ff_results['localized_charges']
            polarizabilities = pe_ff_results['localized_polarizabilities']

            assert np.max(np.abs(charges - ref_charges)) < 1.0e-4
            assert np.max(np.abs(polarizabilities -
                                 ref_polarizabilities)) < 1.0e-4

    def test_loprop_water(self):

        # vlxtag: RHF, PEForceFieldGenerator

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'loprop_water.inp')

        ref_charges = np.array([-0.66166, 0.33083, 0.33083])
        ref_polarizabilities = np.array(
            [[4.10442, 0., 0.08114, 2.90629, -0., 4.10442],
             [3.17099, 0., -0.0896, 1.22756, 0., 0.55135],
             [0.55135, 0., -0.0896, 1.22756, -0., 3.17099]])

        self.run_loprop(inpfile, ref_charges, ref_polarizabilities)

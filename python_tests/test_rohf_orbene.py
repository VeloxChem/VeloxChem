from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.molecularbasis import MolecularBasis


class TestROHFOrbitalEnergies:

    def test_rohf_orbene(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')

        basis_label = 'def2-svp'

        ref_scf_energy = -75.701357737185

        ref_mo_energies = np.array([
            -20.75515824, -1.48908233, -0.87164463, -0.76012964, -0.38050080,
            0.01498698, 0.20467495, 0.72182555, 0.79976285, 1.03362958,
            1.05388501, 1.14447063, 1.22736808, 1.52047675, 1.55255890,
            1.71539586, 2.01836707, 2.48929115, 2.51761915, 3.16493758,
            3.21447072, 3.43655757, 3.74699641, 4.09952729
        ])

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        task.ao_basis = MolecularBasis.read(task.molecule,
                                            basis_label,
                                            ostream=None)

        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            scf_energy = scf_drv.scf_energy
            mo_energies = scf_drv.scf_tensors['E_alpha']
            assert abs(ref_scf_energy - scf_energy) < 1.0e-8
            assert np.max(np.abs(ref_mo_energies - mo_energies)) < 1.0e-4

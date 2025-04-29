from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver


@pytest.mark.solvers
class Test3PASmall:

    def test_3pa_small(self):

        xyzstr = """3
        xyz
        O   0.000000    0.000000    0.000000
        H   0.758602    0.000000   -0.504284
        H   0.758602    0.000000    0.504284
        """

        molecule = Molecule.read_xyz_string(xyzstr)
        basis = MolecularBasis.read(molecule, 'cc-pvdz')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        three_pa_drv = ThreePATransitionDriver()
        three_pa_drv.nstates = 2
        three_pa_drv.ostream.mute()
        three_pa_results = three_pa_drv.compute(molecule, basis, scf_results)

        ref_lin_3pa_str = np.array([1211.4, 3308.2])
        ref_cir_3pa_str = np.array([934.9, 8270.5])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            lin_3pa_str = np.array([
                val.real for key, val in three_pa_results['3pa_strengths']
                ['linear'].items()
            ])
            cir_3pa_str = np.array([
                val.real for key, val in three_pa_results['3pa_strengths']
                ['circular'].items()
            ])
            assert np.max(np.abs(lin_3pa_str - ref_lin_3pa_str)) < 0.1
            assert np.max(np.abs(cir_3pa_str - ref_cir_3pa_str)) < 0.1

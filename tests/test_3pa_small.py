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

    @pytest.mark.parametrize('ri', [False, True])
    def test_3pa_small(self, ri):

        xyzstr = """3
        xyz
        O   0.000000    0.000000    0.000000
        H   0.758602    0.000000   -0.504284
        H   0.758602    0.000000    0.504284
        """

        molecule = Molecule.read_xyz_string(xyzstr)
        basis = MolecularBasis.read(molecule, 'cc-pvdz')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        three_pa_drv = ThreePATransitionDriver()
        three_pa_drv.nstates = 2
        three_pa_drv.ostream.mute()
        three_pa_results = three_pa_drv.compute(molecule, basis, scf_results)

        ref_lin_3pa_str = {
            False: np.array([1211.4, 3308.2]),
            True: np.array([1211.760, 3309.227])
        }[ri]
        ref_cir_3pa_str = {
            False: np.array([934.9, 8270.5]),
            True: np.array([935.158, 8273.067])
        }[ri]

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

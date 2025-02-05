from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.mpitask import MpiTask
from veloxchem.aoindices import get_basis_function_indices_of_atoms
from veloxchem.aoindices import ao_matrix_to_dalton, ao_matrix_to_veloxchem


class TestAOIndices:

    def test_basis_function_indices(self):

        here = Path(__file__).parent
        inpfile = here / 'data' / 'h2se.inp'

        task = MpiTask([str(inpfile), None])
        molecule = task.molecule
        basis = task.ao_basis

        bf_indices = get_basis_function_indices_of_atoms(molecule, basis)

        ref_indices = [
            0, 1, 2, 3, 4, 21, 9, 15, 22, 10, 16, 23, 11, 17, 24, 12, 18, 27,
            30, 33, 36, 39, 28, 31, 34, 37, 40, 29, 32, 35, 38, 41, 5, 6, 25,
            13, 19, 7, 8, 26, 14, 20
        ]

        assert bf_indices == ref_indices

        ovl_drv = OverlapDriver()
        smat = ovl_drv.compute(molecule, basis).to_numpy()

        sdal = ao_matrix_to_dalton(smat, basis, molecule)
        sdal2 = smat[bf_indices, :][:, bf_indices]

        assert np.max(np.abs(sdal - sdal2)) < 1.0e-12

        reverse_indices = [(x, i) for i, x in enumerate(bf_indices)]
        reverse_indices = [x[1] for x in sorted(reverse_indices)]

        svlx = ao_matrix_to_veloxchem(sdal, basis, molecule)
        svlx2 = sdal[reverse_indices, :][:, reverse_indices]

        assert np.max(np.abs(svlx - svlx2)) < 1.0e-12
        assert np.max(np.abs(svlx - smat)) < 1.0e-12

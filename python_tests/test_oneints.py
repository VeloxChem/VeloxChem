from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import mpi_master

import h5py
import math
import numpy as np
import unittest


class TestOneInts(unittest.TestCase):

    def test_overlap_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = OverlapMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_get_ortho_matrix(self):

        arr = np.array([[ 1.0, 0.2, 0.1 ],
                        [ 0.2, 2.0, 0.3 ],
                        [ 0.1, 0.3, 3.0 ]])

        mat = OverlapMatrix(arr)
        ortho_1 = mat.get_ortho_matrix(1.0e-12).to_numpy()

        evals, evecs = np.linalg.eigh(arr)
        evals_sqrt_inv = np.diag([1.0 / math.sqrt(x) for x in evals])
        ortho_2 = np.dot(evecs, np.dot(evals_sqrt_inv, evecs.T))

        diff = np.max(np.abs(ortho_1 - ortho_2))
        self.assertAlmostEqual(0., diff, 13)

    def test_kinetic_energy_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = KineticEnergyMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = KineticEnergyMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_nuclear_potential_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = NuclearPotentialMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = NuclearPotentialMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_1e_integrals(self):

        task = MpiTask(["inputs/h2se.inp", "inputs/h2se.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        basis = task.ao_basis

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver(rank, size, comm)
        S = ovldrv.compute(molecule, basis, comm)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver(rank, size, comm)
        T = kindrv.compute(molecule, basis, comm)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver(rank, size, comm)
        V = npotdrv.compute(molecule, basis, comm)
        V1 = V.to_numpy()

        # compare with reference

        hf = h5py.File("inputs/h2se.onee.h5", 'r')
        S2 = np.array(hf.get("overlap"))
        T2 = np.array(hf.get("kinetic_energy"))
        V2 = np.array(hf.get("nuclear_potential"))
        hf.close()

        if rank == mpi_master():

            dS = np.max(np.abs(S1 - S2))
            dT = np.max(np.abs(T1 - T2))
            dV = np.max(np.abs(V1 - V2))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-11)
            self.assertTrue(dV < 1.0e-11)


if __name__ == "__main__":
    unittest.main()

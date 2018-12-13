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
import numpy as np
import unittest


class TestOneInts(unittest.TestCase):

    def test_overlap_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = OverlapMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

    def test_kinetic_energy_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = KineticEnergyMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = KineticEnergyMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

    def test_nuclear_potential_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = NuclearPotentialMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = NuclearPotentialMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

    def test_1e_integrals(self):

        task = MpiTask(["inputs/h2se.inp", "inputs/h2se.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        basis = task.ao_basis
        ostream = task.ostream

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver.create(rank, size, comm)
        S = ovldrv.compute(molecule, basis, ostream, comm)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver.create(rank, size, comm)
        T = kindrv.compute(molecule, basis, ostream, comm)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver.create(rank, size, comm)
        V = npotdrv.compute(molecule, basis, ostream, comm)
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

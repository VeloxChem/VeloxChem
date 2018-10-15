from mpi4py import MPI
from HelperClass import Task
from veloxchem.VeloxChemLib import OverlapMatrix
from veloxchem.VeloxChemLib import KineticEnergyMatrix
from veloxchem.VeloxChemLib import NuclearPotentialMatrix
from veloxchem.VeloxChemLib import OverlapIntegralsDriver
from veloxchem.VeloxChemLib import KineticEnergyIntegralsDriver
from veloxchem.VeloxChemLib import NuclearPotentialIntegralsDriver
from veloxchem.VeloxChemLib import mpi_master

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

        # mpi settings

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # process input file on master node

        if (rank == mpi_master()):

            task = Task("inputs/h2se.inp", "inputs/h2se.out")
            molecule = task.molecule
            basis = task.ao_basis
            ostream = task.ostream

        else:

            molecule = Molecule()
            basis = MolecularBasis()
            ostream = OutputStream("")

        # broadcast molecule and basis

        molecule.broadcast(rank, comm)
        basis.broadcast(rank, comm)

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

        S_diff = np.max(np.abs(S1 - S2))
        T_diff = np.max(np.abs(T1 - T2))
        V_diff = np.max(np.abs(V1 - V2))

        self.assertTrue(S_diff < 1.0e-13)
        self.assertTrue(T_diff < 1.0e-11)
        self.assertTrue(V_diff < 1.0e-11)


if __name__ == "__main__":
    unittest.main()

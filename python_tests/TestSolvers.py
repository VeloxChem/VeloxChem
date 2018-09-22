from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *

import numpy as np
import unittest

class TestSolvers(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestSolvers " + label)

    def test_sad_guess(self):

        self.print_title("sad_guess")

        # mpi settings

        comm = MPI.COMM_WORLD

        rank = comm.Get_rank()

        size = comm.Get_size()

        # process input file on master node

        if (rank == mpi_master()):

            task = Task("inputs/water.inp", "inputs/water.out")

            molecule = task.molecule

            ao_basis = task.ao_basis

            min_basis = task.min_basis

            ostream = task.ostream

        else:

            molecule = Molecule()

            ao_basis = MolecularBasis()

            min_basis = MolecularBasis()

            ostream = OutputStream("")

        # broadcast molecule and basis

        molecule.broadcast(rank, comm)

        ao_basis.broadcast(rank, comm)

        min_basis.broadcast(rank, comm)

        # compute overlap

        ovldrv = OverlapIntegralsDriver.create(rank, size, comm)

        S12 = ovldrv.compute(molecule, min_basis, ao_basis, ostream, comm)

        S22 = ovldrv.compute(molecule, ao_basis, ostream, comm)

        # compute initial guess

        saddrv = SADGuessDriver.create(rank, size, comm)

        D = saddrv.compute(molecule, min_basis, ao_basis,
                           S12, S22, ostream, comm)

        # matrix to numpy

        overlap = S22.to_numpy()

        density = D.total_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(density.ndim, 2)

            self.assertEqual(density.shape[0], 24)

            self.assertEqual(density.shape[1], 24)

            # number of electrons

            DS = density.dot(overlap)

            nelec = 0.0

            for i in range(DS.shape[0]):

                nelec += DS[i][i]

            nelec *= 2.0

            self.assertAlmostEqual(nelec, 10., 13)

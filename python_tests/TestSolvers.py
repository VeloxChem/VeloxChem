from mpi4py import MPI
from HelperClass import Task
from veloxchem.VeloxChemLib import Molecule
from veloxchem.VeloxChemLib import MolecularBasis
from veloxchem.VeloxChemLib import OutputStream
from veloxchem.VeloxChemLib import OverlapIntegralsDriver
from veloxchem.VeloxChemLib import SADGuessDriver
from veloxchem.VeloxChemLib import mpi_master

import numpy as np
import unittest


class TestSolvers(unittest.TestCase):

    def test_sad_guess(self):

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
        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        # matrix to numpy

        overlap = S22.to_numpy()
        density = D.total_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(density.ndim, 2)
            self.assertEqual(density.shape[0], 24)
            self.assertEqual(density.shape[1], 24)

            # number of electrons

            self.assertAlmostEqual(2.0 * np.sum(density * overlap), 10., 13)

        # compute other initial guess

        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # closed-shell cation initial guess

        molecule.set_charge(charge + 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        density_a = D.total_to_numpy(0)

        if (rank == mpi_master()):

            self.assertAlmostEqual(np.sum(density_a * overlap), 4., 13)

        # closed-shell anion initial guess

        molecule.set_charge(charge - 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        density_a = D.total_to_numpy(0)

        if (rank == mpi_master()):

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)

        # open-shell cation initial guess

        molecule.set_charge(charge + 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertAlmostEqual(np.sum(density_a * overlap), 5., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 4., 13)

        # open-shell anion initial guess

        molecule.set_charge(charge - 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 5., 13)

        # open-shell triplet initial guess

        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity + 2)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                           comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 4., 13)


if __name__ == "__main__":
    unittest.main()

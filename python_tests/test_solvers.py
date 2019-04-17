from mpi4py import MPI
import numpy as np
import unittest

from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import SADGuessDriver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask


class TestSolvers(unittest.TestCase):

    def test_sad_guess(self):

        task = MpiTask(["inputs/water.inp", "inputs/water.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # compute overlap

        ovldrv = OverlapIntegralsDriver(rank, size, comm)
        S12 = ovldrv.compute(molecule, min_basis, ao_basis, comm)
        S22 = ovldrv.compute(molecule, ao_basis, comm)

        # compute initial guess

        saddrv = SADGuessDriver(rank, size, comm)
        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        # matrix to numpy

        overlap = S22.to_numpy()
        density = D.alpha_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(density.ndim, 2)
            self.assertEqual(density.shape[0], 24)
            self.assertEqual(density.shape[1], 24)

            # number of electrons

            self.assertEqual(molecule.number_of_electrons(), 10)
            self.assertEqual(molecule.number_of_alpha_electrons(), 5)
            self.assertEqual(molecule.number_of_beta_electrons(), 5)

            self.assertAlmostEqual(2.0 * np.sum(density * overlap), 10., 13)

        # compute other initial guess

        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # closed-shell cation initial guess

        molecule.set_charge(charge + 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        density_a = D.alpha_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(molecule.number_of_electrons(), 8)
            self.assertEqual(molecule.number_of_alpha_electrons(), 4)
            self.assertEqual(molecule.number_of_beta_electrons(), 4)

            self.assertAlmostEqual(np.sum(density_a * overlap), 4., 13)

        # closed-shell anion initial guess

        molecule.set_charge(charge - 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        density_a = D.alpha_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(molecule.number_of_electrons(), 12)
            self.assertEqual(molecule.number_of_alpha_electrons(), 6)
            self.assertEqual(molecule.number_of_beta_electrons(), 6)

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)

        # open-shell cation initial guess

        molecule.set_charge(charge + 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(molecule.number_of_electrons(), 9)
            self.assertEqual(molecule.number_of_alpha_electrons(), 5)
            self.assertEqual(molecule.number_of_beta_electrons(), 4)

            self.assertAlmostEqual(np.sum(density_a * overlap), 5., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 4., 13)

        # open-shell anion initial guess

        molecule.set_charge(charge - 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(molecule.number_of_electrons(), 11)
            self.assertEqual(molecule.number_of_alpha_electrons(), 6)
            self.assertEqual(molecule.number_of_beta_electrons(), 5)

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 5., 13)

        # open-shell triplet initial guess

        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity + 2)

        D = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, comm)

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if (rank == mpi_master()):

            self.assertEqual(molecule.number_of_electrons(), 10)
            self.assertEqual(molecule.number_of_alpha_electrons(), 6)
            self.assertEqual(molecule.number_of_beta_electrons(), 4)

            self.assertAlmostEqual(np.sum(density_a * overlap), 6., 13)
            self.assertAlmostEqual(np.sum(density_b * overlap), 4., 13)


if __name__ == "__main__":
    unittest.main()

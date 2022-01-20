from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import SADGuessDriver
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask


class TestInitialGuess:

    def test_sad_guess(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis

        # compute overlap

        ovldrv = OverlapIntegralsDriver(task.mpi_comm)
        S = ovldrv.compute(molecule, ao_basis)

        # compute initial guess

        saddrv = SADGuessDriver(task.mpi_comm)
        D = saddrv.compute(molecule, min_basis, ao_basis, 'restricted')

        # matrix to numpy

        overlap = S.to_numpy()
        density = D.alpha_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert density.ndim == 2
            assert density.shape[0] == 41
            assert density.shape[1] == 41

            # number of electrons

            assert molecule.number_of_electrons() == 10
            assert molecule.number_of_alpha_electrons() == 5
            assert molecule.number_of_beta_electrons() == 5

            assert abs(2.0 * np.sum(density * overlap) - 10.) < 1.0e-13

        # compute other initial guess

        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # closed-shell cation initial guess

        molecule.set_charge(charge + 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, 'restricted')

        density_a = D.alpha_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert molecule.number_of_electrons() == 8
            assert molecule.number_of_alpha_electrons() == 4
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density_a * overlap) - 4.) < 1.0e-13

        # closed-shell anion initial guess

        molecule.set_charge(charge - 2)
        molecule.set_multiplicity(multiplicity)

        D = saddrv.compute(molecule, min_basis, ao_basis, 'restricted')

        density_a = D.alpha_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert molecule.number_of_electrons() == 12
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 6

            assert abs(np.sum(density_a * overlap) - 6.) < 1.0e-13

        # open-shell cation initial guess

        molecule.set_charge(charge + 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, 'unrestricted')

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert molecule.number_of_electrons() == 9
            assert molecule.number_of_alpha_electrons() == 5
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density_a * overlap) - 5.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 4.) < 1.0e-13

        # open-shell anion initial guess

        molecule.set_charge(charge - 1)
        molecule.set_multiplicity(multiplicity + 1)

        D = saddrv.compute(molecule, min_basis, ao_basis, 'unrestricted')

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert molecule.number_of_electrons() == 11
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 5

            assert abs(np.sum(density_a * overlap) - 6.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 5.) < 1.0e-13

        # open-shell triplet initial guess

        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity + 2)

        D = saddrv.compute(molecule, min_basis, ao_basis, 'unrestricted')

        density_a = D.alpha_to_numpy(0)
        density_b = D.beta_to_numpy(0)

        if is_mpi_master(task.mpi_comm):

            assert molecule.number_of_electrons() == 10
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density_a * overlap) - 6.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 4.) < 1.0e-13

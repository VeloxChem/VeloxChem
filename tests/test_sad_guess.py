from pathlib import Path
import pytest
import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.sadguessdriver import SadGuessDriver
from veloxchem.mpitask import MpiTask


@pytest.mark.solvers
class TestSadGuess:

    def test_sad_guess(self):

        here = Path(__file__).parent
        inpfile = here / 'data' / 'water.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis

        # compute overlap

        ovl_drv = OverlapDriver()
        S = ovl_drv.compute(molecule, ao_basis)

        # compute initial guess

        sad_drv = SadGuessDriver()
        density = sad_drv.compute(molecule, min_basis, ao_basis,
                                  'restricted')[0]

        # matrix to numpy

        if task.mpi_rank == mpi_master():

            overlap = S.to_numpy()

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

        density = sad_drv.compute(molecule, min_basis, ao_basis,
                                  'restricted')[0]

        if task.mpi_rank == mpi_master():

            assert molecule.number_of_electrons() == 8
            assert molecule.number_of_alpha_electrons() == 4
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density * overlap) - 4.) < 1.0e-13

        # closed-shell anion initial guess

        molecule.set_charge(charge - 2)
        molecule.set_multiplicity(multiplicity)

        density = sad_drv.compute(molecule, min_basis, ao_basis,
                                  'restricted')[0]

        if task.mpi_rank == mpi_master():

            assert molecule.number_of_electrons() == 12
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 6

            assert abs(np.sum(density * overlap) - 6.) < 1.0e-13

        # open-shell cation initial guess

        molecule.set_charge(charge + 1)
        molecule.set_multiplicity(multiplicity + 1)

        density_a, density_b = sad_drv.compute(molecule, min_basis, ao_basis,
                                               'unrestricted')

        if task.mpi_rank == mpi_master():

            assert molecule.number_of_electrons() == 9
            assert molecule.number_of_alpha_electrons() == 5
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density_a * overlap) - 5.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 4.) < 1.0e-13

        # open-shell anion initial guess

        molecule.set_charge(charge - 1)
        molecule.set_multiplicity(multiplicity + 1)

        density_a, density_b = sad_drv.compute(molecule, min_basis, ao_basis,
                                               'unrestricted')

        if task.mpi_rank == mpi_master():

            assert molecule.number_of_electrons() == 11
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 5

            assert abs(np.sum(density_a * overlap) - 6.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 5.) < 1.0e-13

        # open-shell triplet initial guess

        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity + 2)

        density_a, density_b = sad_drv.compute(molecule, min_basis, ao_basis,
                                               'unrestricted')

        if task.mpi_rank == mpi_master():

            assert molecule.number_of_electrons() == 10
            assert molecule.number_of_alpha_electrons() == 6
            assert molecule.number_of_beta_electrons() == 4

            assert abs(np.sum(density_a * overlap) - 6.) < 1.0e-13
            assert abs(np.sum(density_b * overlap) - 4.) < 1.0e-13

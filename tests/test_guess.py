from mpi4py import MPI
import numpy as np
import os

from veloxchem.veloxchemlib import OverlapDriver, GpuDevices
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.sadguessdriver import SadGuessDriver


class TestInitialGuess:

    def test_sad_guess(self):

        if MPI.COMM_WORLD.Get_rank() == 0:

            molstr = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
            """
            mol = Molecule.read_molecule_string(molstr, 'au')

            bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
            min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

            ovl_drv = OverlapDriver()
            ovl_mat = ovl_drv.compute(mol, bas)
            S = ovl_mat.get_full_matrix().to_numpy()

            if 'VLX_NUM_GPUS_PER_NODE' in os.environ:
                num_gpus_per_node = int(os.environ['VLX_NUM_GPUS_PER_NODE'])
            else:
                devices = GpuDevices()
                num_gpus_per_node = devices.get_number_devices()

            saddrv = SadGuessDriver()
            D = saddrv.compute(mol, min_bas, bas, 'restricted', num_gpus_per_node)

            assert D.ndim == 2
            assert D.shape[0] == 24
            assert D.shape[1] == 24

            assert mol.number_of_electrons() == 10
            assert mol.number_of_alpha_electrons() == 5
            assert mol.number_of_beta_electrons() == 5

            assert abs(np.sum(D * S) - 5.) < 1.0e-13

            # compute other initial guess

            charge = mol.get_charge()
            multiplicity = mol.get_multiplicity()

            # closed-shell cation initial guess

            mol.set_charge(charge + 2)
            mol.set_multiplicity(multiplicity)

            D = saddrv.compute(mol, min_bas, bas, 'restricted', num_gpus_per_node)

            assert mol.number_of_electrons() == 8
            assert mol.number_of_alpha_electrons() == 4
            assert mol.number_of_beta_electrons() == 4

            assert abs(np.sum(D * S) - 4.) < 1.0e-13

            # closed-shell anion initial guess

            mol.set_charge(charge - 2)
            mol.set_multiplicity(multiplicity)

            D = saddrv.compute(mol, min_bas, bas, 'restricted', num_gpus_per_node)

            assert mol.number_of_electrons() == 12
            assert mol.number_of_alpha_electrons() == 6
            assert mol.number_of_beta_electrons() == 6

            assert abs(np.sum(D * S) - 6.) < 1.0e-13

            # open-shell cation initial guess

            mol.set_charge(charge + 1)
            mol.set_multiplicity(multiplicity + 1)

            Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted', num_gpus_per_node)

            assert mol.number_of_electrons() == 9
            assert mol.number_of_alpha_electrons() == 5
            assert mol.number_of_beta_electrons() == 4

            assert abs(np.sum(Da * S) - 5.) < 1.0e-13
            assert abs(np.sum(Db * S) - 4.) < 1.0e-13

            # open-shell anion initial guess

            mol.set_charge(charge - 1)
            mol.set_multiplicity(multiplicity + 1)

            Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted', num_gpus_per_node)

            assert mol.number_of_electrons() == 11
            assert mol.number_of_alpha_electrons() == 6
            assert mol.number_of_beta_electrons() == 5

            assert abs(np.sum(Da * S) - 6.) < 1.0e-13
            assert abs(np.sum(Db * S) - 5.) < 1.0e-13

            # open-shell triplet initial guess

            mol.set_charge(charge)
            mol.set_multiplicity(multiplicity + 2)

            Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted', num_gpus_per_node)

            assert mol.number_of_electrons() == 10
            assert mol.number_of_alpha_electrons() == 6
            assert mol.number_of_beta_electrons() == 4

            assert abs(np.sum(Da * S) - 6.) < 1.0e-13
            assert abs(np.sum(Db * S) - 4.) < 1.0e-13

        MPI.COMM_WORLD.barrier()

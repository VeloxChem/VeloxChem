from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.sadguessdriver import SadGuessDriver


class TestInitialGuess:

    def test_sad_guess(self):

        molstr = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(molstr, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bas = MolecularBasis.read(mol, 'def2-svp', basis_path, ostream=None)
        min_bas = MolecularBasis.read(mol,
                                      'ao-start-guess',
                                      basis_path,
                                      ostream=None)

        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol, bas)
        S = ovl_mat.get_full_matrix().to_numpy()

        saddrv = SadGuessDriver()
        D = saddrv.compute(mol, min_bas, bas, 'restricted')

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

        D = saddrv.compute(mol, min_bas, bas, 'restricted')

        assert mol.number_of_electrons() == 8
        assert mol.number_of_alpha_electrons() == 4
        assert mol.number_of_beta_electrons() == 4

        assert abs(np.sum(D * S) - 4.) < 1.0e-13

        # closed-shell anion initial guess

        mol.set_charge(charge - 2)
        mol.set_multiplicity(multiplicity)

        D = saddrv.compute(mol, min_bas, bas, 'restricted')

        assert mol.number_of_electrons() == 12
        assert mol.number_of_alpha_electrons() == 6
        assert mol.number_of_beta_electrons() == 6

        assert abs(np.sum(D * S) - 6.) < 1.0e-13

        # open-shell cation initial guess

        mol.set_charge(charge + 1)
        mol.set_multiplicity(multiplicity + 1)

        Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted')

        assert mol.number_of_electrons() == 9
        assert mol.number_of_alpha_electrons() == 5
        assert mol.number_of_beta_electrons() == 4

        assert abs(np.sum(Da * S) - 5.) < 1.0e-13
        assert abs(np.sum(Db * S) - 4.) < 1.0e-13

        # open-shell anion initial guess

        mol.set_charge(charge - 1)
        mol.set_multiplicity(multiplicity + 1)

        Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted')

        assert mol.number_of_electrons() == 11
        assert mol.number_of_alpha_electrons() == 6
        assert mol.number_of_beta_electrons() == 5

        assert abs(np.sum(Da * S) - 6.) < 1.0e-13
        assert abs(np.sum(Db * S) - 5.) < 1.0e-13

        # open-shell triplet initial guess

        mol.set_charge(charge)
        mol.set_multiplicity(multiplicity + 2)

        Da, Db = saddrv.compute(mol, min_bas, bas, 'unrestricted')

        assert mol.number_of_electrons() == 10
        assert mol.number_of_alpha_electrons() == 6
        assert mol.number_of_beta_electrons() == 4

        assert abs(np.sum(Da * S) - 6.) < 1.0e-13
        assert abs(np.sum(Db * S) - 4.) < 1.0e-13

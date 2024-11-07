from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import NuclearPotentialErfDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix


class TestNuclearPotentialErfDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def test_nuclear_potential_erf_co_qzvp(self):

        mol_co, bas_qzvp = self.get_data()

        # external charges, coordinates
        charges = [
            1.0,
            2.0,
        ]
        coords = [
            [0.1, 0.2, 0.3],
            [1.0, 1.2, 1.4],
        ]

        # compute range separated nuclear potential matrix
        npot_drv = NuclearPotentialErfDriver()
        npot_mat = npot_drv.compute(mol_co, bas_qzvp, charges, coords, 0.64)

    # TODO: Need test data for range-separated nuclear potential

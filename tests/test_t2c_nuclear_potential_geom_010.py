from pathlib import Path
import numpy as np

from veloxchem import NuclearPotentialGeom010Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestNuclearPotentialGeom010Driver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def test_nuclear_potential_geom_010_co_qzvp(self):

        mol_co, bas_qzvp = self.get_data()
        
        # external dipoles, coordinates
        dipoles = [1.0, 2.0, 3.0, 1.5, 2.5, 3.5,]
        coords = [[0.1, 0.2, 0.3], [1.0, 1.2, 1.4], ]
        
        # compute nuclear potential matrix
        geom_drv = NuclearPotentialGeom010Driver()
        geom_mat = geom_drv.compute(mol_co, bas_qzvp, dipoles, coords)

       # TODO: Need test data for electric field integrals

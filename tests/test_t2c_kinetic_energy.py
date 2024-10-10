from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import KineticEnergyDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix


class TestKineticEnergyDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def test_kinetic_energy_co_qzvp(self):

        mol, bas = self.get_data()

        # compute kinetic energy matrix
        kin_drv = KineticEnergyDriver()
        kin_mat = kin_drv.compute(mol, bas)

        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.kinetic.energy.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # check individual kinetic energy submatrices
        for i, j in zip(indexes[0], indexes[1]):
            sbra = basdims[i]
            ebra = basdims[i + 1]
            sket = basdims[j]
            eket = basdims[j + 1]
            cmat = kin_mat.submatrix((i, j))
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            assert cmat == rmat

        # check full kinetic energy matrix
        fmat = kin_mat.full_matrix()
        fref = SubMatrix([0, 0, 114, 114])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

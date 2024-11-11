from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import NuclearPotentialDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix


class TestNuclearPotentialDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def test_nuclear_potential_co_qzvp(self):

        mol_co, bas_qzvp = self.get_data()

        # compute nuclear potential matrix
        npot_drv = NuclearPotentialDriver()
        npot_mat = npot_drv.compute(mol_co, bas_qzvp)

        # load reference nuclear potential data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.nuclear.potential.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # check individual nuclear potential submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = npot_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full nuclear potential matrix
        fmat = npot_mat.full_matrix()
        fref = SubMatrix([0, 0, 114, 114])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

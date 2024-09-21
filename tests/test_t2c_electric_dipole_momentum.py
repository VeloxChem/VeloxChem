from pathlib import Path
import numpy as np

from veloxchem import ElectricDipoleMomentDriver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestElectricDipoleMomentDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP')

        return mol, bas

    def test_electric_dipole_moment_co_qzvp(self):

        mol, bas = self.get_data()

        # compute overlap matrix
        dip_drv = ElectricDipoleMomentDriver()
        dip_mats = dip_drv.compute(mol, bas, [0.0, 0.0, 0.0])

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.electric.dipole.momentum.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = dip_mats.matrix(label)
            for i, j in zip(indexes[0], indexes[1]):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fmat.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 114, 114])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref

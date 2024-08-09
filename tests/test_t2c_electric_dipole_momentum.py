from pathlib import Path
import numpy as np

from veloxchem import ElectricDipoleMomentumDriver
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
        dip_drv = ElectricDipoleMomentumDriver()
        dip_mats = dip_drv.compute(mol, bas, [0.0, 0.0, 0.0])

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.electric.dipole.momentum.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        fmat_x = dip_mats.matrix('X')
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fmat_x.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(
                np.ascontiguousarray(ref_mat[0, sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        fmat_y = dip_mats.matrix('Y')
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fmat_y.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(
                np.ascontiguousarray(ref_mat[1, sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        fmat_z = dip_mats.matrix('Z')
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fmat_z.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(
                np.ascontiguousarray(ref_mat[2, sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        fmat_x = dip_mats.matrix('X').full_matrix()
        fref_x = SubMatrix([0, 0, 114, 114])
        fref_x.set_values(np.ascontiguousarray(ref_mat[0]))
        print(fmat_x.to_numpy())
        print(fref_x.to_numpy())
        assert fmat_x == fref_x

        fmat_y = dip_mats.matrix('Y').full_matrix()
        fref_y = SubMatrix([0, 0, 114, 114])
        fref_y.set_values(np.ascontiguousarray(ref_mat[1]))
        assert fmat_y == fref_y

        fmat_z = dip_mats.matrix('Z').full_matrix()
        fref_z = SubMatrix([0, 0, 114, 114])
        fref_z.set_values(np.ascontiguousarray(ref_mat[2]))
        assert fmat_z == fref_z

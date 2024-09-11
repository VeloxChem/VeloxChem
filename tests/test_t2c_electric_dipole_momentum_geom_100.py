from pathlib import Path
import numpy as np

from veloxchem import ElectricDipoleMomentumGeom100Driver
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

    def test_electric_dipole_moment_co_qzvp_for_c(self):

        mol, bas = self.get_data()

        # compute electric dipole derivatives matrix
        dip_drv = ElectricDipoleMomentumGeom100Driver()
        dip_mats = dip_drv.compute(mol, bas, [0.0, 0.0, 0.0], 0)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.electric.dipole.moment.geom.100.c.npy')
        ref_mat = np.load(npyfile)
        ref_mat = ref_mat.transpose(0, 2, 1)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'Y_X', 'Z_X', 'X_Y', 'Y_Y', 'Z_Y', 'X_Z', 'Y_Z', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = dip_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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

    def test_electric_dipole_moment_co_qzvp_for_o(self):

        mol, bas = self.get_data()

        # compute electric dipole derivatives matrix
        dip_drv = ElectricDipoleMomentumGeom100Driver()
        dip_mats = dip_drv.compute(mol, bas, [0.0, 0.0, 0.0], 1)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.electric.dipole.moment.geom.100.o.npy')
        ref_mat = np.load(npyfile)
        ref_mat = ref_mat.transpose(0, 2, 1)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'Y_X', 'Z_X', 'X_Y', 'Y_Y', 'Z_Y', 'X_Z', 'Y_Z', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = dip_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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

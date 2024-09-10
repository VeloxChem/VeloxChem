from pathlib import Path
import numpy as np

from veloxchem import OverlapGeom200Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestOverlapGeom200Driver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP')

        return mol, bas

    def test_overlap_co_qzvp_for_cc(self):

        mol, bas = self.get_data()

        # compute overlap hessian matrix
        hess_drv = OverlapGeom200Driver()
        hess_mats = hess_drv.compute(mol, bas, 0)

        # load reference overlap hessian for C,C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.overlap.geom.200.cc.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
        matids = [0, 1, 2, 4, 5, 8]

        for k, label in zip(matids, labels):
            fmat = hess_mats.matrix(label)
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

    def test_hessian_co_qzvp_for_oo(self):

        mol, bas = self.get_data()

        # compute overlap gradient matrix
        hess_drv = OverlapGeom200Driver()
        hess_mats = hess_drv.compute(mol, bas, 1)

        # load reference overlap hessian for O,O atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.overlap.geom.200.oo.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
        matids = [0, 1, 2, 4, 5, 8]

        for k, label in zip(matids, labels):
            fmat = hess_mats.matrix(label)
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

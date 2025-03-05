from pathlib import Path
import numpy as np

from veloxchem import TwoCenterElectronRepulsionGeom100Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestTwoCenterElectronRepulsionGeom100Driver:

    def get_data(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def test_electron_repulsion_h2o_svp_for_o1(self):

        mol, bas = self.get_data()

        # compute electron repulsion gradient matrix
        grad_drv = TwoCenterElectronRepulsionGeom100Driver()
        grad_mats = grad_drv.compute(mol, bas, 0)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int2c2e.geom.100.o1.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 7, 19, 24]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = grad_mats.matrix(label)
            for i in range(0, 3):
                for j in range(0, 3):
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
            fref = SubMatrix([0, 0, 24, 24])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
            
    def test_electron_repulsion_h2o_svp_for_h2(self):

        mol, bas = self.get_data()

        # compute electron repulsion gradient matrix
        grad_drv = TwoCenterElectronRepulsionGeom100Driver()
        grad_mats = grad_drv.compute(mol, bas, 1)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int2c2e.geom.100.h2.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 7, 19, 24]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = grad_mats.matrix(label)
            for i in range(0, 3):
                for j in range(0, 3):
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
            fref = SubMatrix([0, 0, 24, 24])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
            
    def test_electron_repulsion_h2o_svp_for_h3(self):

        mol, bas = self.get_data()

        # compute electron repulsion gradient matrix
        grad_drv = TwoCenterElectronRepulsionGeom100Driver()
        grad_mats = grad_drv.compute(mol, bas, 2)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int2c2e.geom.100.h3.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 7, 19, 24]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = grad_mats.matrix(label)
            for i in range(0, 3):
                for j in range(0, 3):
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
            fref = SubMatrix([0, 0, 24, 24])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref

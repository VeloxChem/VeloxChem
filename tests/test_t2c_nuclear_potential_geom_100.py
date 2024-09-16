from pathlib import Path
import numpy as np

from veloxchem import NuclearPotentialGeom100Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestOverlapGeom100Driver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP')

        return mol, bas

    def test_nuclear_potential_co_qzvp_for_c(self):

        mol, bas = self.get_data()

        # compute nuclear potential gradient matrix
        grad_drv = NuclearPotentialGeom100Driver()
        grad_mats = grad_drv.compute(mol, bas, 0)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.nuclear.potential.geom.100.c.npy')
        ref_mat = np.load(npyfile)
        ref_mat = -ref_mat

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = grad_mats.matrix(label)
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

    def test_nuclear_potential_co_qzvp_for_o(self):

        mol, bas = self.get_data()

        # compute nuclear potential gradient matrix
        grad_drv = NuclearPotentialGeom100Driver()
        grad_mats = grad_drv.compute(mol, bas, 1)

        # load reference overlap gradient for O atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.nuclear.potential.geom.100.o.npy')
        ref_mat = np.load(npyfile)
        ref_mat = -ref_mat

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = grad_mats.matrix(label)
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
                    # assert cmat == rmat
                    assert np.max(np.abs(cmat.to_numpy() -
                                         rmat.to_numpy())) < 1.0e-11
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 114, 114])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            # assert smat == fref
            assert np.max(np.abs(smat.to_numpy() - fref.to_numpy())) < 1.0e-11

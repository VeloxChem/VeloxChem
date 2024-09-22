from pathlib import Path
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockGeom1000Driver
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t


class TestFockDriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas

    def test_h2o_fock_2jk_grad_h3_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "2jk", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h3.npy')
        ref_mat = 2.0 * np.load(npyfile)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h3.npy')
        ref_mat = ref_mat - np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_2jk_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 1, "2jk", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h2.npy')
        ref_mat = 2.0 * np.load(npyfile)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h2.npy')
        ref_mat = ref_mat - np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_2jkx_grad_h3_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "2jkx", 0.64, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h3.npy')
        ref_mat = 2.0 * np.load(npyfile)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h3.npy')
        ref_mat = ref_mat - 0.64 * np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_2jkx_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 1, "2jkx", 0.21, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h2.npy')
        ref_mat = 2.0 * np.load(npyfile)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h2.npy')
        ref_mat = ref_mat - 0.21 * np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_j_grad_h3_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "j", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h3.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_j_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 1, "j", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.geom.1000.h2.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_k_grad_h3_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "k", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h3.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_k_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 1, "k", 0.0, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h2.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_kx_grad_h3_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "kx", 0.38, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h3.npy')
        ref_mat = 0.38 * np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        
    def test_h2o_fock_k_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockGeom1000Driver()
        fock_mats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 1, "kx", 0.21, 0.0)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.geom.1000.h2.npy')
        ref_mat = 0.21 * np.load(npyfile)
        
        # dimension of molecular basis
        basdims = [0, 4, 7]
        
        # check individual submatrices of X matrix
        fock_mat_x = fock_mats.matrix("X")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_x.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[0][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_x.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[0]))
        
        # check individual submatrices of Y matrix
        fock_mat_y = fock_mats.matrix("Y")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_y.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[1][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_y.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[1]))
        
        # check individual submatrices of Z matrix
        fock_mat_z = fock_mats.matrix("Z")
        for i in range(2):
            for j in range(2):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fock_mat_z.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[2][sbra:ebra,
                                                                sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat_z.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat[2]))
        

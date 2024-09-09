from pathlib import Path
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockDriver
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t

class TestFockDriver:

    def get_data_h4(self):

        h4str = """
            H 1.0 0.8  0.9
            H 2.1 2.0  1.7
            H 0.1 0.3 -0.3
            H 2.1 3.1  2.0
        """
        mol = Molecule.read_str(h4str, 'au')
        bas = MolecularBasis.read(mol, 'def2-sv(p)')

        return mol, bas
        
    def get_data_co(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP')

        return mol, bas
        
    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas

    def test_h4_fock_2jk_svp(self):

        mol_h4, bas_svp = self.get_data_h4()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h4.sv(p).density.npy')
        den_mat = make_matrix(bas_svp, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svp, mol_h4, den_mat, "2jk", 0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h4.sv(p).fock.2jk.npy')
        ref_mat = np.load(npyfile)

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 8, 8])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref

    def test_co_fock_2jk_qzvp(self):

        mol, bas = self.get_data_co()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.density.npy')
        den_mat = make_matrix(bas, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))
        
        ## compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas, mol, den_mat, "2jk", 0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.fock.2j-k.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 114, 114])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref
        
    def test_h2o_fock_2jk_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "2jk", 0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.2jk.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(2)
        basdims = [0, 4, 7]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref
        
    def test_h2o_fock_2jkx_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "2jkx", 0.68, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.npy')
        ref_mat = 2.0 * np.load(npyfile)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.npy')
        ref_mat -= 0.68 * np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(2)
        basdims = [0, 4, 7]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref
        
    def test_h2o_fock_j_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "j", 0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(2)
        basdims = [0, 4, 7]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref
        
    def test_h2o_fock_k_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "k", 0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(2)
        basdims = [0, 4, 7]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref
        
    def test_h2o_fock_kx_sto3g(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "kx", 0.68, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.npy')
        ref_mat = 0.68 * np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(2)
        basdims = [0, 4, 7]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = fock_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))
        
        assert fmat == fref

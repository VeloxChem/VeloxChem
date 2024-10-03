from pathlib import Path
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockDriver
from veloxchem import T4CScreener
from veloxchem import SubMatrix
from veloxchem import make_matrix
from veloxchem import mat_t


class TestFockErfDriver:

    def get_data_h2o_dimer(self):

        h2o_dimer_str = """
            O -1.464  0.099  0.300
            H -1.956  0.624 -0.340
            H -1.797 -0.799  0.206
            O  1.369  0.146 -0.395
            H  1.894  0.486  0.335
            H  0.451  0.163 -0.083
        """
        mol = Molecule.read_str(h2o_dimer_str, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-tzvp')

        return mol, bas

    def test_h2o_dimer_fock_k_tzvp_with_screener(self):

        mol_h2o_dimer, bas_tzvp = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.density.npy')
        den_mat = make_matrix(bas_tzvp, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_tzvp, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute_full_fock_serial(t4c_drv, den_mat, "k_rs",
                                                     0.0, 0.63, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.erf.k.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(4)
        basdims = [0, 22, 52, 72, 86]

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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 86, 86])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_kx_tzvp_with_screener(self):

        mol_h2o_dimer, bas_tzvp = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.density.npy')
        den_mat = make_matrix(bas_tzvp, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_tzvp, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute_full_fock_serial(t4c_drv, den_mat, "kx_rs",
                                                     0.21, 0.63, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.erf.k.npy')
        ref_mat = 0.21 * np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(4)
        basdims = [0, 22, 52, 72, 86]

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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 86, 86])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

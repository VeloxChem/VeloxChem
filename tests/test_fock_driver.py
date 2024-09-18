from pathlib import Path
import numpy as np

from mpi4py import MPI

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockDriver
from veloxchem import T4CScreener
from veloxchem import SubMatrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t


class TestFockDriver:

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
        bas = MolecularBasis.read(mol, 'def2-svpd')

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

    def test_h2o_dimer_focks_svpd_with_mpi(self):

        comm = MPI.COMM_WORLD

        t4c_drv = None
        den_mats = None

        if comm.Get_rank() == 0:
            mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

            # load density matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
            dmat = make_matrix(bas_svpd, mat_t.general)
            dmat.set_values(np.load(npyfile))
            den_mats = Matrices()
            den_mats.add(dmat, "0")
            den_mats.add(dmat, "1")

            # screen basis function pairs
            t4c_drv = T4CScreener()
            t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # broadcast data
        t4c_drv = comm.bcast(t4c_drv, 0)
        den_mats = Matrices.bcast(den_mats, comm, 0)

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mats = fock_drv.mpi_compute(comm, t4c_drv, den_mats, ["j", "k"],
                                         0.0, 0.0, 15)

        if comm.Get_rank() == 0:
            # dimension of molecular basis
            basdims = [0, 16, 58, 78]

            # load reference Fock matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.gen.npy')
            ref_mat = np.load(npyfile)

            # retrieve first Fock matrix
            fock_mat = fock_mats.matrix("0")

            # check individual overlap submatrices
            for i in range(3):
                for j in range(3):
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
                    rmat.set_values(
                        np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                    # compare submatrices
                    assert cmat == rmat

            # check full Fock matrix
            fmat = fock_mat.full_matrix()
            fref = SubMatrix([0, 0, 78, 78])
            fref.set_values(np.ascontiguousarray(ref_mat))

            # load reference Fock matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
            ref_mat = np.load(npyfile)
            ref_mat = ref_mat.T

            # retrieve second Fock matrix
            fock_mat = fock_mats.matrix("1")

            # check individual overlap submatrices
            for i in range(3):
                for j in range(3):
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
                    rmat.set_values(
                        np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                    # compare submatrices
                    assert cmat == rmat

            # check full Fock matrix
            fmat = fock_mat.full_matrix()
            fref = SubMatrix([0, 0, 78, 78])
            fref.set_values(np.ascontiguousarray(ref_mat))

            assert fmat == fref

    def test_h2o_dimer_fock_2jk_svpd_with_mpi(self):

        comm = MPI.COMM_WORLD

        t4c_drv = None
        den_mat = None
        if comm.Get_rank() == 0:
            mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()
            # load density matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
            den_mat = make_matrix(bas_svpd, mat_t.symmetric)
            den_mat.set_values(np.load(npyfile))
            # screen basis function pairs
            t4c_drv = T4CScreener()
            t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # broadcast data
        t4c_drv = comm.bcast(t4c_drv, 0)
        den_mat = comm.bcast(den_mat, 0)

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.mpi_compute(comm, t4c_drv, den_mat, "2jk", 0.0, 0.0,
                                        15)

        if comm.Get_rank() == 0:
            # load reference Fock matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
            ref_mat = 2.0 * np.load(npyfile)

            # load reference Fock matrix
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
            ref_mat -= np.load(npyfile)

            # dimension of molecular basis
            indexes = np.triu_indices(3)
            basdims = [0, 16, 58, 78]

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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

            # check full Fock matrix
            fmat = fock_mat.full_matrix()
            fref = SubMatrix([0, 0, 78, 78])
            fref.set_values(np.ascontiguousarray(ref_mat))

            assert fmat == fref

    def test_h2o_dimer_focks_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        dmat = make_matrix(bas_svpd, mat_t.general)
        dmat.set_values(np.load(npyfile))

        den_mats = Matrices()
        den_mats.add(dmat, "0")
        den_mats.add(dmat, "1")

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mats = fock_drv.compute(t4c_drv, den_mats, ["j", "k"], 0.0, 0.0,
                                     15)

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.gen.npy')
        ref_mat = np.load(npyfile)

        # retrieve first Fock matrix
        fock_mat = fock_mats.matrix("0")

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
        ref_mat = np.load(npyfile)
        ref_mat = ref_mat.T

        # retrieve second Fock matrix
        fock_mat = fock_mats.matrix("1")

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_gen_2jk_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        den_mat = make_matrix(bas_svpd, mat_t.general)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "2jk", 0.0, 0.0, 15)
        
        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.gen.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
        ref_kmat = np.load(npyfile)
        ref_mat = ref_mat - ref_kmat.T

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_gen_2jkx_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        den_mat = make_matrix(bas_svpd, mat_t.general)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "2jkx", 0.28, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.gen.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
        ref_kmat = np.load(npyfile)
        ref_mat = ref_mat - 0.28 * ref_kmat.T

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_gen_j_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        den_mat = make_matrix(bas_svpd, mat_t.general)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "j", 0.0, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.gen.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_gen_k_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        den_mat = make_matrix(bas_svpd, mat_t.general)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "k", 0.0, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
        ref_mat = np.load(npyfile)
        ref_mat = ref_mat.T

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_gen_kx_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy')
        den_mat = make_matrix(bas_svpd, mat_t.general)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "kx", 0.73, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.gen.npy')
        ref_mat = np.load(npyfile)
        ref_mat = 0.73 * ref_mat.T

        # dimension of molecular basis
        basdims = [0, 16, 58, 78]

        # check individual overlap submatrices
        for i in range(3):
            for j in range(3):
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
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_2jk_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "2jk", 0.0, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat -= np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_2jkx_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "2jkx", 0.38, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat -= 0.38 * np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_j_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "j", 0.0, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_k_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "k", 0.0, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_kx_svpd_with_screener(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svpd, mol_h2o_dimer, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(t4c_drv, den_mat, "kx", 0.23, 0.0, 15)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat = 0.23 * np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_2jk_svpd(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svpd, mol_h2o_dimer, den_mat, "2jk",
                                    0.0, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat -= np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_2jkx_svpd(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svpd, mol_h2o_dimer, den_mat, "2jkx",
                                    0.38, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat -= 0.38 * np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_j_svpd(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svpd, mol_h2o_dimer, den_mat, "j", 0.0,
                                    0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.j.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_k_svpd(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svpd, mol_h2o_dimer, den_mat, "k", 0.0,
                                    0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

    def test_h2o_dimer_fock_kx_svpd(self):

        mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.density.npy')
        den_mat = make_matrix(bas_svpd, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))

        # compute Fock matrix
        fock_drv = FockDriver()
        fock_mat = fock_drv.compute(bas_svpd, mol_h2o_dimer, den_mat, "kx",
                                    0.22, 0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.svpd.k.npy')
        ref_mat = 0.22 * np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 16, 58, 78]

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
        fref = SubMatrix([0, 0, 78, 78])
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
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "2jk", 0.0,
                                    0.0)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.j.npy')
        ref_mat = 2.0 * np.load(npyfile)

        # load reference Fock matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.k.npy')
        ref_mat -= np.load(npyfile)

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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
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
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "2jkx", 0.68,
                                    0.0)

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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
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
        fock_mat = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, "kx", 0.68,
                                    0.0)

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
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full Fock matrix
        fmat = fock_mat.full_matrix()
        fref = SubMatrix([0, 0, 7, 7])
        fref.set_values(np.ascontiguousarray(ref_mat))

        assert fmat == fref

import pickle

from mpi4py import MPI
from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import AtomBasis
from veloxchem.mpitools import is_master
from tester import Tester


class TestAtomBasis:

    def get_hydrogen_svp_1s(self):

        bf = BasisFunction(
            [1.301070100000e+01, 1.962257200000e+00, 4.445379600000e-01],
            [1.968215800000e-02, 1.379652400000e-01, 4.783193500000e-01], 0)
        bf.normalize()

        return bf

    def get_hydrogen_svp_2s(self):

        bf = BasisFunction([1.219496200000e-01], [1.000000000000e+00], 0)
        bf.normalize()

        return bf

    def get_hydrogen_svp_1p(self):

        bf = BasisFunction([8.000000000000e-01], [1.000000000000e+00], 1)
        bf.normalize()

        return bf

    def get_hydrogen_svp(self):

        return AtomBasis([
            self.get_hydrogen_svp_1s(),
            self.get_hydrogen_svp_2s(),
            self.get_hydrogen_svp_1p()
        ], 1, 'DEF2-SVP', '')

    def get_lithium_svp_1s(self):

        bf = BasisFunction([
            2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00,
            2.450300905100e+00, 7.220957185500e-01
        ], [
            6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01,
            4.860657481700e-01, 4.362697795500e-01
        ], 0)
        bf.normalize()

        return bf

    def get_lithium_svp_2s(self):

        bf = BasisFunction([5.281088472100e-02], [1.000000000000e+00], 0)
        bf.normalize()

        return bf

    def get_lithium_svp_3s(self):

        bf = BasisFunction([2.096094879800e-02], [1.000000000000e+00], 0)
        bf.normalize()

        return bf

    def get_lithium_svp_1p(self):

        bf = BasisFunction([1.450000000000e+00, 3.000000000000e-01],
                           [2.586000000000e-01, 1.000000000000e+00], 1)
        bf.normalize()

        return bf

    def get_lithium_svp_2p(self):

        bf = BasisFunction([8.200000000000e-02], [1.000000000000e+00], 1)
        bf.normalize()

        return bf

    def get_lithium_svp(self):

        return AtomBasis([
            self.get_lithium_svp_1s(),
            self.get_lithium_svp_2s(),
            self.get_lithium_svp_3s(),
            self.get_lithium_svp_1p(),
            self.get_lithium_svp_2p()
        ], 3, 'DEF2-SVP', '')

    def test_pickle(self):

        bas_a = self.get_lithium_svp()

        # test pickling
        bobj = pickle.dumps(bas_a)
        bas_b = pickle.loads(bobj)
        Tester.compare_atom_basis(bas_a, bas_b)

    def test_set_and_get_identifier(self):

        bas = self.get_hydrogen_svp()

        assert bas.get_identifier() == 1

        bas.set_identifier(5)
        assert bas.get_identifier() == 5

    def test_set_and_get_name(self):

        bas = self.get_hydrogen_svp()

        assert bas.get_name() == 'DEF2-SVP'

        bas.set_name('DEF2-SVP(NEW)')
        assert bas.get_name() == 'DEF2-SVP(NEW)'

    def test_set_and_get_ecp_label(self):

        bas = self.get_hydrogen_svp()

        assert bas.get_ecp_label() == ''

        bas.set_ecp_label('DEF2-ECP')
        assert bas.get_ecp_label() == 'DEF2-ECP'

    def test_need_ecp(self):

        bas = self.get_hydrogen_svp()

        assert bas.need_ecp() is False

        bas.set_ecp_label('DEF2-ECP')
        assert bas.need_ecp() is True

    def test_add(self):

        bas_a = self.get_hydrogen_svp()

        bas_b = AtomBasis([], 1, 'DEF2-SVP', '')
        bas_b.add(self.get_hydrogen_svp_1s())
        bas_b.add(self.get_hydrogen_svp_2s())
        bas_b.add(self.get_hydrogen_svp_1p())

        Tester.compare_atom_basis(bas_a, bas_b)

    def test_reduce_to_valence_basis(self):

        bas_f = self.get_hydrogen_svp()
        bas_v = bas_f.reduce_to_valence_basis()

        bas_b = AtomBasis(
            [self.get_hydrogen_svp_1s(),
             self.get_hydrogen_svp_2s()], 1, 'DEF2-SVP(Valence)', '')

        Tester.compare_atom_basis(bas_b, bas_v)

    def test_get_basis_functions(self):

        bas = self.get_hydrogen_svp()

        bf1s = self.get_hydrogen_svp_1s()
        bf2s = self.get_hydrogen_svp_2s()
        bf1p = self.get_hydrogen_svp_1p()

        bfs_all = bas.get_basis_functions()
        assert len(bfs_all) == 3
        Tester.compare_basis_functions(bfs_all[0], bf1s)
        Tester.compare_basis_functions(bfs_all[1], bf2s)
        Tester.compare_basis_functions(bfs_all[2], bf1p)

        bfs_s = bas.get_basis_functions(0)
        assert len(bfs_s) == 2
        Tester.compare_basis_functions(bfs_s[0], bf1s)
        Tester.compare_basis_functions(bfs_s[1], bf2s)

        bfs_p = bas.get_basis_functions(1)
        assert len(bfs_p) == 1
        Tester.compare_basis_functions(bfs_p[0], bf1p)

        bfs_s_3 = bas.get_basis_functions(0, 3)
        assert len(bfs_s_3) == 1
        Tester.compare_basis_functions(bfs_s_3[0], bf1s)

        bfs_s_2 = bas.get_basis_functions(0, 2)
        assert len(bfs_s_2) == 0

        bfs_s_1 = bas.get_basis_functions(0, 1)
        assert len(bfs_s_1) == 1
        Tester.compare_basis_functions(bfs_s_1[0], bf2s)

    def test_max_angular_momentum(self):

        bas = self.get_hydrogen_svp()
        assert bas.max_angular_momentum() == 1

        bf1f = BasisFunction([1.219496200000e-01], [1.000000000000e+00], 4)
        bf1f.normalize()
        bas.add(bf1f)

        assert bas.max_angular_momentum() == 4

    def test_number_of_basis_functions(self):

        bas = self.get_lithium_svp()
        assert bas.number_of_basis_functions(0) == 3
        assert bas.number_of_basis_functions(1) == 2
        assert bas.number_of_basis_functions(2) == 0

        assert bas.number_of_basis_functions(0, 5) == 1
        assert bas.number_of_basis_functions(0, 3) == 0
        assert bas.number_of_basis_functions(0, 1) == 2

        assert bas.number_of_basis_functions(1, 3) == 0
        assert bas.number_of_basis_functions(1, 2) == 1
        assert bas.number_of_basis_functions(1, 1) == 1

        assert bas.number_of_basis_functions(2, 5) == 0
        assert bas.number_of_basis_functions(2, 3) == 0
        assert bas.number_of_basis_functions(2, 1) == 0

    def test_number_of_primitive_basis_functions(self):

        bas = self.get_lithium_svp()
        assert bas.number_of_primitive_basis_functions(0) == 7
        assert bas.number_of_primitive_basis_functions(1) == 3
        assert bas.number_of_primitive_basis_functions(2) == 0

    def test_contraction_depths(self):

        bas = self.get_lithium_svp()
        assert bas.contraction_depths(0) == {1, 5}
        assert bas.contraction_depths(1) == {1, 2}
        assert bas.contraction_depths(2) == set()

    def test_contraction_str(self):

        bas = self.get_lithium_svp()
        assert bas.contraction_str() == '(3S,2P)'

    def test_primitives_str(self):

        bas = self.get_lithium_svp()
        assert bas.primitives_str() == '(7S,3P)'

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if is_master(rank):
            bas_a = self.get_lithium_svp()
        else:
            bas_a = None
        bas_a = comm.bcast(bas_a)
        bas_b = self.get_lithium_svp()
        Tester.compare_atom_basis(bas_a, bas_b)

import pickle

from mpi4py import MPI
from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import AtomBasis
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.mpitools import is_master
from tester import Tester


class TestMolecularBasis:

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

    def get_hydrogen_svp_2p(self):

        bf = BasisFunction([1.170409905000e-01], [1.000000000000e+00], 1)
        bf.normalize()

        return bf

    def get_hydrogen_svp(self):

        return AtomBasis([
            self.get_hydrogen_svp_1s(),
            self.get_hydrogen_svp_2s(),
            self.get_hydrogen_svp_1p()
        ], 1, 'DEF2-SVP', '')

    def get_oxygen_svp_1s(self):

        bf = BasisFunction([
            2.266176778500e+03, 3.408701019100e+02, 7.736313516700e+01,
            2.147964494000e+01, 6.658943312400e+00
        ], [
            -5.343180992600e-03, -3.989003923000e-02, -1.785391198500e-01,
            -4.642768495900e-01, -4.430974517200e-01
        ], 0)
        bf.normalize()

        return bf

    def get_oxygen_svp_2s(self):

        bf = BasisFunction([8.097597566800e-01], [1.000000000000e+00], 0)
        bf.normalize()

        return bf

    def get_oxygen_svp_3s(self):

        bf = BasisFunction([2.553077223400e-01], [1.000000000000e+00], 0)
        bf.normalize()

        return bf

    def get_oxygen_svp_1p(self):

        bf = BasisFunction(
            [1.772150431700e+01, 3.863550544000e+00, 1.048092088300e+00],
            [4.339457319300e-02, 2.309412076500e-01, 5.137531106400e-01], 1)
        bf.normalize()

        return bf

    def get_oxygen_svp_2p(self):

        bf = BasisFunction([2.764154441100e-01], [1.000000000000e+00], 1)
        bf.normalize()

        return bf

    def get_oxygen_svp_1d(self):

        bf = BasisFunction([1.200000000000e+00], [1.000000000000e+00], 2)
        bf.normalize()

        return bf

    def get_oxygen_svp(self):

        return AtomBasis([
            self.get_oxygen_svp_1s(),
            self.get_oxygen_svp_2s(),
            self.get_oxygen_svp_3s(),
            self.get_oxygen_svp_1p(),
            self.get_oxygen_svp_2p(),
            self.get_oxygen_svp_1d()
        ], 8, 'DEF2-SVP', '')

    def get_hydrogen_svp_red(self):

        return AtomBasis(
            [self.get_hydrogen_svp_1s(),
             self.get_hydrogen_svp_2s()], 1, 'DEF2-SVP(Valence)', '')

    def get_oxygen_svp_red(self):

        return AtomBasis([
            self.get_oxygen_svp_1s(),
            self.get_oxygen_svp_2s(),
            self.get_oxygen_svp_3s(),
            self.get_oxygen_svp_1p(),
            self.get_oxygen_svp_2p()
        ], 8, 'DEF2-SVP(Valence)', '')

    def get_hydrogen_svpd(self):

        return AtomBasis([
            self.get_hydrogen_svp_1s(),
            self.get_hydrogen_svp_2s(),
            self.get_hydrogen_svp_1p(),
            self.get_hydrogen_svp_2p()
        ], 1, 'DEF2-SVPD', '')

    def get_h2o(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        return Molecule.read_str(h2ostr, 'au')

    def get_h2o_svp(self):

        o_bas = self.get_oxygen_svp()
        h_bas = self.get_hydrogen_svp()

        return MolecularBasis([o_bas, h_bas], [0, 1, 1])

    def get_h2o_mixed(self):

        o_bas = self.get_oxygen_svp()
        h_bas = self.get_hydrogen_svp()
        h_basp = self.get_hydrogen_svpd()

        return MolecularBasis([o_bas, h_bas, h_basp], [0, 1, 2])

    def test_constructor(self):

        h2o = self.get_h2o()

        a_basis = MolecularBasis.read(h2o, 'DEF2-SVP', 'basis')
        b_basis = self.get_h2o_svp()
        Tester.compare_molecular_basis(a_basis, b_basis)

    def test_constructor_dict(self):

        h2o = self.get_h2o()

        bdict = {'3': 'DEF2-SVPD'}

        a_basis = MolecularBasis.read_dict(h2o, 'DEF2-SVP', bdict, 'basis')
        b_basis = self.get_h2o_mixed()
        Tester.compare_molecular_basis(a_basis, b_basis)

    def test_pickle(self):

        a_basis = self.get_h2o_svp()

        # test pickling
        bobj = pickle.dumps(a_basis)
        b_basis = pickle.loads(bobj)
        Tester.compare_molecular_basis(a_basis, b_basis)

    def test_add(self):

        a_basis = self.get_h2o_svp()

        # test add
        o_bas = self.get_oxygen_svp()
        h_bas = self.get_hydrogen_svp()

        b_basis = MolecularBasis()
        b_basis.add(o_bas)
        b_basis.add(h_bas)
        b_basis.add(h_bas)

        Tester.compare_molecular_basis(a_basis, b_basis)

    def test_reduce_to_valence_basis(self):

        basis = self.get_h2o_svp()
        a_basis = basis.reduce_to_valence_basis()

        o_bas = self.get_oxygen_svp_red()
        h_bas = self.get_hydrogen_svp_red()

        b_basis = MolecularBasis()
        b_basis.add(o_bas)
        b_basis.add(h_bas)
        b_basis.add(h_bas)

        Tester.compare_molecular_basis(a_basis, b_basis)

    def test_get_basis_sets(self):

        basis = self.get_h2o_svp()
        bsets = basis.get_basis_sets()

        assert len(bsets) == 2

        o_bas = self.get_oxygen_svp()
        Tester.compare_atom_basis(o_bas, bsets[0])

        h_bas = self.get_hydrogen_svp()
        Tester.compare_atom_basis(h_bas, bsets[1])

    def test_get_basis_sets_indexes(self):

        basis = self.get_h2o_svp()
        indexes = basis.get_basis_sets_indexes()

        assert indexes == [0, 1, 1]

    def test_max_angular_momentum(self):

        basis = self.get_h2o_svp()

        assert basis.max_angular_momentum() == 2

        assert basis.max_angular_momentum([0, 1]) == 2

        assert basis.max_angular_momentum([0, 2]) == 2

        assert basis.max_angular_momentum([1, 2]) == 1

    def test_get_basis_functions(self):

        # oxygen basis functions
        o1s = self.get_oxygen_svp_1s()
        o2s = self.get_oxygen_svp_2s()
        o3s = self.get_oxygen_svp_3s()
        o1p = self.get_oxygen_svp_1p()
        o2p = self.get_oxygen_svp_2p()
        o1d = self.get_oxygen_svp_1d()

        # hydrogen basis functions
        h1s = self.get_hydrogen_svp_1s()
        h2s = self.get_hydrogen_svp_2s()
        h1p = self.get_hydrogen_svp_1p()

        # set up basis set for water
        basis = self.get_h2o_svp()

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions()
        b_bfs = [o1s, o2s, o3s, o1p, o2p, o1d, h1s, h2s, h1p, h1s, h2s, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2])
        b_bfs = [h1s, h2s, h1p, h1s, h2s, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2])
        b_bfs = [o1s, o2s, o3s, o1p, o2p, o1d, h1s, h2s, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(0)
        b_bfs = [o1s, o2s, o3s, h1s, h2s, h1s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 0)
        b_bfs = [h1s, h2s, h1s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 0)
        b_bfs = [o1s, o2s, o3s, h1s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(1)
        b_bfs = [o1p, o2p, h1p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 1)
        b_bfs = [h1p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 1)
        b_bfs = [o1p, o2p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(2)
        b_bfs = [o1d]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 2)
        assert len(a_bfs) == 0

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 2)
        b_bfs = [o1d]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(0, 5)
        b_bfs = [o1s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 0, 5)
        assert len(a_bfs) == 0
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 0, 5)
        b_bfs = [o1s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(0, 3)
        b_bfs = [h1s, h1s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 0, 3)
        b_bfs = [h1s, h1s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 0, 3)
        b_bfs = [h1s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(0, 1)
        b_bfs = [o2s, o3s, h2s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 0, 1)
        b_bfs = [h2s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 0, 1)
        b_bfs = [o2s, o3s, h2s]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(1, 3)
        b_bfs = [o1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 1, 3)
        assert len(a_bfs) == 0
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 1, 3)
        b_bfs = [o1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions(1, 1)
        b_bfs = [o2p, h1p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([1, 2], 1, 1)
        b_bfs = [h1p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

        # test generic getter for basis functions
        a_bfs = basis.get_basis_functions([0, 2], 1, 1)
        b_bfs = [o2p, h1p]
        assert len(a_bfs) == len(b_bfs)
        for bfa, bfb in zip(a_bfs, b_bfs):
            Tester.compare_basis_functions(bfa, bfb)

    def test_get_atomic_indexes(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes()
        b_indexes = [0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2])
        b_indexes = [1, 1, 1, 2, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2])
        b_indexes = [0, 0, 0, 0, 0, 0, 2, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(0)
        b_indexes = [0, 0, 0, 1, 1, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 0)
        b_indexes = [1, 1, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 0)
        b_indexes = [0, 0, 0, 2, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(1)
        b_indexes = [0, 0, 1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 1)
        b_indexes = [1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 1)
        b_indexes = [0, 0, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(2)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 2)
        assert len(a_indexes) == 0

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 2)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(0, 5)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 0, 5)
        assert len(a_indexes) == 0

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 0, 5)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(0, 3)
        b_indexes = [1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 0, 3)
        b_indexes = [1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 0, 3)
        b_indexes = [2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(0, 1)
        b_indexes = [0, 0, 1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 0, 1)
        b_indexes = [1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 0, 1)
        b_indexes = [0, 0, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(1, 3)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 1, 3)
        assert len(a_indexes) == 0

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 1, 3)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(1, 1)
        b_indexes = [0, 1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 1, 1)
        b_indexes = [1, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 1, 1)
        b_indexes = [0, 2]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes(2, 1)
        b_indexes = [0]
        assert a_indexes == b_indexes

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([1, 2], 2, 1)
        assert len(a_indexes) == 0

        # test getter for atomic indexes
        a_indexes = basis.get_atomic_indexes([0, 2], 2, 1)
        b_indexes = [0]
        assert a_indexes == b_indexes

    def test_number_of_basis_functions(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test number of basis functions
        assert basis.number_of_basis_functions(0) == 7
        assert basis.number_of_basis_functions(1) == 4
        assert basis.number_of_basis_functions(2) == 1
        assert basis.number_of_basis_functions(3) == 0

        # test number of basis functions
        assert basis.number_of_basis_functions([1, 2], 0) == 4
        assert basis.number_of_basis_functions([1, 2], 1) == 2
        assert basis.number_of_basis_functions([1, 2], 2) == 0
        assert basis.number_of_basis_functions([1, 2], 3) == 0

        # test number of basis functions
        assert basis.number_of_basis_functions([0, 2], 0) == 5
        assert basis.number_of_basis_functions([0, 2], 1) == 3
        assert basis.number_of_basis_functions([0, 2], 2) == 1
        assert basis.number_of_basis_functions([0, 2], 3) == 0

        # test number of basis functions
        assert basis.number_of_basis_functions(0, 5) == 1
        assert basis.number_of_basis_functions(0, 3) == 2
        assert basis.number_of_basis_functions(0, 1) == 4
        assert basis.number_of_basis_functions(1, 3) == 1
        assert basis.number_of_basis_functions(1, 1) == 3
        assert basis.number_of_basis_functions(2, 1) == 1

        # test number of basis functions
        assert basis.number_of_basis_functions([1, 2], 0, 5) == 0
        assert basis.number_of_basis_functions([1, 2], 0, 3) == 2
        assert basis.number_of_basis_functions([1, 2], 0, 1) == 2
        assert basis.number_of_basis_functions([1, 2], 1, 3) == 0
        assert basis.number_of_basis_functions([1, 2], 1, 1) == 2
        assert basis.number_of_basis_functions([1, 2], 2, 1) == 0

        # test number of basis functions
        assert basis.number_of_basis_functions([0, 2], 0, 5) == 1
        assert basis.number_of_basis_functions([0, 2], 0, 3) == 1
        assert basis.number_of_basis_functions([0, 2], 0, 1) == 3
        assert basis.number_of_basis_functions([0, 2], 1, 3) == 1
        assert basis.number_of_basis_functions([0, 2], 1, 1) == 2
        assert basis.number_of_basis_functions([0, 2], 2, 1) == 1

    def test_number_of_primitive_basis_functions(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test number of primitive basis functions
        assert basis.number_of_primitive_basis_functions(0) == 15
        assert basis.number_of_primitive_basis_functions(1) == 6
        assert basis.number_of_primitive_basis_functions(2) == 1
        assert basis.number_of_primitive_basis_functions(3) == 0

        # test number of primitive basis functions
        assert basis.number_of_primitive_basis_functions([1, 2], 0) == 8
        assert basis.number_of_primitive_basis_functions([1, 2], 1) == 2
        assert basis.number_of_primitive_basis_functions([1, 2], 2) == 0

        # test number of primitive basis functions
        assert basis.number_of_primitive_basis_functions([0, 2], 0) == 11
        assert basis.number_of_primitive_basis_functions([0, 2], 1) == 5
        assert basis.number_of_primitive_basis_functions([0, 2], 2) == 1

    def test_contraction_depths(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test contraction depths
        assert basis.contraction_depths(0) == {1, 3, 5}
        assert basis.contraction_depths(1) == {1, 3}
        assert basis.contraction_depths(2) == {1}

        # test contraction depths
        assert basis.contraction_depths([1, 2], 0) == {1, 3}
        assert basis.contraction_depths([1, 2], 1) == {1}
        assert basis.contraction_depths([1, 2], 2) == set()

        # test contraction depths
        assert basis.contraction_depths([0, 2], 0) == {1, 3, 5}
        assert basis.contraction_depths([0, 2], 1) == {1, 3}
        assert basis.contraction_depths([0, 2], 2) == {1}

    def test_get_dimensions_of_basis(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test basis dimensions
        assert basis.get_dimensions_of_basis() == 24

        # test partial basis dimensions
        assert basis.get_dimensions_of_basis(0) == 0
        assert basis.get_dimensions_of_basis(1) == 7
        assert basis.get_dimensions_of_basis(2) == 19
        assert basis.get_dimensions_of_basis(3) == 24

    def test_get_dimensions_of_primitive_basis(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test basis dimensions
        assert basis.get_dimensions_of_primitive_basis() == 38

    def test_get_index_map(self):

        # set up basis for water
        basis = self.get_h2o_svp()

        # test compressed indexes map
        assert basis.get_index_map(0, 5) == [7, 0]
        assert basis.get_index_map(0, 3) == [7, 3, 5]
        assert basis.get_index_map(0, 1) == [7, 1, 2, 4, 6]

        # test compressed indexes map
        assert basis.get_index_map(1, 5) == [4]
        assert basis.get_index_map(1, 3) == [4, 7]
        assert basis.get_index_map(1, 1) == [4, 8, 9, 10]

        # test compressed indexes map
        assert basis.get_index_map(2, 5) == [1]
        assert basis.get_index_map(2, 3) == [1]
        assert basis.get_index_map(2, 1) == [1, 19]

        # test compressed indexes map
        assert basis.get_index_map([1, 2], 0, 5) == [7]
        assert basis.get_index_map([1, 2], 0, 3) == [7, 3, 5]
        assert basis.get_index_map([1, 2], 0, 1) == [7, 4, 6]

        # test compressed indexes map
        assert basis.get_index_map([1, 2], 1, 5) == [4]
        assert basis.get_index_map([1, 2], 1, 3) == [4]
        assert basis.get_index_map([1, 2], 1, 1) == [4, 9, 10]

        # test compressed indexes map
        assert basis.get_index_map([1, 2], 2, 5) == [1]
        assert basis.get_index_map([1, 2], 2, 3) == [1]
        assert basis.get_index_map([1, 2], 2, 1) == [1]

        # test compressed indexes map
        assert basis.get_index_map([0, 2], 0, 5) == [7, 0]
        assert basis.get_index_map([0, 2], 0, 3) == [7, 5]
        assert basis.get_index_map([0, 2], 0, 1) == [7, 1, 2, 6]

        # test compressed indexes map
        assert basis.get_index_map([0, 2], 1, 5) == [4]
        assert basis.get_index_map([0, 2], 1, 3) == [4, 7]
        assert basis.get_index_map([0, 2], 1, 1) == [4, 8, 10]

        # test compressed indexes map
        assert basis.get_index_map([0, 2], 2, 5) == [1]
        assert basis.get_index_map([0, 2], 2, 3) == [1]
        assert basis.get_index_map([0, 2], 2, 1) == [1, 19]

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if is_master(rank):
            a_basis = self.get_h2o_svp()
        else:
            a_basis = None
        a_basis = comm.bcast(a_basis)
        b_basis = self.get_h2o_svp()
        Tester.compare_molecular_basis(a_basis, b_basis)

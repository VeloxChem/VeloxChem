import numpy as np
import math as mt


class Tester:

    @staticmethod
    def compare_molecules(lhs, rhs):

        tol = 1.0e-12

        # check multiplicities
        m_a = lhs.get_multiplicity()
        m_b = rhs.get_multiplicity()
        assert m_a == m_b

        # check charges
        q_a = lhs.get_charge()
        q_b = rhs.get_charge()
        assert mt.isclose(q_a, q_b, rel_tol=tol, abs_tol=tol)

        # check elemental identifiers
        elem_a = lhs.get_identifiers()
        elem_b = rhs.get_identifiers()
        assert elem_a == elem_b

        # check coordinates
        coords_a = np.array(lhs.get_coordinates())
        coords_b = np.array(rhs.get_coordinates())
        assert np.allclose(coords_a, coords_b, tol, tol, False)

    @staticmethod
    def compare_basis_functions(lhs, rhs):

        tol = 1.0e-12

        # check angular momentum
        l_a = lhs.get_angular_momentum()
        l_b = rhs.get_angular_momentum()
        assert l_a == l_b

        # check exponents
        fe_a = np.array(lhs.get_exponents())
        fe_b = np.array(rhs.get_exponents())
        assert np.allclose(fe_a, fe_b, tol, tol, False)

        # check normalization factors
        fn_a = np.array(lhs.get_normalization_factors())
        fn_b = np.array(rhs.get_normalization_factors())
        assert np.allclose(fn_a, fn_b, tol, tol, False)

    @staticmethod
    def compare_atom_basis(lhs, rhs):

        # check elemental identifier
        id_a = lhs.get_identifier()
        id_b = rhs.get_identifier()
        assert id_a == id_b

        # check name
        label_a = lhs.get_name()
        label_b = rhs.get_name()
        assert label_a == label_b

        # check ecp labels
        label_a = lhs.get_ecp_label()
        label_b = rhs.get_ecp_label()
        assert label_a == label_b

        # check basis functions
        bfs_a = lhs.get_basis_functions()
        bfs_b = rhs.get_basis_functions()
        assert len(bfs_a) == len(bfs_b)
        for bf_a, bf_b in zip(bfs_a, bfs_b):
            Tester.compare_basis_functions(bf_a, bf_b)

    @staticmethod
    def compare_molecular_basis(lhs, rhs):

        # check atom basis sets
        bas_a = lhs.get_basis_sets()
        bas_b = rhs.get_basis_sets()
        assert len(bas_a) == len(bas_b)
        for ba, bb in zip(bas_a, bas_b):
            Tester.compare_atom_basis(ba, bb)

        # check basis sets indexes
        bidx_a = lhs.get_basis_sets_indexes()
        bidx_b = rhs.get_basis_sets_indexes()
        assert bidx_a == bidx_b

    @staticmethod
    def compare_gto_blocks(lhs, rhs):

        tol = 1.0e-12

        # check coordinates
        coords_a = np.array(lhs.get_coordinates())
        coords_b = np.array(rhs.get_coordinates())
        assert np.allclose(coords_a, coords_b, tol, tol, False)

        # check exponents
        fexps_a = np.array(lhs.get_exponents())
        fexps_b = np.array(rhs.get_exponents())
        assert np.allclose(fexps_a, fexps_b, tol, tol, False)

        # check normalization factors
        fnorms_a = np.array(lhs.get_normalization_factors())
        fnorms_b = np.array(rhs.get_normalization_factors())
        assert np.allclose(fnorms_a, fnorms_b, tol, tol, False)

        # check orbital indexes
        orb_idx_a = lhs.get_orbital_indexes()
        orb_idx_b = rhs.get_orbital_indexes()
        assert orb_idx_a == orb_idx_b

        # check atomic indexes
        atm_idx_a = lhs.get_atomic_indexes()
        atm_idx_b = rhs.get_atomic_indexes()
        assert atm_idx_a == atm_idx_b

        # check angular momentum
        angmom_a = lhs.get_angular_momentum()
        angmom_b = rhs.get_angular_momentum()
        assert angmom_a == angmom_b

        # check number of primitives
        nprims_a = lhs.number_of_primitives()
        nprims_b = rhs.number_of_primitives()
        assert nprims_a == nprims_b

    @staticmethod
    def compare_submatrices(lhs, rhs):

        tol = 1.0e-12

        # check dimensions
        dims_a = lhs.get_dimensions()
        dims_b = rhs.get_dimensions()
        assert dims_a == dims_b

        # check submatrix values
        mat_a = lhs.to_numpy()
        mat_b = rhs.to_numpy()
        assert np.allclose(mat_a, mat_b, tol, tol, False)

    @staticmethod
    def compare_matrices(lhs, rhs):

        # check list of submatrices
        keys_a = lhs.get_angular_pairs()
        keys_b = rhs.get_angular_pairs()
        assert len(keys_a) == len(keys_b)
        for key_a, key_b in zip(keys_a, keys_b):
            assert key_a == key_b
            ma = lhs.get_submatrix(key_a)
            mb = rhs.get_submatrix(key_b)
            Tester.compare_submatrices(ma, mb)

        # check matrix type
        assert lhs.get_type() == rhs.get_type()

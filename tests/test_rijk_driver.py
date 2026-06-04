from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import T4CScreener
from veloxchem.veloxchemlib import make_matrix, mat_t, SubMatrix
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.rijkfockdriver import RIJKFockDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.fockdriver import FockDriver


class TestRIMODriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jkfit')

        return mol, bas, aux_bas

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='single MPI process only')
    def test_h2o_compute_fock(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol_h2o, bas_sto3g)

        dmat = make_matrix(bas_sto3g, mat_t.symmetric)
        dmat.set_values(scf_results['D_alpha'])

        molorbs = scf_drv.molecular_orbitals

        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_aux, verbose=False)
        ri_fock_drv.compute_bq_vectors(mol_h2o,
                                       bas_sto3g,
                                       bas_aux,
                                       verbose=False)

        jmat = ri_fock_drv.compute_j_fock(dmat, 'j', verbose=False)
        kmat = ri_fock_drv.compute_k_fock(dmat, molorbs, verbose=False)

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        jmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "j", 0.0, 0.0, 15)
        kmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "k", 0.0, 0.0, 15)

        assert np.max(np.abs(kmat_ref.to_numpy() - kmat.to_numpy())) < 1e-3
        assert np.max(np.abs(jmat_ref.to_numpy() - jmat.to_numpy())) < 1e-3

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='single MPI process only')
    def test_h2o_compute_mo_bq_vectors(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol_h2o, bas_sto3g)

        # Compute reference MO transformed Fock matrices (vv block)
        dmat = make_matrix(bas_sto3g, mat_t.symmetric)
        dmat.set_values(scf_results['D_alpha'])

        molorbs = scf_drv.molecular_orbitals

        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_aux, verbose=False)
        ri_fock_drv.compute_bq_vectors(mol_h2o,
                                       bas_sto3g,
                                       bas_aux,
                                       verbose=False)

        jmat = ri_fock_drv.compute_j_fock(dmat, 'j', verbose=False)
        kmat = ri_fock_drv.compute_k_fock(dmat, molorbs, verbose=False)

        screener = T4CScreener()
        screener.partition(bas_sto3g, mol_h2o, 'eri')
        ri_fock_drv.compute_screened_bq_vectors(screener,
                                                mol_h2o,
                                                bas_aux,
                                                12,
                                                verbose=False)

        no = mol_h2o.number_of_alpha_occupied_orbitals(bas_sto3g)
        nv = bas_sto3g.get_dimensions_of_basis() - no
        co = scf_results['C_alpha'][:, :no]
        cv = scf_results['C_alpha'][:, no:]

        jmat_ref = np.linalg.multi_dot([cv.T, jmat.to_numpy(), cv])
        kmat_ref = np.linalg.multi_dot([cv.T, kmat.to_numpy(), cv])

        # Compute MO transformed Fock matrix using the MO transformed Bq vectors
        n_aux = bas_aux.get_dimensions_of_basis()
        co_matrix = SubMatrix([0, 0, molorbs.number_aos(), no])
        co_matrix.set_values(co)
        cv_matrix = SubMatrix([0, 0, molorbs.number_aos(), nv])
        cv_matrix.set_values(cv)

        # Batch over n_aux
        batches = 3
        batch_dim = n_aux // batches
        aux_start = 0

        jmat = np.zeros_like(jmat_ref)
        kmat = np.zeros_like(jmat_ref)

        for batch in range(batches):

            if batch == batches - 1:
                batch_dim = n_aux - (batches - 1) * batch_dim

            bq_oo_buf = ri_fock_drv.compute_mo_bq_vectors(
                co_matrix, co_matrix, aux_start, aux_start + batch_dim)
            bq_ov_buf = ri_fock_drv.compute_mo_bq_vectors(
                co_matrix, cv_matrix, aux_start, aux_start + batch_dim)
            bq_vv_buf = ri_fock_drv.compute_mo_bq_vectors(
                cv_matrix, cv_matrix, aux_start, aux_start + batch_dim)

            bq_oo = np.zeros((batch_dim, no, no))
            bq_ov = np.zeros((batch_dim, no, nv))
            bq_vv = np.zeros((batch_dim, nv, nv))

            for i in range(batch_dim):
                bq_oo[i] = bq_oo_buf[i].to_numpy().reshape(no, no)
                bq_vv[i] = bq_vv_buf[i].to_numpy().reshape(nv, nv)
                bq_ov[i] = bq_ov_buf[i].to_numpy().reshape(no, nv)

            jmat += np.einsum('xii,xab->ab', bq_oo, bq_vv)
            kmat += np.einsum('xia,xib->ab', bq_ov, bq_ov)

            aux_start += batch_dim

        np.testing.assert_allclose(jmat, jmat_ref, atol=1e-12)
        np.testing.assert_allclose(kmat, kmat_ref, atol=1e-12)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='single MPI process only')
    def test_h2o_local_compute_mo_bq_vectors(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol_h2o, bas_sto3g)

        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_aux, verbose=False)
        screener = T4CScreener()
        screener.partition(bas_sto3g, mol_h2o, 'eri')
        ri_fock_drv.compute_screened_bq_vectors(screener,
                                                mol_h2o,
                                                bas_aux,
                                                12,
                                                verbose=False)

        no = mol_h2o.number_of_alpha_occupied_orbitals(bas_sto3g)
        nv = bas_sto3g.get_dimensions_of_basis() - no
        co = scf_results['C_alpha'][:, :no]
        cv = scf_results['C_alpha'][:, no:]

        n_aux = bas_aux.get_dimensions_of_basis()
        co_matrix = SubMatrix(
            [0, 0, scf_drv.molecular_orbitals.number_aos(), no])
        co_matrix.set_values(co)
        cv_matrix = SubMatrix(
            [0, 0, scf_drv.molecular_orbitals.number_aos(), nv])
        cv_matrix.set_values(cv)

        # Indices to compute
        indices = [4, 7, 12, 52, 104]

        # Compute reference MO bq vectors
        bq_oo_buf = ri_fock_drv.compute_mo_bq_vectors(co_matrix, co_matrix, 0,
                                                      n_aux)
        bq_ov_buf = ri_fock_drv.compute_mo_bq_vectors(co_matrix, cv_matrix, 0,
                                                      n_aux)
        bq_vv_buf = ri_fock_drv.compute_mo_bq_vectors(cv_matrix, cv_matrix, 0,
                                                      n_aux)

        bq_oo_ref = np.zeros((len(indices), no, no))
        bq_ov_ref = np.zeros((len(indices), no, nv))
        bq_vv_ref = np.zeros((len(indices), nv, nv))

        for i, index in enumerate(indices):
            bq_oo_ref[i] = bq_oo_buf[index].to_numpy().reshape(no, no)
            bq_vv_ref[i] = bq_vv_buf[index].to_numpy().reshape(nv, nv)
            bq_ov_ref[i] = bq_ov_buf[index].to_numpy().reshape(no, nv)

        # Compute local MO bq vectors, i.e. only the requested indices
        ri_fock_drv.local_compute_screened_bq_vectors(screener,
                                                      mol_h2o,
                                                      bas_aux,
                                                      indices,
                                                      12,
                                                      verbose=True)

        bq_oo_buf = ri_fock_drv.compute_mo_bq_vectors(co_matrix, co_matrix, 0,
                                                      len(indices))
        bq_ov_buf = ri_fock_drv.compute_mo_bq_vectors(co_matrix, cv_matrix, 0,
                                                      len(indices))
        bq_vv_buf = ri_fock_drv.compute_mo_bq_vectors(cv_matrix, cv_matrix, 0,
                                                      len(indices))

        bq_oo = np.zeros((len(indices), no, no))
        bq_ov = np.zeros((len(indices), no, nv))
        bq_vv = np.zeros((len(indices), nv, nv))

        for i in range(len(indices)):
            bq_oo[i] = bq_oo_buf[i].to_numpy().reshape(no, no)
            bq_vv[i] = bq_vv_buf[i].to_numpy().reshape(nv, nv)
            bq_ov[i] = bq_ov_buf[i].to_numpy().reshape(no, nv)

        np.testing.assert_allclose(bq_oo, bq_oo_ref, atol=1e-12)
        np.testing.assert_allclose(bq_vv, bq_vv_ref, atol=1e-12)
        np.testing.assert_allclose(bq_ov, bq_ov_ref, atol=1e-12)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='single MPI process only')
    def test_h2o_bq_memory_estimate(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        ri_fock_drv = RIJKFockDriver()
        screener = T4CScreener()
        screener.partition(bas_sto3g, mol_h2o, 'eri')

        estimate = ri_fock_drv.estimate_memory_for_bq_vectors(
            screener, mol_h2o, bas_aux, 12)

        assert estimate == 25312

from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import T4CScreener
from veloxchem.veloxchemlib import PackedMatrix, make_matrix, mat_t
from veloxchem.veloxchemlib import SimdRIJKFockDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
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

    def get_simd_fock(self, mol, bas, aux_bas, density, coeffs):
        """
        Returns the Coulomb and the exchange the simd driver forms for the given
        density and coefficients, taken apart from the composed matrix it
        answers with.
        """

        nao = bas.get_dimensions_of_basis()

        packed_density = PackedMatrix(nao, nao, mat_t.symmetric)
        packed_density.from_numpy(np.ascontiguousarray(density))

        packed_coeffs = PackedMatrix(nao, coeffs.shape[1], mat_t.general)
        packed_coeffs.from_numpy(np.ascontiguousarray(coeffs))

        ri_fock_drv = SimdRIJKFockDriver()
        ri_fock_drv.prepare(mol, bas, aux_bas, 0.0, 1 << 30)

        # NOTE: the driver answers with twice the Coulomb less the scaled
        # exchange, so a factor of zero leaves twice the Coulomb and a factor of
        # one removes the exchange from it.

        two_jmat = ri_fock_drv.compute(packed_density, packed_coeffs,
                                       0.0).to_numpy(max_memory=8.0)
        fock_mat = ri_fock_drv.compute(packed_density, packed_coeffs,
                                       1.0).to_numpy(max_memory=8.0)

        return 0.5 * two_jmat, two_jmat - fock_mat

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
        nocc = int(np.sum(molorbs.occa_to_numpy()))

        jmat, kmat = self.get_simd_fock(mol_h2o, bas_sto3g, bas_aux,
                                        scf_results['D_alpha'],
                                        molorbs.alpha_to_numpy()[:, :nocc])

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        jmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "j", 0.0, 0.0, 15)
        kmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "k", 0.0, 0.0, 15)

        assert np.max(np.abs(kmat_ref.to_numpy() - kmat)) < 1e-3
        assert np.max(np.abs(jmat_ref.to_numpy() - jmat)) < 1e-3

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='single MPI process only')
    def test_density_factor_matches_exact_k_for_fractional_density(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(mol_h2o, bas_sto3g)

        C = scf_drv.molecular_orbitals.alpha_to_numpy()
        occ = np.zeros(C.shape[1])
        occ[:5] = 0.9
        occ[5] = 0.1
        d_np = C @ np.diag(occ) @ C.T
        d_np = 0.5 * (d_np + d_np.T)

        dmat = make_matrix(bas_sto3g, mat_t.symmetric)
        dmat.set_values(d_np)

        # NOTE: the exchange is formed from a factor of the density, which the
        # SCF driver builds when the density is not idempotent.

        density_factor = scf_drv._simd_ri_jk_occupied(
            d_np, 'alpha', use_density_factor=True)

        _, kmat = self.get_simd_fock(mol_h2o, bas_sto3g, bas_aux, d_np,
                                     density_factor)

        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")

        fock_drv = FockDriver()
        kmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "k", 0.0, 0.0, 15)

        assert np.max(np.abs(kmat_ref.to_numpy() - kmat)) < 1e-3

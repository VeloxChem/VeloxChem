import numpy as np

from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfUnrestrictedDriver(ScfDriver):
    """
    Implements spin unrestricted open shell SCF method with C2-DIIS and
    two-level C2-DIIS convergence accelerators.
    """

    def __init__(self, comm, ostream):
        """
        Initializes spin unrestricted open shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        super().__init__(comm, ostream)

        self.restricted = False

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin unrestricted open shell electronic gradient using
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param den_mat:
            The density matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The electronic gradient.
        """

        if self.rank == mpi_master():
            smat = ovl_mat.to_numpy()
            tmat = oao_mat.to_numpy()

            dmat_a = den_mat.alpha_to_numpy(0)
            dmat_b = den_mat.beta_to_numpy(0)

            fmat_a = fock_mat.alpha_to_numpy(0)
            fmat_b = fock_mat.beta_to_numpy(0)

            fds_a = np.matmul(fmat_a, np.matmul(dmat_a, smat))
            fds_b = np.matmul(fmat_b, np.matmul(dmat_b, smat))

            e_grad_a = np.linalg.norm(
                np.matmul(tmat.T, np.matmul(fds_a - fds_a.T, tmat)))

            e_grad_b = np.linalg.norm(
                np.matmul(tmat.T, np.matmul(fds_b - fds_b.T, tmat)))

            e_grad = e_grad_a + e_grad_b
        else:
            e_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())

        return e_grad

    def comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of spin unrestricted open shell density change between
        two density matrices. Overloaded base class method.

        :param den_mat:
            The current density matrix.
        :param old_den_mat:
            The previous density matrix.

        :return:
            The norm of change between two density matrices.
        """

        if self.rank == mpi_master():
            diff_mat = den_mat.sub(old_den_mat)
            ddmat_a = diff_mat.alpha_to_numpy(0)
            ddmat_b = diff_mat.beta_to_numpy(0)

            diff_den_a = np.linalg.norm(ddmat_a)
            diff_den_b = np.linalg.norm(ddmat_b)

            diff_den = max(diff_den_a, diff_den_b)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def store_diis_data(self, i, fock_mat, den_mat):
        """
        Stores spin unrestricted open shell Fock/Kohn-Sham and density matrices
        for current iteration. Overloaded base class method.

        :param i:
            The number of current SCF iteration.
        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param den_mat:
            The density matrix.
        """

        if self.rank == mpi_master():

            if not self.skip_iter:

                if len(self.fock_matrices) == self.max_err_vecs:

                    self.fock_matrices.popleft()
                    self.den_matrices.popleft()

                    self.fock_matrices_beta.popleft()
                    self.den_matrices_beta.popleft()

                self.fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self.den_matrices.append(den_mat.alpha_to_numpy(0))

                self.fock_matrices_beta.append(fock_mat.beta_to_numpy(0))
                self.den_matrices_beta.append(den_mat.beta_to_numpy(0))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        in OAO basis by applying Lowdin or canonical orthogonalization to AO
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The effective Fock/Kohn-Sham matrices.
        """

        if self.rank == mpi_master():

            if len(self.fock_matrices) == 1:

                return (np.copy(self.fock_matrices[0]),
                        np.copy(self.fock_matrices_beta[0]))

            if len(self.fock_matrices) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_unrestricted_error_vectors(
                    self.fock_matrices, self.fock_matrices_beta,
                    self.den_matrices, self.den_matrices_beta, ovl_mat, oao_mat)

                weights = acc_diis.compute_weights()

                return self.get_scaled_fock(weights)

            return fock_mat.alpha_to_numpy(0), fock_mat.beta_to_numpy(0)

        return None, None

    def get_scaled_fock(self, weights):
        """
        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrices.
        """

        effmat_a = np.zeros(self.fock_matrices[0].shape)
        effmat_b = np.zeros(self.fock_matrices_beta[0].shape)

        for w, fa, fb in zip(weights, self.fock_matrices,
                             self.fock_matrices_beta):

            effmat_a += w * fa
            effmat_b += w * fb

        return effmat_a, effmat_b

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """
        Generates spin unrestricted molecular orbital by diagonalizing
        spin unrestricted open shell Fock/Kohn-Sham matrix. Overloaded base
        class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The molecular orbitals.
        """

        if self.rank == mpi_master():

            tmat = oao_mat.to_numpy()

            fmo_a = np.matmul(tmat.T, np.matmul(fock_mat[0], tmat))
            fmo_b = np.matmul(tmat.T, np.matmul(fock_mat[1], tmat))

            eigs_a, evecs_a = np.linalg.eigh(fmo_a)
            eigs_b, evecs_b = np.linalg.eigh(fmo_b)

            orb_coefs_a = np.matmul(tmat, evecs_a)
            orb_coefs_b = np.matmul(tmat, evecs_b)

            orb_coefs_a, eigs_a = self.delete_mos(orb_coefs_a, eigs_a)
            orb_coefs_b, eigs_b = self.delete_mos(orb_coefs_b, eigs_b)

            return MolecularOrbitals([orb_coefs_a, orb_coefs_b],
                                     [eigs_a, eigs_b], molorb.unrest)

        return MolecularOrbitals()

    def compute_s2(self, molecule, smat, mol_orbs):
        """
        Computes expectation value of the S**2 operator.

        :param molecule:
            The molecule.
        :param smat:
            The overlap matrix (numpy array).
        :param mol_orbs:
            The molecular orbitals.

        :return:
            Expectation value <S**2>.
        """

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()

        a_b = float(nalpha - nbeta) / 2.0
        s2_exact = a_b * (a_b + 1.0)

        Cocc_a = mol_orbs.alpha_to_numpy()[:, :nalpha]
        Cocc_b = mol_orbs.beta_to_numpy()[:, :nbeta]

        ovl_a_b = np.matmul(Cocc_a.T, np.matmul(smat, Cocc_b))
        s2 = s2_exact + nbeta - np.sum(ovl_a_b**2)

        return s2

    def get_scf_type(self):
        """
        Gets string for spin unrestricted open shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin unrestricted open shell SCF calculation.
        """

        if self.dft:
            return "Spin-Unrestricted Kohn-Sham"
        
        return "Spin-Unrestricted Hartree-Fock"

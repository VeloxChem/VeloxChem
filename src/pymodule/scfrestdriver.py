import numpy as np

from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfRestrictedDriver(ScfDriver):
    """
    Implements spin restricted closed shell SCF method with C2-DIIS and
    two-level C2-DIIS convergence accelerators.
    """

    def __init__(self, comm, ostream):
        """
        Initializes spin restricted closed shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        super().__init__(comm, ostream)

        self.restricted = True

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin restricted closed shell electronic gradient using
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

            dmat = den_mat.alpha_to_numpy(0)
            fmat = fock_mat.to_numpy(0)

            fds = np.matmul(fmat, np.matmul(dmat, smat))

            e_grad = 2.0 * np.linalg.norm(
                np.matmul(tmat.T, np.matmul(fds - fds.T, tmat)))
        else:
            e_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())

        return e_grad

    def comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of spin restricted closed shell density change between
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
            ddmat = diff_mat.alpha_to_numpy(0)

            diff_den = np.linalg.norm(ddmat)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def store_diis_data(self, i, fock_mat, den_mat):
        """
        Stores spin restricted closed shell Fock/Kohn-Sham and density matrices
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

                self.fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self.den_matrices.append(den_mat.alpha_to_numpy(0))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        in OAO basis by applying Lowdin or canonical orthogonalization to AO
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The effective Fock/Kohn-Sham matrix.
        """

        if self.rank == mpi_master():

            if len(self.fock_matrices) == 1:

                return np.copy(self.fock_matrices[0])

            if len(self.fock_matrices) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_error_vectors(self.fock_matrices,
                                               self.den_matrices, ovl_mat,
                                               oao_mat)

                weights = acc_diis.compute_weights()

                return self.get_scaled_fock(weights)

            return fock_mat.alpha_to_numpy(0)

        return None

    def get_scaled_fock(self, weights):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self.fock_matrices[0].shape, dtype=float)

        for w, fmat in zip(weights, self.fock_matrices):

            effmat = effmat + w * fmat

        return effmat

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted closed shell Fock/Kohn-Sham matrix. Overloaded base
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

            fmo = np.matmul(tmat.transpose(), np.matmul(fock_mat, tmat))

            eigs, evecs = np.linalg.eigh(fmo)

            orb_coefs = np.matmul(tmat, evecs)

            orb_coefs, eigs = self.delete_mos(orb_coefs, eigs)

            return MolecularOrbitals([orb_coefs], [eigs], molorb.rest)

        return MolecularOrbitals()

    def get_scf_type(self):
        """
        Gets string for spin restricted closed shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin restricted closed shell SCF calculation.
        """

        return "Spin-Restricted Hartree-Fock"

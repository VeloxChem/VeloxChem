import numpy as np

from .aodensitymatrix import AODensityMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .veloxchemlib import denmat
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfRestrictedOpenDriver(ScfDriver):
    """
    Implements spin restricted open shell SCF
    """

    def __init__(self, comm, ostream):
        """
        Initializes spin restricted closed shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.
        """

        super().__init__(comm, ostream)

        self.restricted = True
        self.restricted_open = True

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin restricted open shell electronic gradient using
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

            e_mat_a = np.matmul(tmat.T, np.matmul(fds_a - fds_a.T, tmat))
            e_mat_b = np.matmul(tmat.T, np.matmul(fds_b - fds_b.T, tmat))

            e_mat = e_mat_a + e_mat_b
            e_mat *= np.sqrt(2)

            e_grad = np.linalg.norm(e_mat)
            max_grad = np.max(np.abs(e_mat))
        else:
            e_grad = 0.0
            max_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        return e_grad, max_grad

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

                    self.fock_matrices_proj.popleft()

                self.fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self.den_matrices.append(den_mat.alpha_to_numpy(0))

                self.fock_matrices_beta.append(fock_mat.beta_to_numpy(0))
                self.den_matrices_beta.append(den_mat.beta_to_numpy(0))

                self.fock_matrices_proj.append(
                    self.get_projected_fock(
                        fock_mat.alpha_to_numpy(0),
                        fock_mat.beta_to_numpy(0),
                        den_mat.alpha_to_numpy(0),
                        den_mat.beta_to_numpy(0),
                        self.scf_tensors['S'],
                    )
                )

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
                return self.fock_matrices_proj[0]

            if len(self.fock_matrices) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_restricted_open_error_vectors(
                    self.fock_matrices, self.fock_matrices_beta,
                    self.den_matrices, self.den_matrices_beta, ovl_mat, oao_mat)

                weights = acc_diis.compute_weights()

                return self.get_scaled_fock(weights)

        return None, None

    def get_scaled_fock(self, weights):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self.fock_matrices_proj[0].shape, dtype=float)

        for w, fmat in zip(weights, self.fock_matrices_proj):
            effmat = effmat + w * fmat

        return effmat

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted projected open shell Fock/Kohn-Sham matrix. Overloaded base
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
            fock_mat = tmat.T @ fock_mat @ tmat

            eigs, evecs = np.linalg.eigh(fock_mat)

            orb_coefs = tmat @ evecs
            orb_coefs, eigs = self.delete_mos(orb_coefs, eigs)

            return MolecularOrbitals([orb_coefs], [eigs], molorb.rest)

        return MolecularOrbitals()

    def get_projected_fock(self, fa, fb, da, db, s):

        naos = len(s)

        f0 = 0.5 * (fa + fb)

        ga = s @ da @ fa - fa @ da @ s
        gb = s @ db @ fb - fb @ db @ s
        g = ga + gb

        inactive = s @ db
        active = s @ (da - db)
        virtual = np.eye(naos) - s @ da

        fcorr = inactive @ (g - f0) @ active.T - active @ (g + f0) @ inactive.T \
            + active @ (g - f0) @ virtual.T - virtual @ (g + f0) @ active.T

        fproj = f0 + fcorr

        return fproj

    def get_scf_type(self):
        """
        Gets string for spin restricted open shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin unrestricted open shell SCF calculation.
        """

        pe_type = " with PE" if self.pe else ""

        if self.dft:
            return "Spin-Restricted Open-Shell Kohn-Sham" + pe_type

        return "Spin-Restricted Open-Shell Hartree-Fock" + pe_type

    def gen_new_density(self, molecule):
        c = self.mol_orbs.alpha_to_numpy()
        na = molecule.number_of_alpha_electrons()
        nb = molecule.number_of_beta_electrons()
        da = c[:, :na] @ c[:, :na].T
        db = c[:, :nb] @ c[:, :nb].T
        ao_density = AODensityMatrix([da, db], denmat.unrest)
        return ao_density

import numpy as np

from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .aodensitymatrix import AODensityMatrix
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfUnrestrictedDriver(ScfDriver):
    """Implements spin unrestricted open shell SCF method (derrived class).

        Implements spin unrestricted open shell SCF method with C2-DIIS and
        two-level C2-DIIS convergence accelerators.
    """

    def __init__(self, comm, ostream):
        """Initializes spin unrestricted open shell SCF driver.

        Initializes spin unrestricted open shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.
        """

        super().__init__(comm, ostream)

    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        """Computes spin unrestricted open shell SCF energy components.

        Computes spin unrestricted open shell SCF energy components: electronic
        energy, kinetic energy, and nuclear potential energy. Overloaded base
        class method.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix (only 2e-part).
        kin_mat
            The kinetic energy matrix.
        npot_mat
            The nuclear potential matrix.
        den_mat
            The density matrix.
        Returns
        -------
            The tuple (electronic energy, kinetic energy, nuclear potential
            energy).
        """

        if self.rank == mpi_master():
            # electronic, kinetic, nuclear energy
            e_ee = fock_mat.get_energy(0, den_mat, 0)
            e_kin = 2.0 * kin_mat.get_energy(den_mat, 0)
            e_en = -2.0 * npot_mat.get_energy(den_mat, 0)
        else:
            e_ee = 0.0
            e_kin = 0.0
            e_en = 0.0

        e_ee = self.comm.bcast(e_ee, root=mpi_master())
        e_kin = self.comm.bcast(e_kin, root=mpi_master())
        e_en = self.comm.bcast(e_en, root=mpi_master())

        return (e_ee, e_kin, e_en)

    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        """Computes full spin unrestricted open shell Fock/Kohn-Sham matrix.

        Computes full spin unrestricted open shell Fock/Kohn-Sham matrix by
        adding to 2e-part of Fock/Kohn-Sham matrix the kinetic energy and
        nuclear potential matrices. Overloaded base class method.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix (2e-part).
        kin_mat
            The kinetic energy matrix.
        npot_mat
            The nuclear potential matrix.
        """

        if self.rank == mpi_master():
            fock_mat.add_hcore(kin_mat, npot_mat, 0)

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """Computes spin unrestricted open shell electronic gradient.

        Computes spin unrestricted open shell electronic gradient using
        Fock/Kohn-Sham matrix. Overloaded base class method.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        ovl_mat
            The overlap matrix.
        den_mat
            The density matrix.
        oao_mat
            The orthogonalization matrix.
        Returns
        -------
            The electronic gradient.
        """

        if self.rank == mpi_master():
            smat = ovl_mat.to_numpy()

            dmat_a = den_mat.alpha_to_numpy(0)
            dmat_b = den_mat.beta_to_numpy(0)

            fmat_a = fock_mat.alpha_to_numpy(0)
            fmat_b = fock_mat.beta_to_numpy(0)

            tmat = oao_mat.to_numpy()

            fds_a = np.matmul(fmat_a, np.matmul(dmat_a, smat))
            sdf_a = np.matmul(smat, np.matmul(dmat_a, fmat_a))

            fds_b = np.matmul(fmat_b, np.matmul(dmat_b, smat))
            sdf_b = np.matmul(smat, np.matmul(dmat_b, fmat_b))

            # fx = np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat))

            e_grad_a = 2.0 * np.linalg.norm(
                np.matmul(tmat.transpose(), np.matmul(fds_a - sdf_a, tmat)))

            e_grad_b = 2.0 * np.linalg.norm(
                np.matmul(tmat.transpose(), np.matmul(fds_b - sdf_b, tmat)))

            e_grad = max(e_grad_a, e_grad_b)
        else:
            e_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())

        return e_grad

    def comp_density_change(self, den_mat, old_den_mat):
        """Computes norm of spin unrestricted open shell density change.

        Computes norm of spin unrestricted open shell density change between
        two density matrices. Overloaded base class method.

        Parameters
        ----------
        den_mat
            The current density matrix.
        old_den_mat
            The previous density matrix.
        Returns
        -------
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
        """Stores spin unrestricted open shell Fock/Kohn-Sham and density
        matrices for current iteration.

        Stores spin unrestricted open shell Fock/Kohn-Sham and density matrices
        for current iteration. Overloaded base class method.

        Parameters
        ----------
        i
            The number of current SCF iteration.
        fock_mat
            The Fock/Kohn-Sham matrix.
        den_mat
            The density matrix.
        """

        if self.rank == mpi_master():

            if not self.skip_iter:

                if len(self.fock_matrices) == self.max_err_vecs:

                    self.fock_matrices.popleft()
                    self.den_matrices.popleft()

                    self.fock_matrices_beta.popleft()
                    self.den_matrices_beta.popleft()

                self.fock_matrices.append(np.copy(fock_mat.alpha_to_numpy(0)))
                self.den_matrices.append(np.copy(den_mat.alpha_to_numpy(0)))

                self.fock_matrices_beta.append(
                    np.copy(fock_mat.beta_to_numpy(0)))
                self.den_matrices_beta.append(np.copy(den_mat.beta_to_numpy(0)))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """Computes effective spin unrestricted open shell Fock/Kohn-Sham
        matrix in OAO basis.

        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        in OAO basis by applying Lowdin or canonical orthogonalization to AO
        Fock/Kohn-Sham matrix. Overloaded base class method.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        ovl_mat
            The overlap matrix.
        oao_mat
            The orthogonalization matrix.
        """

        if self.rank == mpi_master():

            if len(self.fock_matrices) == 1:

                return (np.copy(self.fock_matrices[0]),
                        np.copy(self.fock_matrices_beta[0]))

            if len(self.fock_matrices) > 1:

                acc_diis_a = CTwoDiis()
                acc_diis_b = CTwoDiis()

                acc_diis_a.compute_error_vectors(self.fock_matrices,
                                                 self.den_matrices, ovl_mat,
                                                 oao_mat)
                acc_diis_b.compute_error_vectors(self.fock_matrices_beta,
                                                 self.den_matrices_beta,
                                                 ovl_mat, oao_mat)

                err_a = acc_diis_a.comp_bmatrix()
                err_b = acc_diis_b.comp_bmatrix()

                bmat = 0.5 * (err_a + err_b)
                beigs, bvecs = np.linalg.eigh(bmat)

                weights_a = acc_diis_a.pick_weights(
                    acc_diis_a.norm_bvectors(bvecs))
                weights_b = acc_diis_b.pick_weights(
                    acc_diis_b.norm_bvectors(bvecs))

                return self.get_scaled_fock((weights_a, weights_b))

            return (np.copy(fock_mat.alpha_to_numpy(0)),
                    np.copy(fock_mat.beta_to_numpy(0)))

        return None

    def get_scaled_fock(self, weights):
        """Computes scaled spin unrestricted open shell Fock/Kohn-Sham matrix.

        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        Parameters
        ----------
        weights
            The weights of Fock/Kohn-Sham matrices.
        Returns
        -------
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat_a = np.zeros(self.fock_matrices[0].shape, dtype=float)
        effmat_b = np.zeros(self.fock_matrices[1].shape, dtype=float)

        for ind, (fa, fb) in enumerate(
                zip(self.fock_matrices, self.fock_matrices_beta)):

            effmat_a += weights[0][ind] * fa
            effmat_b += weights[1][ind] * fb

        return effmat_a, effmat_b

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """Generates spin unrestricted molecular orbitals.

        Generates spin unrestricted molecular orbital by diagonalizing
        spin unrestricted open shell Fock/Kohn-Sham matrix. Overloaded base
        class method.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        oao_mat
            The orthogonalization matrix.
        Returns
        -------
            The molecular orbitals.
        """

        if self.rank == mpi_master():

            tmat = oao_mat.to_numpy()

            fmo_a = np.matmul(tmat.transpose(), np.matmul(fock_mat[0], tmat))
            fmo_b = np.matmul(tmat.transpose(), np.matmul(fock_mat[1], tmat))

            eigs_a, evecs_a = np.linalg.eigh(fmo_a)
            eigs_b, evecs_b = np.linalg.eigh(fmo_b)

            orb_coefs_a = np.matmul(tmat, evecs_a)
            orb_coefs_b = np.matmul(tmat, evecs_b)

            orb_coefs_a, eigs_a = self.delete_mos(orb_coefs_a, eigs_a)
            orb_coefs_b, eigs_b = self.delete_mos(orb_coefs_b, eigs_b)

            return MolecularOrbitals([orb_coefs_a, orb_coefs_b],
                                     [eigs_a, eigs_b], molorb.unrest)

        return MolecularOrbitals()

    def gen_new_density(self, molecule):
        """Generates spin unrestricted open shell density matrix.

        Generates spin unrestricted open shell density matrix from current
        spin unrestricted molecular orbitals. Overloaded base class method.

        Parameters
        ----------
        molecule
            The molecule.
        Returns
        -------
            The density matrix.
        """

        if self.rank == mpi_master():

            return self.mol_orbs.get_density(molecule)

        return AODensityMatrix()

    def print_scf_energy(self):
        """Prints SCF energy information to output stream.

        Prints SCF energy information to output stream.

        Parameters
        ----------
        molecule
            The molecule.
        """

        self.ostream.print_header("Spin-Unrestricted Hatree-Fock:".ljust(92))
        self.ostream.print_header("------------------------------".ljust(92))
        self.print_energy_components()

        return

    def get_scf_type(self):
        """Gets string for spin unrestricted open shell SCF calculation.

        Gets string for spin unrestricted open shell SCF calculation.
        Overloaded base class method.

        Returns
        -------
            The string for spin unrestricted open shell SCF calculation.
        """

        return "Spin-Unrestricted Hatree-Fock"

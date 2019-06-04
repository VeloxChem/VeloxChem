import numpy as np

from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .aodensitymatrix import AODensityMatrix
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfRestrictedDriver(ScfDriver):
    """Implements spin restricted closed shell SCF method (derrived class).

        Implements spin restricted closed shell SCF method with C2-DIIS and
        two-level C2-DIIS convergence accelerators.
    """

    def __init__(self, comm, ostream):
        """Initializes spin restricted closed shell SCF driver.

        Initializes spin restricted closed shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.
        """

        super().__init__(comm, ostream)

    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        """Computes spin restricted closed shell SCF energy components.

        Computes spin restricted closed shell SCF energy components: electronic
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
        """Computes full spin restricted closed shell Fock/Kohn-Sham matrix.

        Computes full spin restricted closed shell Fock/Kohn-Sham matrix by
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
        """Computes spin restricted closed shell electronic gradient.

        Computes spin restricted closed shell electronic gradient using
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
            dmat = den_mat.alpha_to_numpy(0)
            fmat = fock_mat.to_numpy(0)
            tmat = oao_mat.to_numpy()

            fa = np.matmul(fmat, np.matmul(dmat, smat))
            fb = np.matmul(smat, np.matmul(dmat, fmat))

            # fx = np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat))

            e_grad = 2.0 * np.linalg.norm(
                np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat)))
        else:
            e_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())

        return e_grad

    def comp_density_change(self, den_mat, old_den_mat):
        """Computes norm of spin restricted closed shell density change.

        Computes norm of spin restricted closed shell density change between
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
            ddmat = diff_mat.alpha_to_numpy(0)

            diff_den = np.linalg.norm(ddmat)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def store_diis_data(self, i, fock_mat, den_mat):
        """Stores spin restricted closed shell Fock/Kohn-Sham and density
        matrices for current iteration.

        Stores spin restricted closed shell Fock/Kohn-Sham and density matrices
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

                self.fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self.den_matrices.append(den_mat.alpha_to_numpy(0))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """Computes effective spin restricted closed shell Fock/Kohn-Sham
        matrix in OAO basis.

        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
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
        """Computes scaled spin restricted closed shell Fock/Kohn-Sham matrix.

        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        Parameters
        ----------
        weights
            The weights of Fock/Kohn-Sham matrices.
        Returns
        -------
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self.fock_matrices[0].shape, dtype=float)

        for w, fmat in zip(weights, self.fock_matrices):

            effmat = effmat + w * fmat

        return effmat

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """Generates spin restricted molecular orbitals.

        Generates spin restricted molecular orbital by diagonalizing
        spin restricted closed shell Fock/Kohn-Sham matrix. Overloaded base
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

            fmo = np.matmul(tmat.transpose(), np.matmul(fock_mat, tmat))

            eigs, evecs = np.linalg.eigh(fmo)

            orb_coefs = np.matmul(tmat, evecs)

            orb_coefs, eigs = self.delete_mos(orb_coefs, eigs)

            return MolecularOrbitals([orb_coefs], [eigs], molorb.rest)

        return MolecularOrbitals()

    def gen_new_density(self, molecule):
        """Generates spin restricted closed shell density matrix.

        Generates spin restricted closed shell density matrix from current
        spin restricted molecular orbitals. Overloaded base class method.

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

        self.ostream.print_header("Spin-Restricted Hatree-Fock:".ljust(92))
        self.ostream.print_header("----------------------------".ljust(92))
        self.print_energy_components()

        return

    def get_scf_type(self):
        """Gets string for spin restricted closed shell SCF calculation.

        Gets string for spin restricted closed shell SCF calculation.
        Overloaded base class method.

        Returns
        -------
            The string for spin restricted closed shell SCF calculation.
        """

        return "Spin-Restricted Hatree-Fock"

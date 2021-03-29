import numpy as np

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .orbitalresponse import OrbitalResponse
from .qqscheme import get_qq_scheme


class TdaOrbitalResponse(OrbitalResponse):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the Tamm-Dancoff Approximation (TDA)
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        super().__init__(comm, ostream)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param rsp_dict:
            The dictionary of response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(rsp_dict, method_dict)

    def compute_rhs(self, molecule, basis, scf_tensors, tda_results, profiler):
        """
        Computes the right-hand side (RHS) of the TDA orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from the converged SCF calculation.
        :param tda_results:
            The results from the converged TDA calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        profiler.start_timer(0, 'RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Write the linear operator for matrix-vector product => in parent class
        # 5) Run the conjugate gradient => in parent class

        if self.rank == mpi_master():

            # 1) Calculate unrelaxed one-particle and transition density matrix
            ovlp = scf_tensors['S']
            mo = scf_tensors['C']

            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            # Take vector of interest and convert to matrix form
            exc_vec = tda_results['eigenvectors'][:, self.n_state_deriv]
            exc_vec = exc_vec.reshape(nocc, nvir).copy()

            # Transform the excitation vectors to the AO basis
            exc_vec_ao = np.linalg.multi_dot([mo_occ, exc_vec, mo_vir.T])

            # Calcuate the unrelaxed one-particle density matrix in MO basis
            dm_oo = -np.matmul(exc_vec, exc_vec.T)
            dm_vv = np.matmul(exc_vec.T, exc_vec)

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = (np.linalg.multi_dot([mo_occ, dm_oo, mo_occ.T]) +
                           np.linalg.multi_dot([mo_vir, dm_vv, mo_vir.T]))

            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix([unrel_dm_ao, exc_vec_ao], denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_ao_rhs.set_fock_type(fockmat.rgenjk, 1)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)
		# TODO: probably at this point we will have to do another calculation
		# for the additional contributions from the xc kernel derivative in DFT
        fock_ao_rhs.reduce_sum(self.rank, self.nodes, self.comm)

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            fock_ao_rhs_0 = fock_ao_rhs.alpha_to_numpy(0)
            fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1)

            fmo_rhs_0 = np.linalg.multi_dot(
                [mo_occ.T, 0.5 * fock_ao_rhs_0, mo_vir])

            sdp_pds = (
                np.linalg.multi_dot([ovlp, exc_vec_ao, 0.5 * fock_ao_rhs_1.T]) -
                np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, exc_vec_ao, ovlp]))

            rhs_mo = fmo_rhs_0 + np.linalg.multi_dot(
                [mo_occ.T, sdp_pds, mo_vir])

        profiler.stop_timer(0, 'RHS')

        if self.rank == mpi_master():
            return {
                'rhs_mo': rhs_mo,
                'dm_oo': dm_oo,
                'dm_vv': dm_vv,
                'unrel_dm_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
            }
        else:
            return {}

    def compute_omega(self, ovlp, mo_occ, mo_vir, epsilon_dm_ao, tda_results,
                      fock_ao_rhs, fock_lambda):
        """
        Calculates the TDA Lagrange multipliers for the overlap matrix.

        :param ovlp:
            The overlap matrix.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.
        :param epsilon_dm_ao:
            The energy-weighted relaxed density matrix.
        :param tda_results:
            The results from the TDA calculation.
        :param fock_ao_rhs:
            The AOFockMatrix from the right-hand side of the orbital response eq.
        :param fock_lambda:
            The Fock matrix from Lagrange multipliers.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        # Get the excitation vector of interest and transform it to AO
        nocc = mo_occ.shape[1]
        nvir = mo_vir.shape[1]
        exc_vec = tda_results['eigenvectors'][:, self.n_state_deriv]
        exc_vec = exc_vec.reshape(nocc, nvir).copy()
        exc_vec_ao = np.linalg.multi_dot([mo_occ, exc_vec, mo_vir.T])

        # The density matrix; only alpha block;
        # Only works for the restricted case
        D_occ = np.matmul(mo_occ, mo_occ.T)
        D_vir = np.matmul(mo_vir, mo_vir.T)

        # Because the excitation vector is not symmetric,
        # we need both the matrix (OO block in omega, and probably VO)
        # and its transpose (VV, OV blocks)
        # this comes from the transformation of the 2PDM contribution
        # from MO to AO basis
        fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1)
        Ft = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, exc_vec_ao, ovlp])
        F = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, exc_vec_ao.T, ovlp.T])

        # Compute the contributions from the 2PDM and the relaxed 1PDM
        # to the omega Lagrange multipliers:
        fmat = (fock_lambda.alpha_to_numpy(0) +
                fock_lambda.alpha_to_numpy(0).T +
                0.5 * fock_ao_rhs.alpha_to_numpy(0))

        omega_1pdm_2pdm_contribs = (np.linalg.multi_dot([D_occ, F, D_occ]) +
                                    np.linalg.multi_dot([D_occ, Ft, D_vir]) +
                                    np.linalg.multi_dot([D_occ, Ft, D_vir]).T +
                                    np.linalg.multi_dot([D_vir, Ft, D_vir]) +
                                    np.linalg.multi_dot([D_occ, fmat, D_occ]))

        omega = -epsilon_dm_ao - omega_1pdm_2pdm_contribs

        return omega

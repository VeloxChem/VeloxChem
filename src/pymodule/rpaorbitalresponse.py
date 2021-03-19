import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .profiler import Profiler
from .orbitalresponse import OrbitalResponse
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme


class RpaOrbitalResponse(OrbitalResponse):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the random phase approximation (RPA)
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - n_state_deriv: The number of the excited state of interest.
        - solver: The linear equations solver.
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

    def compute(self, molecule, basis, scf_tensors, rpa_results):
        """
        Performs orbital response Lagrange multipliers
        calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
        :param rpa_results:
            The results from converged RPA calculation.

        :return:
            A dictionary containing the Lagrange multipliers and relaxed
            one-particle density.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'OrbitalResponse: not implemented for unrestricted case')

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
            ea = scf_tensors['E']

            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            eocc = ea[:nocc]
            evir = ea[nocc:]

            # Take vector of interest and convert to matrix form
            exc_vec = rpa_results['eigenvectors'][:nocc * nvir, self.n_state_deriv]
            deexc_vec = rpa_results['eigenvectors'][nocc * nvir:, self.n_state_deriv]
            exc_vec = exc_vec.reshape(nocc, nvir).copy()
            deexc_vec = deexc_vec.reshape(nocc, nvir).copy()

            # TODO: remove factor once RPA vectors are properly normalized
            #exc_vec *= np.sqrt(2)
            #deexc_vec *= np.sqrt(2)

            # Construct plus/minus combinations of excitation and de-excitation part
            xpy = exc_vec + deexc_vec
            xmy = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            xpy_ao = np.linalg.multi_dot([mo_occ, xpy, mo_vir.T])
            xmy_ao = np.linalg.multi_dot([mo_occ, xmy, mo_vir.T])

            # Calcuate the unrelaxed one-particle density matrix in MO basis
            dm_oo = -0.5 * (np.matmul(xpy, xpy.T) + np.matmul(xmy, xmy.T))
            dm_vv =  0.5 * (np.matmul(xpy.T, xpy) + np.matmul(xmy.T, xmy))

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = (np.linalg.multi_dot([mo_occ, dm_oo, mo_occ.T]) +
                           np.linalg.multi_dot([mo_vir, dm_vv, mo_vir.T]))

            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix([unrel_dm_ao, xpy_ao, xmy_ao], denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_ao_rhs.set_fock_type(fockmat.rgenjk, 1)
        fock_ao_rhs.set_fock_type(fockmat.rgenjk, 2)

        # TODO: make sure that ERI settings are consistent with response solver

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)
        fock_ao_rhs.reduce_sum(self.rank, self.nodes, self.comm)

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            fock_ao_rhs_1dm = fock_ao_rhs.alpha_to_numpy(0)
            fock_ao_rhs_xpy = fock_ao_rhs.alpha_to_numpy(1)
            fock_ao_rhs_xmy = fock_ao_rhs.alpha_to_numpy(2)

            fmo_rhs_1dm = np.linalg.multi_dot(
                [mo_occ.T, 0.5 * fock_ao_rhs_1dm, mo_vir])

            # TODO: factor out overlap matrix?
            sdp_pds = 0.5 * (
                np.linalg.multi_dot([ovlp, xpy_ao, 0.5 * fock_ao_rhs_xpy.T]) +
                np.linalg.multi_dot([ovlp, xmy_ao, 0.5 * fock_ao_rhs_xmy.T]) -
                np.linalg.multi_dot([0.5 * fock_ao_rhs_xpy.T, xpy_ao, ovlp]) -
                np.linalg.multi_dot([0.5 * fock_ao_rhs_xmy.T, xmy_ao, ovlp]) -
                np.linalg.multi_dot([ovlp, xpy_ao, 0.5 * fock_ao_rhs_xpy]) +
                np.linalg.multi_dot([ovlp, xmy_ao, 0.5 * fock_ao_rhs_xmy]) +
                np.linalg.multi_dot([0.5 * fock_ao_rhs_xpy, xpy_ao, ovlp]) -
                np.linalg.multi_dot([0.5 * fock_ao_rhs_xmy, xmy_ao, ovlp])
                )

            rhs_mo = fmo_rhs_1dm + np.linalg.multi_dot(
                [mo_occ.T, sdp_pds, mo_vir])
        else:
            rhs_mo = None
        rhs_mo = self.comm.bcast(rhs_mo, root=mpi_master())

        profiler.stop_timer(0, 'RHS')

        # Calculate the lambda multipliers in the parent class
        lambda_multipliers = self.compute_lambda(molecule, basis, scf_tensors,
                                                 rhs_mo, profiler)

        profiler.start_timer(0, 'omega')

        # Calculate the overlap matrix multipliers
        if self.rank == mpi_master():

            # 1. compute an energy-weighted density matrix
            # (needed to compute omega)
            eo_diag = np.diag(eocc)
            ev_diag = np.diag(evir)

            epsilon_dm_ao = np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, 0.5 * dm_oo) + eo_diag, mo_occ.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_vir, np.matmul(ev_diag, 0.5 * dm_vv), mo_vir.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, lambda_multipliers), mo_vir.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, lambda_multipliers), mo_vir.T]).T

            # 2. compute the omega multipliers in AO basis:
            lambda_ao = np.linalg.multi_dot(
                [mo_occ, lambda_multipliers, mo_vir.T])
            ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)
        else:
            ao_density_lambda = AODensityMatrix()
        ao_density_lambda.broadcast(self.rank, self.comm)

        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_lambda.set_fock_type(fockmat.rgenjk, 0)

        eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                        screening)
        fock_lambda.reduce_sum(self.rank, self.nodes, self.comm)

        if self.rank == mpi_master():
            omega_ao = self.compute_omega(ovlp, mo_occ, mo_vir, epsilon_dm_ao,
                                          xpy_ao, xmy_ao, fock_ao_rhs, fock_lambda)

        self.ostream.print_blank()
        self.ostream.flush()

        profiler.stop_timer(0, 'omega')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of CG')
        profiler.print_memory_usage(self.ostream)

        if self.rank == mpi_master():
            if self.is_converged:
                # Calculate the relaxed one-particle density matrix
                # Factor 4: (ov + vo)*(alpha + beta)
                rel_dm_ao = unrel_dm_ao + 4 * lambda_ao
                return {
                    'lambda_ao': lambda_ao,
                    'omega_ao': omega_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'rel_dm_ao': rel_dm_ao,
                }
            else:
                # return only unrelaxed density matrix if not converged
                return {'unrel_dm_ao': unrel_dm_ao}
        else:
            return {}

    def compute_omega(self, ovlp, mo_occ, mo_vir, epsilon_dm_ao,
                      xpy_ao, xmy_ao, fock_ao_rhs, fock_lambda):
        """
        Calculates the Lagrange multipliers for the overlap matrix.

        :param ovlp:
            The overlap matrix.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.
        :param epsilon_dm_ao:
            The energy-weighted relaxed density matrix.
        :param xpy_ao:
            The sum of excitation and deexcitation vectors in AO basis.
        :param xmy_ao:
            The difference of excitation and deexcitation vectors in AO basis.
        :param fock_ao_rhs:
            The AOFockMatrix from the right-hand side of the orbital response eq.
        :param fock_lambda:
            The Fock matrix from Lagrange multipliers.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        # The density matrix; only alpha block;
        # Only works for the restricted case
        D_occ = np.matmul(mo_occ, mo_occ.T)
        D_vir = np.matmul(mo_vir, mo_vir.T)

        # Because the excitation vector is not symmetric,
        # we need both the matrix (OO block in omega, and probably VO)
        # and its transpose (VV, OV blocks)
        # this comes from the transformation of the 2PDM contribution
        # from MO to AO basis
        fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1) # xpy
        fock_ao_rhs_2 = fock_ao_rhs.alpha_to_numpy(2) # xmy

        Fp1_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao, ovlp.T])
        Fm1_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao, ovlp.T])
        Fp2_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao, ovlp.T])
        Fm2_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao, ovlp.T])
        #Fp1_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao, ovlp.T])
        #Fm1_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao, ovlp.T])
        #Fp2_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao, ovlp.T])
        #Fm2_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao, ovlp.T])
        Fp1_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao.T, ovlp.T])
        Fm1_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao.T, ovlp.T])
        Fp2_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao.T, ovlp.T])
        Fm2_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao.T, ovlp.T])
        # We see that:
        # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
        # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

        #Ft = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao, ovlp])
        #F = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao.T, ovlp.T])

        # Compute the contributions from the 2PDM and the relaxed 1PDM
        # to the omega Lagrange multipliers:
        fmat = (fock_lambda.alpha_to_numpy(0) +
                fock_lambda.alpha_to_numpy(0).T +
                0.5 * fock_ao_rhs.alpha_to_numpy(0))

        omega_1pdm_2pdm_contribs = 0.5 * (np.linalg.multi_dot([D_vir, Fp1_vv + Fm1_vv
                                      - Fp2_vv + Fm2_vv, D_vir]) +
                                    np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv
                                      - Fp2_vv + Fm2_vv, D_vir]) +
                                    np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv
                                      - Fp2_vv + Fm2_vv, D_vir]).T +
                                    np.linalg.multi_dot([D_occ, Fp1_oo + Fm1_oo
                                      - Fp2_vv + Fm2_vv, D_occ]) +
                                    2 * np.linalg.multi_dot([D_occ, fmat, D_occ])
                                    )


        #omega_1pdm_2pdm_contribs = (np.linalg.multi_dot([D_occ, F, D_occ]) +
        #                            np.linalg.multi_dot([D_occ, Ft, D_vir]) +
        #                            np.linalg.multi_dot([D_occ, Ft, D_vir]).T +
        #                            np.linalg.multi_dot([D_vir, Ft, D_vir]) +
        #                            np.linalg.multi_dot([D_occ, fmat, D_occ]))

        omega = -epsilon_dm_ao - omega_1pdm_2pdm_contribs

        return omega

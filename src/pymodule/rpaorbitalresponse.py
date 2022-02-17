import numpy as np

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .orbitalresponse import OrbitalResponse
from .qqscheme import get_qq_scheme


class RpaOrbitalResponse(OrbitalResponse):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the random phase approximation (RPA)
    level of theory.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes orbital response computation driver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        super().__init__(comm, ostream)

    def update_settings(self, orbrsp_dict, rsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param rsp_dict:
            The dictionary of response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(orbrsp_dict, rsp_dict, method_dict)

    def compute_rhs(self, molecule, basis, scf_tensors, rpa_results, dft_dict, profiler):
        """
        Computes the right-hand side (RHS) of the RPA orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
        :param rpa_results:
            The results from converged RPA calculation.
        :param dft_dict:
            The dictionary containing DFT information.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        profiler.start_timer('RHS')

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
            exc_vec = rpa_results['eigenvectors'][:nocc * nvir,
                                                  self.state_deriv_index]
            deexc_vec = rpa_results['eigenvectors'][nocc * nvir:,
                                                    self.state_deriv_index]
            exc_vec = exc_vec.reshape(nocc, nvir).copy()
            deexc_vec = deexc_vec.reshape(nocc, nvir).copy()

            # Construct plus/minus combinations of excitation and de-excitation part
            xpy = exc_vec + deexc_vec
            xmy = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            xpy_ao = np.linalg.multi_dot([mo_occ, xpy, mo_vir.T])
            xmy_ao = np.linalg.multi_dot([mo_occ, xmy, mo_vir.T])

            # Calcuate the unrelaxed one-particle density matrix in MO basis
            dm_oo = -0.5 * (np.matmul(xpy, xpy.T) + np.matmul(xmy, xmy.T))
            dm_vv = 0.5 * (np.matmul(xpy.T, xpy) + np.matmul(xmy.T, xmy))

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = (np.linalg.multi_dot([mo_occ, dm_oo, mo_occ.T]) +
                           np.linalg.multi_dot([mo_vir, dm_vv, mo_vir.T]))

            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix([unrel_dm_ao, xpy_ao, xmy_ao],
                                        denmat.rest)
            if self.dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                perturbed_dm_ao = AODensityMatrix([xmy_ao, 0*xpy_ao, xmy_ao, 0*xmy_ao],
                                                   denmat.rest)

                # corresponds to rho^{omega_b,omega_c} in quadratic response,
                # which is zero for TDDFT orbital response
                zero_dm_ao = AODensityMatrix([0*xpy_ao, 0*xpy_ao],
                                              denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()
            if self.dft:
                perturbed_dm_ao = AODensityMatrix()
                zero_dm_ao =  AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_ao_rhs.set_fock_type(fockmat.rgenjk, 1)
        fock_ao_rhs.set_fock_type(fockmat.rgenjk, 2)
        if self.dft:
            perturbed_dm_ao.broadcast(self.rank, self.comm)
            zero_dm_ao.broadcast(self.rank, self.comm)
            # Fock matrix for computing the TDDFT E[3] term g^xc
            fock_gxc_ao = AOFockMatrix(zero_dm_ao)
            if self.xcfun.is_hybrid():
                fact_xc = self.xcfun.get_frac_exact_exchange()
                for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                    fock_ao_rhs.set_scale_factor(fact_xc, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_scale_factor(fact_xc, ifock)
                fock_ao_rhs.set_fock_type(fockmat.restjkx, 0)
                fock_ao_rhs.set_fock_type(fockmat.rgenjkx, 1)
                fock_ao_rhs.set_fock_type(fockmat.rgenjkx, 2)
                fock_gxc_ao.set_fock_type(fockmat.rgenjkx, 0)
                fock_gxc_ao.set_fock_type(fockmat.rgenjkx, 1)
            else:
                fock_ao_rhs.set_fock_type(fockmat.restj, 0)
                fock_ao_rhs.set_fock_type(fockmat.rgenj, 1)
                fock_ao_rhs.set_fock_type(fockmat.rgenj, 2)
                fock_gxc_ao.set_fock_type(fockmat.rgenj, 0)
                fock_gxc_ao.set_fock_type(fockmat.rgenj, 1)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)
        if self.dft:
            if not self.xcfun.is_hybrid():
                for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                    fock_ao_rhs.scale(2.0, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.scale(2.0, ifock)
            xc_drv = XCIntegrator(self.comm)
            molgrid.distribute(self.rank, self.nodes, self.comm)
            # Linear response routine for f^xc
            xc_drv.integrate(fock_ao_rhs, dm_ao_rhs, gs_density,
                             molecule, basis, molgrid, self.xcfun.get_func_label())
            # Quadratic response routine for TDDFT E[3] term g^xc
            xc_drv.integrate(fock_gxc_ao, perturbed_dm_ao, zero_dm_ao,
                             gs_density, molecule, basis, molgrid,
                             self.xcfun.get_func_label(), "quadratic")

            fock_gxc_ao.reduce_sum(self.rank, self.nodes, self.comm)

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
                np.linalg.multi_dot([0.5 * fock_ao_rhs_xmy, xmy_ao, ovlp]))

            rhs_mo = fmo_rhs_1dm + np.linalg.multi_dot(
                [mo_occ.T, sdp_pds, mo_vir])

            # Add TDDFT E[3] contribution to the RHS:
            if self.dft:
                gxc_ao = fock_gxc_ao.alpha_to_numpy(0)
                gxc_mo = np.linalg.multi_dot([mo_occ.T, gxc_ao, mo_vir])
                rhs_mo += 0.25 * gxc_mo

        profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'rhs_mo': rhs_mo,
                'dm_oo': dm_oo,
                'dm_vv': dm_vv,
                'xpy_ao': xpy_ao,
                'xmy_ao': xmy_ao,
                'unrel_dm_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
            }
        else:
            return {}

    def compute_omega(self, ovlp, mo_occ, mo_vir, epsilon_dm_ao, rpa_results,
                      fock_ao_rhs, fock_lambda):
        """
        Calculates the RPA Lagrange multipliers for the overlap matrix.

        :param ovlp:
            The overlap matrix.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.
        :param epsilon_dm_ao:
            The energy-weighted relaxed density matrix.
        :param rpa_results:
            The results from the RPA calculation.
        :param fock_ao_rhs:
            The AOFockMatrix from the right-hand side of the orbital response eq.
        :param fock_lambda:
            The Fock matrix from Lagrange multipliers.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        # Get the excitation and deexcitation vector of interest,
        # construct plus/minus combinations and transform them to AO
        nocc = mo_occ.shape[1]
        nvir = mo_vir.shape[1]
        exc_vec = rpa_results['eigenvectors'][:nocc * nvir, self.state_deriv_index]
        deexc_vec = rpa_results['eigenvectors'][nocc * nvir:,
                                                self.state_deriv_index]
        exc_vec = exc_vec.reshape(nocc, nvir).copy()
        deexc_vec = deexc_vec.reshape(nocc, nvir).copy()
        xpy = exc_vec + deexc_vec
        xmy = exc_vec - deexc_vec
        xpy_ao = np.linalg.multi_dot([mo_occ, xpy, mo_vir.T])
        xmy_ao = np.linalg.multi_dot([mo_occ, xmy, mo_vir.T])

        # Transform the vectors to the AO basis
        # The density matrix; only alpha block;
        # Only works for the restricted case
        # (since scf_tensors['C'] only gives alpha block...)
        D_occ = np.matmul(mo_occ, mo_occ.T)
        D_vir = np.matmul(mo_vir, mo_vir.T)

        # Because the excitation vector is not symmetric,
        # we need both the matrix (OO block in omega, and probably VO)
        # and its transpose (VV, OV blocks)
        # this comes from the transformation of the 2PDM contribution
        # from MO to AO basis
        fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1)  # xpy
        fock_ao_rhs_2 = fock_ao_rhs.alpha_to_numpy(2)  # xmy

        Fp1_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao, ovlp.T])
        Fm1_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao, ovlp.T])
        Fp2_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao, ovlp.T])
        Fm2_vv = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao, ovlp.T])
        # Fp1_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao, ovlp.T])
        # Fm1_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao, ovlp.T])
        # Fp2_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao, ovlp.T])
        # Fm2_ov = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao, ovlp.T])
        Fp1_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, xpy_ao.T, ovlp.T])
        Fm1_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_2, xmy_ao.T, ovlp.T])
        Fp2_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, xpy_ao.T, ovlp.T])
        Fm2_oo = np.linalg.multi_dot([0.5 * fock_ao_rhs_2.T, xmy_ao.T, ovlp.T])
        # We see that:
        # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
        # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

        # Compute the contributions from the 2PDM and the relaxed 1PDM
        # to the omega Lagrange multipliers:
        fmat = (fock_lambda.alpha_to_numpy(0) +
                fock_lambda.alpha_to_numpy(0).T +
                0.5 * fock_ao_rhs.alpha_to_numpy(0))

        omega_1pdm_2pdm_contribs = 0.5 * (
            np.linalg.multi_dot([D_vir, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir])
          + np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir])
          + np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir]).T
          + np.linalg.multi_dot([D_occ, Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ])
          + 2 * np.linalg.multi_dot([D_occ, fmat, D_occ]))

        omega = -epsilon_dm_ao - omega_1pdm_2pdm_contribs

        return omega

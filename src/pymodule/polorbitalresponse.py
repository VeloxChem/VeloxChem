import numpy as np
import time as tm

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .inputparser import parse_input
from .cphfsolver import CphfSolver
from .sanitychecks import polgrad_sanity_check


class PolOrbitalResponse(CphfSolver):
    """
    Implements orbital response Lagrange multipliers computation
    for the polarizability gradient.

    Instance variables
        - frequencies: The sequence of  frequencies for which the
            polarizability gradient is computed.
        - vector_components: The components of the response vectors
            corresponding to the operator components in linear response.
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

        self.flag = 'Polarizability Orbital Response'

        self.is_complex = False

        self.frequencies = (0,)
        self.vector_components = 'xyz'
        self.cphf_results = None

        self._input_keywords['orbitalresponse'].update({
            'vector_components':
                ('str_lower', 'Cartesian components of operator'),
            'frequencies': ('seq_range', 'frequencies'),
            'is_complex': ('bool', 'whether the polarizability is complex'),
        })

    def update_settings(self, orbrsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(orbrsp_dict, method_dict)

        orbrsp_keywords = {
            key: val[0]
            for key, val in self._input_keywords['orbitalresponse'].items()
        }

        parse_input(self, orbrsp_keywords, orbrsp_dict)

    def get_full_solution_vector(self, solution):
        """ Gets a full solution vector from a general distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector
        """

        if self.is_complex:
            return self.get_full_solution_vector_complex(solution)
        else:
            return self.get_full_solution_vector_real(solution)

    @staticmethod
    def get_full_solution_vector_complex(solution):
        """
        Gets a full complex solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The real and imaginary parts of the full solution vector.
        """
        x_realger = solution.get_full_vector(0)
        x_realung = solution.get_full_vector(1)
        x_imagung = solution.get_full_vector(2)
        x_imagger = solution.get_full_vector(3)

        if solution.rank == mpi_master():
            x_real = np.hstack((x_realger, x_realger)) + np.hstack(
                (x_realung, -x_realung))
            x_imag = np.hstack((x_imagung, -x_imagung)) + np.hstack(
                (x_imagger, x_imagger))
            return x_real + 1j * x_imag
        else:
            return None

    @staticmethod
    def get_full_solution_vector_real(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None

    def compute_rhs(self, molecule, basis, scf_tensors, lr_results):

        if self.is_complex:
            return self.compute_rhs_complex(molecule, basis, scf_tensors,
                                            lr_results)
        else:
            return self.compute_rhs_real(molecule, basis, scf_tensors,
                                         lr_results)

    def compute_rhs_complex(self, molecule, basis, scf_tensors, lr_results):
        """
        Computes the complex right-hand side (RHS) of the polarizability
        orbital response equation including the necessary density matrices
        using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
        :param lr_results:
            The results from converged CPP calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Run the solver => in parent class

        loop_start_time = tm.time()

        orbrsp_rhs = {}
        for f, w in enumerate(self.frequencies):

            if self.rank == mpi_master():
                self.ostream.print_info(
                    'Building RHS for w = {:4.3f}'.format(w))
                self.ostream.flush()

                # 1) Calculate unrelaxed one-particle and transition density matrix
                ovlp = scf_tensors['S']
                mo = scf_tensors['C']  # only alpha part

                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]

                # TODO: do we keep this factor like that?
                sqrt2 = np.sqrt(2.0)

                # Check if response vectors exist for desired frequency of gradient
                polgrad_sanity_check(self, self.flag, lr_results)

            # TODO: make get_full_solution_vector directly available from the
            # parent class (i.e. lrsolver)?
            full_vec = ([
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ])
            self.ostream.flush()

            if self.rank == mpi_master():
                # Save the number of vector components
                dof = len(self.vector_components)

                # Extract the excitation and de-excitation components
                # from the full solution vector.
                exc_vec = (1 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                # Construct plus/minus combinations of excitation and
                # de-excitation part
                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec
                self.ostream.flush()

                # Transform the vectors to the AO basis
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(x_plus_y.shape[0])
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(x_minus_y.shape[0])
                ])

                # Turn them into a list (for AODensityMatrix)
                xpmy_ao_list_real = list(np.array(x_plus_y_ao.real)) + list(
                    np.array(x_minus_y_ao.real))
                xpmy_ao_list_imag = list(np.array(x_plus_y_ao.imag)) + list(
                    np.array(x_minus_y_ao.imag))

                # Calculate the symmetrized unrelaxed one-particle density matrix
                # in MO basis
                dm_oo = np.zeros((dof, dof, nocc, nocc), dtype=np.complex_)
                dm_vv = np.zeros((dof, dof, nvir, nvir), dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        dm_vv[x, y] = 0.25 * (
                            np.linalg.multi_dot([x_plus_y[x].T, x_plus_y[y]]).T
                            + np.linalg.multi_dot(
                                [x_minus_y[x].T, x_minus_y[y]]).T +
                            np.linalg.multi_dot([x_plus_y[x].T, x_plus_y[y]]) +
                            np.linalg.multi_dot([x_minus_y[x].T, x_minus_y[y]]))

                        dm_oo[x, y] = -0.25 * (
                            np.linalg.multi_dot([x_plus_y[x], x_plus_y[y].T]) +
                            np.linalg.multi_dot([x_minus_y[x], x_minus_y[y].T])
                            + np.linalg.multi_dot([x_plus_y[x], x_plus_y[y].T
                                                  ]).T +
                            np.linalg.multi_dot([x_minus_y[x], x_minus_y[y].T
                                                ]).T)

                # Transform unrelaxed one-particle density matrix to
                # AO basis and create a list
                unrel_dm_ao = np.zeros((dof, dof, nao, nao), dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        unrel_dm_ao[x,
                                    y] = (np.linalg.multi_dot(
                                        [mo_occ, dm_oo[x, y], mo_occ.T]) +
                                          np.linalg.multi_dot(
                                              [mo_vir, dm_vv[x, y], mo_vir.T]))

                dm_ao_list_real = list(
                    np.array(unrel_dm_ao.real).reshape(dof**2, nao, nao))
                dm_ao_list_imag = list(
                    np.array(unrel_dm_ao.imag).reshape(dof**2, nao, nao))

                # 2) Construct the right-hand side
                dm_ao_rhs_real = AODensityMatrix(
                    dm_ao_list_real + xpmy_ao_list_real, denmat.rest)
                dm_ao_rhs_imag = AODensityMatrix(
                    dm_ao_list_imag + xpmy_ao_list_imag, denmat.rest)

                if self._dft:
                    # FIXME split vectors in Re/Im for list
                    # 3) Construct density matrices for E[3] term:
                    # XCIntegrator expects a DM with real and imaginary part,
                    # so we set the imaginary part to zero.
                    # Create lists with the corresponding vector components
                    perturbed_dm_ao_list = []
                    zero_dm_ao_list = []
                    # TODO: only upper triangular matrix and transpose?
                    for x in range(dof):
                        for y in range(dof):
                            perturbed_dm_ao_list.extend([
                                x_minus_y_ao[x], 0 * x_minus_y_ao[x],
                                x_minus_y_ao[y], 0 * x_minus_y_ao[y]
                            ])
                            zero_dm_ao_list.extend(
                                [0 * x_minus_y_ao[x], 0 * x_minus_y_ao[y]])

                    perturbed_dm_ao = AODensityMatrix(perturbed_dm_ao_list,
                                                      denmat.rest)

                    # corresponds to rho^{omega_b,omega_c} in quadratic response,
                    # which is zero for orbital response
                    zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)
            else:
                dof = None
                dm_ao_rhs_real = AODensityMatrix()
                dm_ao_rhs_imag = AODensityMatrix()
                if self._dft:
                    perturbed_dm_ao = AODensityMatrix()
                    zero_dm_ao = AODensityMatrix()

            dof = self.comm.bcast(dof, root=mpi_master())

            dm_ao_rhs_real.broadcast(self.rank, self.comm)
            dm_ao_rhs_imag.broadcast(self.rank, self.comm)

            molgrid = dft_dict['molgrid']
            gs_density = dft_dict['gs_density']

            # Fock matrices with corresponding type
            fock_ao_rhs_real = AOFockMatrix(dm_ao_rhs_real)
            fock_ao_rhs_imag = AOFockMatrix(dm_ao_rhs_imag)
            # Set the vector-related components to general Fock matrix
            # (not 1PDM part)
            for ifock in range(dof**2, dof**2 + 2 * dof):
                fock_ao_rhs_real.set_fock_type(fockmat.rgenjk, ifock)
                fock_ao_rhs_imag.set_fock_type(fockmat.rgenjk, ifock)
            if self._dft:
                perturbed_dm_ao.broadcast(self.rank, self.comm)
                zero_dm_ao.broadcast(self.rank, self.comm)
                fock_gxc_ao = AOFockMatrix(zero_dm_ao)
                if self.xcfun.is_hybrid():
                    fact_xc = self.xcfun.get_frac_exact_exchange()
                    for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                        fock_ao_rhs.set_scale_factor(fact_xc, ifock)
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.set_scale_factor(fact_xc, ifock)
                        fock_gxc_ao.set_fock_type(fockmat.rgenjkx, ifock)
                    for ifock in range(dof**2):
                        fock_ao_rhs.set_fock_type(fockmat.restjkx, ifock)
                    for ifock in range(dof**2, dof**2 + 2 * dof):
                        fock_ao_rhs.set_fock_type(fockmat.rgenjkx, ifock)
                else:
                    for ifock in range(dof**2):
                        fock_ao_rhs.set_fock_type(fockmat.restj, ifock)
                    for ifock in range(dof**2, dof**2 + 2 * dof):
                        fock_ao_rhs.set_fock_type(fockmat.rgenj, ifock)
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.set_fock_type(fockmat.rgenj, ifock)
            else:
                fock_gxc_ao = None  # None if not DFT

            if self._dft:
                if not self.xcfun.is_hybrid():
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.scale(2.0, ifock)
                xc_drv = XCIntegrator(self.comm)
                xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis,
                                          perturbed_dm_ao, zero_dm_ao,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")

                fock_gxc_ao.reduce_sum(self.rank, self.nodes, self.comm)

            self._comp_lr_fock(fock_ao_rhs_real, dm_ao_rhs_real, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)
            self._comp_lr_fock(fock_ao_rhs_imag, dm_ao_rhs_imag, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)

            # Calculate the RHS and transform it to the MO basis
            if self.rank == mpi_master():
                # extract the 1PDM contributions
                fock_ao_rhs_1dm = np.zeros((dof**2, nao, nao),
                                           dtype=np.complex_)
                fock_ao_rhs_1dm_real = np.zeros((dof**2, nao, nao))
                fock_ao_rhs_1dm_imag = np.zeros((dof**2, nao, nao))
                for i in range(dof**2):
                    fock_ao_rhs_1dm_real[i] = fock_ao_rhs_real.alpha_to_numpy(i)
                    fock_ao_rhs_1dm_imag[i] = fock_ao_rhs_imag.alpha_to_numpy(i)
                # combine to complex array
                fock_ao_rhs_1dm = fock_ao_rhs_1dm_real + 1j * fock_ao_rhs_1dm_imag

                # Transform to MO basis
                fock_mo_rhs_1dm = np.array([
                    np.linalg.multi_dot([mo_occ.T, fock_ao_rhs_1dm[x], mo_vir])
                    for x in range(dof**2)
                ])

                # extract the x_plus_y and x_minus_y contributions
                # TODO: extract all Fock matrices at the same time?
                fock_ao_rhs_x_plus_y_real = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y_real = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_plus_y_imag = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y_imag = np.zeros((dof, nao, nao))
                for i in range(dof):
                    fock_ao_rhs_x_plus_y_real[
                        i] = fock_ao_rhs_real.alpha_to_numpy(dof**2 + i)
                    fock_ao_rhs_x_minus_y_real[
                        i] = fock_ao_rhs_real.alpha_to_numpy(dof**2 + dof + i)
                    fock_ao_rhs_x_plus_y_imag[
                        i] = fock_ao_rhs_imag.alpha_to_numpy(dof**2 + i)
                    fock_ao_rhs_x_minus_y_imag[
                        i] = fock_ao_rhs_imag.alpha_to_numpy(dof**2 + dof + i)
                # combine to complex
                fock_ao_rhs_x_plus_y = (fock_ao_rhs_x_plus_y_real +
                                        1j * fock_ao_rhs_x_plus_y_imag)
                fock_ao_rhs_x_minus_y = (fock_ao_rhs_x_minus_y_real +
                                         1j * fock_ao_rhs_x_minus_y_imag)

                fock_mo_rhs_2dm = np.zeros((dof, dof, nocc, nvir),
                                           dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x], x_plus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x], x_minus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y], x_plus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y], x_minus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x].T, x_plus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x].T, x_minus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y].T, x_plus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y].T, x_minus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x], x_plus_y_ao[y].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x], x_minus_y_ao[y].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y], x_plus_y_ao[x].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y], x_minus_y_ao[x].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x].T, x_plus_y_ao[y].T])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x].T, x_minus_y_ao[y].T])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y].T, x_plus_y_ao[x].T])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y].T, x_minus_y_ao[x].T])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                fock_mo_rhs_2dm = 0.25 * fock_mo_rhs_2dm.reshape(
                    dof**2, nocc, nvir)

                # Calculate the dipole contributions to the RHS:
                # Dipole integrals in AO basis
                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_mats = dipole_drv.compute(molecule, basis)
                dipole_ints_ao = np.zeros((dof, nao, nao))
                k = 0
                if 'x' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                    k += 1
                if 'y' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                    k += 1
                if 'z' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.z_to_numpy()

                # Transform them to MO basis (oo and vv blocks only)
                dipole_ints_oo = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
                    for x in range(dof)
                ])
                dipole_ints_vv = np.array([
                    np.linalg.multi_dot([mo_vir.T, dipole_ints_ao[x], mo_vir])
                    for x in range(dof)
                ])

                # Contract with vectors to get dipole contribution to the RHS
                rhs_dipole_contrib = np.zeros((dof, dof, nocc, nvir),
                                              dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        rhs_dipole_contrib[x, y] = (
                            0.5 * (np.linalg.multi_dot(
                                [x_minus_y[x].T, dipole_ints_oo[y]]).T +
                                   np.linalg.multi_dot(
                                       [dipole_ints_oo[x], x_minus_y[y]])) -
                            0.5 * (np.linalg.multi_dot(
                                [x_minus_y[x], dipole_ints_vv[y]]) +
                                   np.linalg.multi_dot(
                                       [dipole_ints_vv[x], x_minus_y[y].T]).T))

                rhs_dipole_contrib = rhs_dipole_contrib.reshape(
                    dof**2, nocc, nvir)

                rhs_mo = fock_mo_rhs_1dm + fock_mo_rhs_2dm + rhs_dipole_contrib
                self.ostream.flush()

                # Add DFT E[3] contribution to the RHS:
                if self._dft:
                    gxc_ao = np.zeros((dof**2, nao, nao))

                    for i in range(dof**2):
                        gxc_ao[i] = fock_gxc_ao.alpha_to_numpy(2 * i)

                    gxc_mo = np.array([
                        np.linalg.multi_dot([mo_occ.T, gxc_ao[x], mo_vir])
                        for x in range(dof**2)
                    ])
                    # different factor compared to TDDFT orbital response
                    # because here vectors are scaled by 1/sqrt(2)
                    rhs_mo += 0.5 * gxc_mo

            self.profiler.stop_timer('RHS')

            if self.rank == mpi_master():
                orbrsp_rhs[(w)] = {
                    'cphf_rhs': rhs_mo,
                    'dm_oo': dm_oo,
                    'dm_vv': dm_vv,
                    'x_plus_y_ao': x_plus_y_ao,
                    'x_minus_y_ao': x_minus_y_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'fock_ao_rhs_real': fock_ao_rhs_real,
                    'fock_ao_rhs_imag': fock_ao_rhs_imag,
                    'fock_gxc_ao': fock_gxc_ao,  # None if not DFT
                }
                if (f == 0):
                    tot_rhs_mo = np.concatenate((rhs_mo.real, rhs_mo.imag))
                else:
                    tot_rhs_mo = np.concatenate(
                        (tot_rhs_mo, rhs_mo.real, rhs_mo.imag))

        valstr = '** Time spent on constructing the orbrsp real RHS for '
        valstr += '{} frequencies: '.format(len(self.frequencies))
        valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        if self.rank == mpi_master():
            orbrsp_rhs['cphf_rhs'] = tot_rhs_mo
            return orbrsp_rhs
        else:
            return {}

    def compute_rhs_real(self, molecule, basis, scf_tensors, lr_results):
        """
        Computes the right-hand side (RHS) of the polarizability
        orbital response equation including the necessary density matrices
        using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
        :param lr_results:
            The results from converged linear response calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Run the solver => in parent class

        loop_start_time = tm.time()

        orbrsp_rhs = {}
        for f, w in enumerate(self.frequencies):

            if self.rank == mpi_master():
                self.ostream.print_info(
                    'Building RHS for w = {:4.3f}'.format(w))
                self.ostream.flush()

                # 1) Calculate unrelaxed one-particle and transition density matrix
                ovlp = scf_tensors['S']
                mo = scf_tensors['C']  # only alpha part

                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]

                # TODO: do we keep this factor like that?
                sqrt2 = np.sqrt(2.0)

                # Check if response vectors exist for desired frequency of gradient
                polgrad_sanity_check(self, self.flag, lr_results)

            # TODO: make get_full_solution_vector directly available from the
            # parent class (i.e. lrsolver)?
            full_vec = ([
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ])

            if self.rank == mpi_master():
                # Save the number of vector components
                dof = len(self.vector_components)

                # Extract the excitation and de-excitation components
                # from the full solution vector.
                exc_vec = (1 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                # Construct plus/minus combinations of excitation and
                # de-excitation part
                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # Transform the vectors to the AO basis
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(x_plus_y.shape[0])
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(x_minus_y.shape[0])
                ])

                # Turn them into a list (for AODensityMatrix)
                xpmy_ao_list = list(x_plus_y_ao) + list(x_minus_y_ao)

                # Calculate the symmetrized unrelaxed one-particle density matrix
                # in MO basis
                dm_oo = np.zeros((dof, dof, nocc, nocc))
                dm_vv = np.zeros((dof, dof, nvir, nvir))
                for x in range(dof):
                    for y in range(dof):
                        dm_vv[x, y] = 0.25 * (
                            np.linalg.multi_dot([x_plus_y[x].T, x_plus_y[y]]).T
                            + np.linalg.multi_dot(
                                [x_minus_y[x].T, x_minus_y[y]]).T +
                            np.linalg.multi_dot([x_plus_y[x].T, x_plus_y[y]]) +
                            np.linalg.multi_dot([x_minus_y[x].T, x_minus_y[y]]))
                        dm_oo[x, y] = -0.25 * (
                            np.linalg.multi_dot([x_plus_y[x], x_plus_y[y].T]) +
                            np.linalg.multi_dot([x_minus_y[x], x_minus_y[y].T])
                            + np.linalg.multi_dot([x_plus_y[x], x_plus_y[y].T
                                                  ]).T +
                            np.linalg.multi_dot([x_minus_y[x], x_minus_y[y].T
                                                ]).T)

                # Transform unrelaxed one-particle density matrix to
                # AO basis and create a list
                unrel_dm_ao = np.zeros((dof, dof, nao, nao))
                for x in range(dof):
                    for y in range(dof):
                        unrel_dm_ao[x,
                                    y] = (np.linalg.multi_dot(
                                        [mo_occ, dm_oo[x, y], mo_occ.T]) +
                                          np.linalg.multi_dot(
                                              [mo_vir, dm_vv[x, y], mo_vir.T]))
                dm_ao_list = list(unrel_dm_ao.reshape(dof**2, nao, nao))

                # 2) Construct the right-hand side
                dm_ao_rhs = AODensityMatrix(dm_ao_list + xpmy_ao_list,
                                            denmat.rest)

                if self._dft:
                    # 3) Construct density matrices for E[3] term:
                    # XCIntegrator expects a DM with real and imaginary part,
                    # so we set the imaginary part to zero.
                    # Create lists with the corresponding vector components
                    perturbed_dm_ao_list = []
                    zero_dm_ao_list = []
                    # TODO: only upper triangular matrix and transpose?
                    for x in range(dof):
                        for y in range(dof):
                            perturbed_dm_ao_list.extend([
                                x_minus_y_ao[x], 0 * x_minus_y_ao[x],
                                x_minus_y_ao[y], 0 * x_minus_y_ao[y]
                            ])
                            zero_dm_ao_list.extend(
                                [0 * x_minus_y_ao[x], 0 * x_minus_y_ao[y]])

                    perturbed_dm_ao = AODensityMatrix(perturbed_dm_ao_list,
                                                      denmat.rest)

                    # corresponds to rho^{omega_b,omega_c} in quadratic response,
                    # which is zero for orbital response
                    zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)
            else:
                dof = None
                dm_ao_rhs = AODensityMatrix()
                if self._dft:
                    perturbed_dm_ao = AODensityMatrix()
                    zero_dm_ao = AODensityMatrix()

            dof = self.comm.bcast(dof, root=mpi_master())

            dm_ao_rhs.broadcast(self.rank, self.comm)

            molgrid = dft_dict['molgrid']
            gs_density = dft_dict['gs_density']

            # Fock matrices with corresponding type
            fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
            # Set the vector-related components to general Fock matrix
            # (not 1PDM part)
            for ifock in range(dof**2, dof**2 + 2 * dof):
                fock_ao_rhs.set_fock_type(fockmat.rgenjk, ifock)
            if self._dft:
                perturbed_dm_ao.broadcast(self.rank, self.comm)
                zero_dm_ao.broadcast(self.rank, self.comm)
                # Fock matrix for computing the DFT E[3] term g^xc
                fock_gxc_ao = AOFockMatrix(zero_dm_ao)
                if self.xcfun.is_hybrid():
                    fact_xc = self.xcfun.get_frac_exact_exchange()
                    for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                        fock_ao_rhs.set_scale_factor(fact_xc, ifock)
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.set_scale_factor(fact_xc, ifock)
                        fock_gxc_ao.set_fock_type(fockmat.rgenjkx, ifock)
                    for ifock in range(dof**2):
                        fock_ao_rhs.set_fock_type(fockmat.restjkx, ifock)
                    for ifock in range(dof**2, dof**2 + 2 * dof):
                        fock_ao_rhs.set_fock_type(fockmat.rgenjkx, ifock)
                else:
                    for ifock in range(dof**2):
                        fock_ao_rhs.set_fock_type(fockmat.restj, ifock)
                    for ifock in range(dof**2, dof**2 + 2 * dof):
                        fock_ao_rhs.set_fock_type(fockmat.rgenj, ifock)
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.set_fock_type(fockmat.rgenj, ifock)
            else:
                fock_gxc_ao = None  # None if not DFT

            if self._dft:
                if not self.xcfun.is_hybrid():
                    for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                        fock_gxc_ao.scale(2.0, ifock)
                xc_drv = XCIntegrator(self.comm)
                xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis,
                                          perturbed_dm_ao, zero_dm_ao,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")

                fock_gxc_ao.reduce_sum(self.rank, self.nodes, self.comm)

            self._comp_lr_fock(fock_ao_rhs, dm_ao_rhs, molecule, basis,
                               eri_dict, dft_dict, pe_dict, self.profiler)

            # Calculate the RHS and transform it to the MO basis
            if self.rank == mpi_master():
                # extract the 1PDM contributions
                fock_ao_rhs_1dm = np.zeros((dof**2, nao, nao))
                for i in range(dof**2):
                    fock_ao_rhs_1dm[i] = fock_ao_rhs.alpha_to_numpy(i)

                # Transform to MO basis
                fock_mo_rhs_1dm = np.array([
                    np.linalg.multi_dot([mo_occ.T, fock_ao_rhs_1dm[x], mo_vir])
                    for x in range(dof**2)
                ])

                # extract the x_plus_y and x_minus_y contributions
                # TODO: extract all Fock matrices at the same time?
                fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
                for i in range(dof):
                    fock_ao_rhs_x_plus_y[i] = fock_ao_rhs.alpha_to_numpy(
                        dof**2 + i)
                    fock_ao_rhs_x_minus_y[i] = fock_ao_rhs.alpha_to_numpy(
                        dof**2 + dof + i)

                fock_mo_rhs_2dm = np.zeros((dof, dof, nocc, nvir))
                for x in range(dof):
                    for y in range(dof):
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x], x_plus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x], x_minus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y], x_plus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y], x_minus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x].T, x_plus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x].T, x_minus_y_ao[y]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y].T, x_plus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y].T, x_minus_y_ao[x]])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, tmp, ovlp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x], x_plus_y_ao[y].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x], x_minus_y_ao[y].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y], x_plus_y_ao[x].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y], x_minus_y_ao[x].T]).T
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[x].T, x_plus_y_ao[y].T])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[x].T, x_minus_y_ao[y].T])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_plus_y[y].T, x_plus_y_ao[x].T])
                        fock_mo_rhs_2dm[x, y] += -1.0 * np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                        tmp = np.linalg.multi_dot(
                            [fock_ao_rhs_x_minus_y[y].T, x_minus_y_ao[x].T])
                        fock_mo_rhs_2dm[x, y] += np.linalg.multi_dot(
                            [mo_occ.T, ovlp.T, tmp.T, mo_vir])
                fock_mo_rhs_2dm = 0.25 * fock_mo_rhs_2dm.reshape(
                    dof**2, nocc, nvir)

                # Calculate the dipole contributions to the RHS:
                # Dipole integrals in AO basis
                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_mats = dipole_drv.compute(molecule, basis)
                dipole_ints_ao = np.zeros((dof, nao, nao))
                k = 0
                if 'x' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                    k += 1
                if 'y' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                    k += 1
                if 'z' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.z_to_numpy()

                # Transform them to MO basis (oo and vv blocks only)
                dipole_ints_oo = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
                    for x in range(dof)
                ])
                dipole_ints_vv = np.array([
                    np.linalg.multi_dot([mo_vir.T, dipole_ints_ao[x], mo_vir])
                    for x in range(dof)
                ])

                # Contract with vectors to get dipole contribution to the RHS
                rhs_dipole_contrib = np.zeros((dof, dof, nocc, nvir))
                for x in range(dof):
                    for y in range(dof):
                        rhs_dipole_contrib[x, y] = (
                            0.5 * (np.linalg.multi_dot(
                                [x_minus_y[x].T, dipole_ints_oo[y]]).T +
                                   np.linalg.multi_dot(
                                       [dipole_ints_oo[x], x_minus_y[y]])) -
                            0.5 * (np.linalg.multi_dot(
                                [x_minus_y[x], dipole_ints_vv[y]]) +
                                   np.linalg.multi_dot(
                                       [dipole_ints_vv[x], x_minus_y[y].T]).T))

                rhs_dipole_contrib = rhs_dipole_contrib.reshape(
                    dof**2, nocc, nvir)

                rhs_mo = fock_mo_rhs_1dm + fock_mo_rhs_2dm + rhs_dipole_contrib

                # Add DFT E[3] contribution to the RHS:
                if self._dft:
                    gxc_ao = np.zeros((dof**2, nao, nao))

                    for i in range(dof**2):
                        gxc_ao[i] = fock_gxc_ao.alpha_to_numpy(2 * i)

                    gxc_mo = np.array([
                        np.linalg.multi_dot([mo_occ.T, gxc_ao[x], mo_vir])
                        for x in range(dof**2)
                    ])
                    # different factor compared to TDDFT orbital response
                    # because here vectors are scaled by 1/sqrt(2)
                    rhs_mo += 0.5 * gxc_mo

            self.profiler.stop_timer('RHS')

            if self.rank == mpi_master():
                orbrsp_rhs[(w)] = {
                    'cphf_rhs': rhs_mo,
                    'dm_oo': dm_oo,
                    'dm_vv': dm_vv,
                    'x_plus_y_ao': x_plus_y_ao,
                    'x_minus_y_ao': x_minus_y_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'fock_ao_rhs': fock_ao_rhs,
                    'fock_gxc_ao': fock_gxc_ao,  # None if not DFT
                }
                if (f == 0):
                    tot_rhs_mo = rhs_mo
                else:
                    tot_rhs_mo = np.append(tot_rhs_mo, rhs_mo, axis=0)

        valstr = '** Time spent on constructing the orbrsp real RHS for '
        valstr += '{} frequencies: '.format(len(self.frequencies))
        valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        if self.rank == mpi_master():
            orbrsp_rhs['cphf_rhs'] = tot_rhs_mo
            return orbrsp_rhs
        else:
            return {}

    # NOTES:
    #   - epsilon_dm_ao not returned from cphfsolver,
    #     to be calculated inside compute_omega
    #   - fock_ao_rhs and fock_gxc_ao come from cphfsolver dictionary
    #   - fock_lambda not returned yet, put in dictionary from cphfsolver
    #     (otherwise needs to be recalculated)

    def compute_omega(self, molecule, basis, scf_tensors, lr_results):
        """
        Guides the calculation of the polarizability Lagrange multipliers
        for the overlap matrix according to is_complex instance variable
        """

        if self.is_complex:
            return self.compute_omega_complex(molecule, basis, scf_tensors,
                                              lr_results)
        else:
            return self.compute_omega_real(molecule, basis, scf_tensors,
                                           lr_results)

    def compute_omega_real(self, molecule, basis, scf_tensors, lr_results):
        """
        Calculates the real  polarizability Lagrange multipliers for the
        overlap matrix.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results from the linear response calculation.
        """

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        n_freqs = len(self.frequencies)

        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):
            self.ostream.print_info('Building omega for w = {:4.3f}'.format(w))
            self.ostream.flush()

            if self.rank == mpi_master():

                # Get overlap, MO coefficients from scf_tensors
                ovlp = scf_tensors['S']
                nocc = molecule.number_of_alpha_electrons()
                mo = scf_tensors['C']
                mo_occ = mo[:, :nocc]
                mo_vir = mo[:, nocc:]
                nocc = mo_occ.shape[1]
                nvir = mo_vir.shape[1]
                nao = mo_occ.shape[0]

                mo_energies = scf_tensors['E']
                eocc = mo_energies[:nocc]
                evir = mo_energies[nocc:]
                eo_diag = np.diag(eocc)
                ev_diag = np.diag(evir)

                # Get fock matrices from cphf_results
                fock_ao_rhs = self.cphf_results[w]['fock_ao_rhs']
                fock_gxc_ao = self.cphf_results[w]['fock_gxc_ao']
                dm_oo = self.cphf_results[w]['dm_oo']
                dm_vv = self.cphf_results[w]['dm_vv']

                # TODO: MPI should this be done before loop over freqs?
                all_cphf_ov = self.cphf_results['cphf_ov']
                dof = int(all_cphf_ov.shape[0] / n_freqs)
                cphf_ov = all_cphf_ov.reshape(n_freqs, dof, nocc, nvir)[f]

                # TODO: do we keep this factor like that?
                sqrt2 = np.sqrt(2.0)

            # TODO: make get_full_solution_vector directly available from the
            # parent class (i.e. lrsolver)?
            full_vec = ([
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ])

            if self.rank == mpi_master():
                # Save the number of vector components
                dof = len(self.vector_components)

                # Extract the excitation and de-excitation components
                # from the full solution vector.
                exc_vec = (1 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # Get dipole moment integrals and transform to MO
                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_mats = dipole_drv.compute(molecule, basis)
                dipole_ints_ao = np.zeros((dof, nao, nao))
                k = 0
                if 'x' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                    k += 1
                if 'y' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                    k += 1
                if 'z' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.z_to_numpy()

                # Transform them to MO basis (oo and ov blocks only)
                dipole_ints_oo = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
                    for x in range(dof)
                ])
                dipole_ints_ov = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_vir])
                    for x in range(dof)
                ])

                # Calculate the dipole moment integrals' contribution to omega
                dipole_ints_contrib_ao = np.zeros((dof, dof, nao, nao))
                for x in range(dof):
                    for y in range(dof):
                        tmp_oo = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x], dipole_ints_ov[y].T
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_ov[x], x_minus_y[y].T]))
                        tmp_ov = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x].T, dipole_ints_oo[y]
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_oo[x].T, x_minus_y[y]]))
                        tmp_vv = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x].T, dipole_ints_ov[y]
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_ov[x].T, x_minus_y[y]]))
                        dipole_ints_contrib_ao[x, y] = (
                            np.linalg.multi_dot([mo_occ, tmp_oo, mo_occ.T]) +
                            np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]) +
                            np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]).T +
                            np.linalg.multi_dot([mo_vir, tmp_vv, mo_vir.T]))

                # Transform the vectors to the AO basis
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(dof)
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(dof)
                ])

                # The density matrix; only alpha block;
                # Only works for the restricted case
                # (since scf_tensors['C'] only gives alpha block...)
                D_occ = np.matmul(mo_occ, mo_occ.T)
                D_vir = np.matmul(mo_vir, mo_vir.T)

                # Construct fock_lambda (or fock_cphf)
                cphf_ao = np.array([
                    np.linalg.multi_dot([mo_occ, cphf_ov[x], mo_vir.T])
                    for x in range(dof**2)
                ])
                cphf_ao_list = list([cphf_ao[x] for x in range(dof**2)])
                ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)
            else:
                dof = None
                ao_density_cphf = AODensityMatrix()

            dof = self.comm.bcast(dof, root=mpi_master())
            ao_density_cphf.broadcast(self.rank, self.comm)
            fock_cphf = AOFockMatrix(ao_density_cphf)

            # TODO: what has to be on MPI master and what not?
            self._comp_lr_fock(fock_cphf, ao_density_cphf, molecule, basis,
                               eri_dict, dft_dict, pe_dict, self.profiler)

            # For now we:
            # - loop over indices m and n
            # - select component m or n in x_plus_y, x_minus_y,
            # fock_ao_rhs and fock_lambda
            # - symmetrize with respect to m and n (only 2PDM?)
            # Notes: fock_ao_rhs is a list with dof**2 matrices corresponding to
            # the contraction of the 1PDMs with ERIs; dof matrices corresponding
            # to the contraction of x_plus_y; and other dof matrices corresponding to
            # the contraction of x_minus_y.

            if self.rank == mpi_master():
                # TODO: what shape should we use: (dof**2, nao, nao)
                #       or (dof, dof, nao, nao)?
                omega = np.zeros((dof * dof, nao, nao))

                # Calculate omega (without for-loops, only the diagonal parts
                # possible for now)
                # Construct epsilon_dm_ao
                epsilon_dm_ao = np.zeros((dof, dof, nao, nao))
                epsilon_cphf_ao = np.zeros((dof, dof, nao, nao))
                for x in range(dof):
                    for y in range(dof):
                        epsilon_dm_ao[x, y] = -1.0 * np.linalg.multi_dot(
                            [mo_occ, eo_diag, dm_oo[x, y], mo_occ.T])
                        epsilon_dm_ao[x, y] -= np.linalg.multi_dot(
                            [mo_vir, ev_diag, dm_vv[x, y], mo_vir.T])
                        epsilon_cphf_ao[x, y] = np.linalg.multi_dot([
                            mo_occ, eo_diag,
                            cphf_ov.reshape(dof, dof, nocc, nvir)[x, y],
                            mo_vir.T
                        ])

                # OV + VO
                epsilon_dm_ao -= (epsilon_cphf_ao +
                                  epsilon_cphf_ao.transpose(0, 1, 3, 2))

                for m in range(dof):
                    for n in range(dof):
                        # TODO: move outside for-loop when all Fock matrices can be
                        # extracted into a numpy array at the same time.

                        # Because the excitation vector is not symmetric,
                        # we need both the matrix (OO block in omega, and probably VO)
                        # and its transpose (VV, OV blocks)
                        # this comes from the transformation of the 2PDM contribution
                        # from MO to AO basis
                        fock_ao_rhs_1_m = fock_ao_rhs.alpha_to_numpy(
                            dof**2 + m)  # x_plus_y
                        fock_ao_rhs_2_m = fock_ao_rhs.alpha_to_numpy(
                            dof**2 + dof + m)  # x_minus_y

                        fock_ao_rhs_1_n = fock_ao_rhs.alpha_to_numpy(
                            dof**2 + n)  # x_plus_y
                        fock_ao_rhs_2_n = fock_ao_rhs.alpha_to_numpy(
                            dof**2 + dof + n)  # x_minus_y

                        Fp1_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m.T, x_plus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n.T, x_plus_y_ao[m], ovlp.T]))

                        Fm1_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m.T, x_minus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n.T, x_minus_y_ao[m], ovlp.T]))

                        Fp2_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m, x_plus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n, x_plus_y_ao[m], ovlp.T]))

                        Fm2_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m, x_minus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n, x_minus_y_ao[m], ovlp.T]))

                        Fp1_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m, x_plus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n, x_plus_y_ao[m].T, ovlp.T]))

                        Fm1_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m, x_minus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n, x_minus_y_ao[m].T, ovlp.T]))

                        Fp2_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m.T, x_plus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n.T, x_plus_y_ao[m].T, ovlp.T]))

                        Fm2_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m.T, x_minus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n.T, x_minus_y_ao[m].T, ovlp.T]))
                        # We see that:
                        # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
                        # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

                        # Compute the contributions from the 2PDM and the relaxed 1PDM
                        # to the omega Lagrange multipliers:
                        fmat = (fock_cphf.alpha_to_numpy(m * dof + n) +
                                fock_cphf.alpha_to_numpy(m * dof + n).T +
                                fock_ao_rhs.alpha_to_numpy(m * dof + n))
                        # dof=3  (0,0), (0,1), (0,2); (1,0), (1,1), (1,2),
                        #        (2,0), (2,1), (2,2) * dof
                        # gamma_{zx} =

                        omega_1pdm_2pdm_contribs = -(np.linalg.multi_dot([
                            D_vir, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]) + np.linalg.multi_dot([
                            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]) + np.linalg.multi_dot([
                            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]).T + np.linalg.multi_dot([
                            D_occ, Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ
                        ]) + np.linalg.multi_dot([D_occ, fmat, D_occ]))

                        omega[m * dof + n] = (epsilon_dm_ao[m, n] +
                                              omega_1pdm_2pdm_contribs +
                                              dipole_ints_contrib_ao[m, n])

                        if self._dft:
                            factor = -0.5
                            omega[m * dof + n] += factor * np.linalg.multi_dot([
                                D_occ,
                                fock_gxc_ao.alpha_to_numpy(2 *
                                                           (m * dof + n)), D_occ
                            ])

                # add omega multipliers in AO basis to cphf_results dictionary
                self.cphf_results[(w)]['omega_ao'] = omega
                self.cphf_results[(w)]['lambda_mo'] = cphf_ov

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing omega multipliers '
            valstr += 'for {} frequencies: '.format(n_freqs)
            valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_omega_complex(self, molecule, basis, scf_tensors, lr_results):
        """
        Calculates the complex polarizability Lagrange multipliers for the
        overlap matrix.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results from the CPP calculation.
        """

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        n_freqs = len(self.frequencies)

        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):
            self.ostream.print_info('Building omega for w = {:4.3f}'.format(w))
            self.ostream.flush()

            if self.rank == mpi_master():

                # Get overlap, MO coefficients from scf_tensors
                ovlp = scf_tensors['S']
                nocc = molecule.number_of_alpha_electrons()
                mo = scf_tensors['C']
                mo_occ = mo[:, :nocc]
                mo_vir = mo[:, nocc:]
                nocc = mo_occ.shape[1]
                nvir = mo_vir.shape[1]
                nao = mo_occ.shape[0]

                mo_energies = scf_tensors['E']
                eocc = mo_energies[:nocc]
                evir = mo_energies[nocc:]
                eo_diag = np.diag(eocc)
                ev_diag = np.diag(evir)

                # Get fock matrices from cphf_results
                fock_ao_rhs_real = self.cphf_results[w]['fock_ao_rhs_real']
                fock_ao_rhs_imag = self.cphf_results[w]['fock_ao_rhs_imag']
                fock_gxc_ao = self.cphf_results[w]['fock_gxc_ao']
                dm_oo = self.cphf_results[w]['dm_oo']  # complex
                dm_vv = self.cphf_results[w]['dm_vv']  # complex

                # TODO: MPI should this be done before loop over freqs?
                all_cphf_ov = self.cphf_results[
                    'cphf_ov']  # FIXME will be complex ?
                n_lambdas = int(all_cphf_ov.shape[0] / n_freqs)
                cphf_ov = all_cphf_ov.reshape(n_freqs, n_lambdas, nocc, nvir)[f]

                # TODO: do we keep this factor like that?
                sqrt2 = np.sqrt(2.0)

            # TODO: make get_full_solution_vector directly available from the
            # parent class (i.e. lrsolver)?
            full_vec = ([
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ])  # complex

            if self.rank == mpi_master():
                # Save the number of vector components
                dof = len(self.vector_components)

                # reshape lambda vector to complex
                cphf_ov = cphf_ov[:dof**2] + 1j * cphf_ov[dof**2:]

                # Extract the excitation and de-excitation components
                # from the full solution vector.
                exc_vec = (1 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # Get dipole moment integrals and transform to MO
                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_mats = dipole_drv.compute(molecule, basis)
                dipole_ints_ao = np.zeros((dof, nao, nao))
                k = 0
                if 'x' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                    k += 1
                if 'y' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                    k += 1
                if 'z' in self.vector_components:
                    dipole_ints_ao[k] = dipole_mats.z_to_numpy()

                # Transform them to MO basis (oo and ov blocks only)
                dipole_ints_oo = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
                    for x in range(dof)
                ])
                dipole_ints_ov = np.array([
                    np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_vir])
                    for x in range(dof)
                ])

                # Calculate the dipole moment integrals' contribution to omega
                dipole_ints_contrib_ao = np.zeros((dof, dof, nao, nao),
                                                  dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        tmp_oo = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x], dipole_ints_ov[y].T
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_ov[x], x_minus_y[y].T]))
                        tmp_ov = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x].T, dipole_ints_oo[y]
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_oo[x].T, x_minus_y[y]]))
                        tmp_vv = 0.5 * (np.linalg.multi_dot([
                            x_minus_y[x].T, dipole_ints_ov[y]
                        ]).T + np.linalg.multi_dot(
                            [dipole_ints_ov[x].T, x_minus_y[y]]))
                        dipole_ints_contrib_ao[x, y] = (
                            np.linalg.multi_dot([mo_occ, tmp_oo, mo_occ.T]) +
                            np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]) +
                            np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]).T +
                            np.linalg.multi_dot([mo_vir, tmp_vv, mo_vir.T]))

                # Transform the vectors to the AO basis
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(dof)
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(dof)
                ])

                # The density matrix; only alpha block;
                # Only works for the restricted case
                # (since scf_tensors['C'] only gives alpha block...)
                D_occ = np.matmul(mo_occ, mo_occ.T)
                D_vir = np.matmul(mo_vir, mo_vir.T)

                # Construct fock_lambda (or fock_cphf)
                cphf_ao = np.array([
                    np.linalg.multi_dot([mo_occ, cphf_ov[x], mo_vir.T])
                    for x in range(dof**2)
                ])
                cphf_ao_list_real = list(
                    np.array([cphf_ao[x].real for x in range(dof**2)]))
                cphf_ao_list_imag = list(
                    np.array([cphf_ao[x].imag for x in range(dof**2)]))
                ao_density_cphf_real = AODensityMatrix(cphf_ao_list_real,
                                                       denmat.rest)
                ao_density_cphf_imag = AODensityMatrix(cphf_ao_list_imag,
                                                       denmat.rest)
            else:
                dof = None
                ao_density_cphf_real = AODensityMatrix()
                ao_density_cphf_imag = AODensityMatrix()

            dof = self.comm.bcast(dof, root=mpi_master())
            ao_density_cphf_real.broadcast(self.rank, self.comm)
            ao_density_cphf_imag.broadcast(self.rank, self.comm)
            fock_cphf_real = AOFockMatrix(ao_density_cphf_real)
            fock_cphf_imag = AOFockMatrix(ao_density_cphf_imag)

            # TODO: what has to be on MPI master and what not?
            self._comp_lr_fock(fock_cphf_real, ao_density_cphf_real, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)
            self._comp_lr_fock(fock_cphf_imag, ao_density_cphf_imag, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)

            # For now we:
            # - loop over indices m and n
            # - select component m or n in x_plus_y, x_minus_y,
            # fock_ao_rhs and fock_lambda
            # - symmetrize with respect to m and n (only 2PDM?)
            # Notes: fock_ao_rhs is a list with dof**2 matrices corresponding to
            # the contraction of the 1PDMs with ERIs; dof matrices corresponding
            # to the contraction of x_plus_y; and other dof matrices corresponding to
            # the contraction of x_minus_y.

            if self.rank == mpi_master():
                # TODO: what shape should we use: (dof**2, nao, nao)
                #       or (dof, dof, nao, nao)?
                omega = np.zeros((dof * dof, nao, nao), dtype=np.complex_)

                # Calculate omega (without for-loops, only the diagonal parts
                # possible for now)
                # Construct epsilon_dm_ao
                epsilon_dm_ao = np.zeros((dof, dof, nao, nao),
                                         dtype=np.complex_)
                epsilon_cphf_ao = np.zeros((dof, dof, nao, nao),
                                           dtype=np.complex_)
                for x in range(dof):
                    for y in range(dof):
                        epsilon_dm_ao[x, y] = -1.0 * np.linalg.multi_dot(
                            [mo_occ, eo_diag, dm_oo[x, y], mo_occ.T])
                        epsilon_dm_ao[x, y] -= np.linalg.multi_dot(
                            [mo_vir, ev_diag, dm_vv[x, y], mo_vir.T])
                        epsilon_cphf_ao[x, y] = np.linalg.multi_dot([
                            mo_occ, eo_diag,
                            cphf_ov.reshape(dof, dof, nocc, nvir)[x, y],
                            mo_vir.T
                        ])

                # OV + VO
                epsilon_dm_ao -= (epsilon_cphf_ao +
                                  epsilon_cphf_ao.transpose(0, 1, 3, 2))

                for m in range(dof):
                    for n in range(dof):
                        # TODO: move outside for-loop when all Fock matrices can be
                        # extracted into a numpy array at the same time.

                        # Because the excitation vector is not symmetric,
                        # we need both the matrix (OO block in omega, and probably VO)
                        # and its transpose (VV, OV blocks)
                        # this comes from the transformation of the 2PDM contribution
                        # from MO to AO basis
                        # complex
                        fock_ao_rhs_1_m = (
                            fock_ao_rhs_real.alpha_to_numpy(dof**2 + m) +
                            1j * fock_ao_rhs_imag.alpha_to_numpy(dof**2 + m)
                        )  # x_plus_y
                        fock_ao_rhs_2_m = (
                            fock_ao_rhs_real.alpha_to_numpy(dof**2 + dof + m) +
                            1j *
                            fock_ao_rhs_imag.alpha_to_numpy(dof**2 + dof + m)
                        )  # x_minus_y

                        fock_ao_rhs_1_n = (
                            fock_ao_rhs_real.alpha_to_numpy(dof**2 + n) +
                            1j * fock_ao_rhs_imag.alpha_to_numpy(dof**2 + n)
                        )  # x_plus_y
                        fock_ao_rhs_2_n = (
                            fock_ao_rhs_real.alpha_to_numpy(dof**2 + dof + n) +
                            1j *
                            fock_ao_rhs_imag.alpha_to_numpy(dof**2 + dof + n)
                        )  # x_minus_y

                        Fp1_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m.T, x_plus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n.T, x_plus_y_ao[m], ovlp.T]))

                        Fm1_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m.T, x_minus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n.T, x_minus_y_ao[m], ovlp.T]))

                        Fp2_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m, x_plus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n, x_plus_y_ao[m], ovlp.T]))

                        Fm2_vv = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m, x_minus_y_ao[n], ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n, x_minus_y_ao[m], ovlp.T]))

                        Fp1_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m, x_plus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n, x_plus_y_ao[m].T, ovlp.T]))

                        Fm1_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m, x_minus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n, x_minus_y_ao[m].T, ovlp.T]))

                        Fp2_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_1_m.T, x_plus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_1_n.T, x_plus_y_ao[m].T, ovlp.T]))

                        Fm2_oo = 0.25 * (np.linalg.multi_dot([
                            fock_ao_rhs_2_m.T, x_minus_y_ao[n].T, ovlp.T
                        ]) + np.linalg.multi_dot(
                            [fock_ao_rhs_2_n.T, x_minus_y_ao[m].T, ovlp.T]))
                        # We see that:
                        # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
                        # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

                        # Compute the contributions from the 2PDM and the relaxed 1PDM
                        # to the omega Lagrange multipliers:
                        fmat = ((fock_cphf_real.alpha_to_numpy(m * dof + n) +
                                 fock_cphf_real.alpha_to_numpy(m * dof + n).T +
                                 fock_ao_rhs_real.alpha_to_numpy(m * dof + n)) +
                                1j *
                                (fock_cphf_imag.alpha_to_numpy(m * dof + n) +
                                 fock_cphf_imag.alpha_to_numpy(m * dof + n).T +
                                 fock_ao_rhs_imag.alpha_to_numpy(m * dof + n)))

                        # dof=3  (0,0), (0,1), (0,2); (1,0), (1,1), (1,2),
                        #        (2,0), (2,1), (2,2) * dof
                        # gamma_{zx} =

                        omega_1pdm_2pdm_contribs = -(np.linalg.multi_dot([
                            D_vir, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]) + np.linalg.multi_dot([
                            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]) + np.linalg.multi_dot([
                            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
                        ]).T + np.linalg.multi_dot([
                            D_occ, Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ
                        ]) + np.linalg.multi_dot([D_occ, fmat, D_occ]))

                        omega[m * dof + n] = (epsilon_dm_ao[m, n] +
                                              omega_1pdm_2pdm_contribs +
                                              dipole_ints_contrib_ao[m, n])

                        if self._dft:
                            factor = -0.5
                            omega[m * dof + n] += factor * np.linalg.multi_dot([
                                D_occ,
                                fock_gxc_ao.alpha_to_numpy(2 *
                                                           (m * dof + n)), D_occ
                            ])

                # add omega multipliers in AO basis to cphf_results dictionary
                self.cphf_results[(w)]['omega_ao'] = omega
                self.cphf_results[(w)]['lambda_mo'] = cphf_ov

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing omega multipliers '
            valstr += 'for {} frequencies: '.format(n_freqs)
            valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def print_cphf_header(self, title):
        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 70

        # print general info
        cur_str = 'Solver Type                     : '
        if self.use_subspace_solver:
            cur_str += 'Iterative Subspace Algorithm'
        else:
            cur_str += 'Conjugate Gradient'
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = ('Number of frequencies           : ' +
                   str(len(self.frequencies)))
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Vector components               : ' + self.vector_components
        self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

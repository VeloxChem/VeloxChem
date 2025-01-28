#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import time as tm

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import XCIntegrator
from .veloxchemlib import mpi_master, denmat
from .cphfsolver import CphfSolver
from .inputparser import parse_input
from .sanitychecks import polgrad_sanity_check
from .oneeints import compute_electric_dipole_integrals
from .distributedarray import DistributedArray
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .profiler import Profiler


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

        #self.sqrt2 = np.sqrt(2.0)

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

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict, lr_results):

        if self.is_complex:
            return self.compute_rhs_complex(molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict,
                                            lr_results)
        else:
            return self.compute_rhs_real(molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict,
                                         lr_results)

    def compute_rhs_complex(self, molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict, lr_results):
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

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Run the solver => in parent class

        loop_start_time = tm.time()

        nocc = molecule.number_of_alpha_electrons()
        nao = basis.get_dimensions_of_basis()

        # number of vector components
        dof = len(self.vector_components)

        if self.rank == mpi_master():
            # check if response vectors exist for desired frequency of gradient
            polgrad_sanity_check(self, self.flag, lr_results)

            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']  # only alpha part

            # reduce dimensions of RHS to unique operator component combinations
            xy_pairs = [(x,y) for x in range(3) for y in range(x,3)]
            dof_red = len(xy_pairs)
        else:
            density = None
            mo = None
            dof_red = None

        mo = self.comm.bcast(mo, root=mpi_master())
        density = self.comm.bcast(density, root=mpi_master())
        dof_red = self.comm.bcast(dof_red, root=mpi_master())

        # TODO: double check dft_dict['gs_density']
        # TODO: maybe rename to mol_grid
        molgrid = dft_dict['molgrid']
        gs_density = [density]

        # MO coefficients
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        orbrsp_rhs = {}

        dist_cphf_rhs = []

        for f, w in enumerate(self.frequencies):

            if self.rank == mpi_master():
                self.ostream.print_info(
                    'Building RHS for w = {:4.3f}'.format(w))
                self.ostream.flush()

            full_vec = [
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ]

            if self.rank == mpi_master():

                # Note: polorbitalresponse uses r instead of mu for dipole operator
                for idx in range(len(full_vec)):
                    full_vec[idx] *= -1.0

                # number of vector components
                dof = len(self.vector_components)

                # extract the excitation and de-excitation components
                # from the full solution vector.
                sqrt2 = np.sqrt(2.0)
                exc_vec = (1.0 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1.0 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                # construct plus/minus combinations of excitation and
                # de-excitation part
                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # transform to AO basis: mi,xia,na->xmn
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(x_plus_y.shape[0])
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(x_minus_y.shape[0])
                ])

                # turn them into a list for AODensityMatrix
                xpmy_ao_list_real = list(np.array(x_plus_y_ao.real)) + list(
                    np.array(x_minus_y_ao.real))
                xpmy_ao_list_imag = list(np.array(x_plus_y_ao.imag)) + list(
                    np.array(x_minus_y_ao.imag))

                # calculate symmetrized unrelaxed one-particle density matrix
                unrel_dm_ao, dm_oo, dm_vv = self.calculate_unrel_dm(molecule, scf_tensors,
                                                      x_plus_y, x_minus_y)
                # create lists
                dm_ao_list_real = list(
                    np.array(unrel_dm_ao.real).reshape(dof**2, nao, nao))
                dm_ao_list_imag = list(
                    np.array(unrel_dm_ao.imag).reshape(dof**2, nao, nao))

                # density matrix for RHS
                dm_ao_rhs_real_list = dm_ao_list_real + xpmy_ao_list_real
                dm_ao_rhs_imag_list = dm_ao_list_imag + xpmy_ao_list_imag

                if self._dft:
                    # construct density matrices for E[3] term
                    #= self.construct_dft_e3_dm_complex(x_minus_y_ao)

                    perturbed_dm_ao_list_rere = []
                    perturbed_dm_ao_list_imim = []
                    perturbed_dm_ao_list_reim = []
                    perturbed_dm_ao_list_imre = []
                    zero_dm_ao_list = []

                    for x in range(dof):
                        for y in range(dof):
                            perturbed_dm_ao_list_rere.extend([
                                np.array(x_minus_y_ao[x].real),
                                np.array(0 * x_minus_y_ao[x].real),
                                np.array(x_minus_y_ao[y].real),
                                np.array(0 * x_minus_y_ao[y].real)])
                            perturbed_dm_ao_list_imim.extend([
                                np.array(x_minus_y_ao[x].imag),
                                np.array(0 * x_minus_y_ao[x].imag),
                                np.array(x_minus_y_ao[y].imag),
                                np.array(0 * x_minus_y_ao[y].imag)])

                            # complex cross-terms
                            perturbed_dm_ao_list_reim.extend([
                                np.array(x_minus_y_ao[x].real),
                                np.array(0 * x_minus_y_ao[x].real),
                                np.array(x_minus_y_ao[y].imag),
                                np.array(0 * x_minus_y_ao[y].imag)])
                            perturbed_dm_ao_list_imre.extend([
                                np.array(x_minus_y_ao[x].imag),
                                np.array(0 * x_minus_y_ao[x].imag),
                                np.array(x_minus_y_ao[y].real),
                                np.array(0 * x_minus_y_ao[y].real)])

                            zero_dm_ao_list.extend([
                                np.array(0 * x_minus_y_ao[x].real),
                                np.array(0 * x_minus_y_ao[y].real)])

            else:
                dof = None
                dm_ao_rhs_real_list = None
                dm_ao_rhs_imag_list = None

                if self._dft:
                    perturbed_dm_ao_list_rere = None
                    perturbed_dm_ao_list_imim = None
                    perturbed_dm_ao_list_reim = None
                    perturbed_dm_ao_list_imre = None
                    zero_dm_ao_list = None

            dof = self.comm.bcast(dof, root=mpi_master())

            dm_ao_rhs_real_list = self.comm.bcast(dm_ao_rhs_real_list, root=mpi_master())
            dm_ao_rhs_imag_list = self.comm.bcast(dm_ao_rhs_imag_list, root=mpi_master())

            if self._dft:
                perturbed_dm_ao_list_rere = self.comm.bcast(perturbed_dm_ao_list_rere, root=mpi_master())
                perturbed_dm_ao_list_imim = self.comm.bcast(perturbed_dm_ao_list_imim, root=mpi_master())
                perturbed_dm_ao_list_reim = self.comm.bcast(perturbed_dm_ao_list_reim, root=mpi_master())
                perturbed_dm_ao_list_imre = self.comm.bcast(perturbed_dm_ao_list_imre, root=mpi_master())
                zero_dm_ao_list = self.comm.bcast(zero_dm_ao_list, root=mpi_master())

            molgrid = dft_dict['molgrid']
            gs_density = dft_dict['gs_density']

            # Fock matrices with corresponding type
            # set the vector-related components to general Fock matrix
            # (not 1PDM part)

            if self._dft:
                # TODO: use deepcopy
                fock_gxc_ao_rere = []
                fock_gxc_ao_imim = []
                fock_gxc_ao_reim = []
                fock_gxc_ao_imre = []
                for i_mat in range(len(zero_dm_ao_list)):
                    fock_gxc_ao_rere.append(zero_dm_ao_list[i_mat].copy())
                    fock_gxc_ao_imim.append(zero_dm_ao_list[i_mat].copy())
                    fock_gxc_ao_reim.append(zero_dm_ao_list[i_mat].copy())
                    fock_gxc_ao_imre.append(zero_dm_ao_list[i_mat].copy())
            else:
                fock_gxc_ao_rere = None
                fock_gxc_ao_imim = None
                fock_gxc_ao_reim = None
                fock_gxc_ao_imre = None

            if self._dft:
                xc_drv = XCIntegrator()

                xc_drv.integrate_kxc_fock(fock_gxc_ao_rere, molecule, basis,
                                          perturbed_dm_ao_list_rere, zero_dm_ao_list,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")
                xc_drv.integrate_kxc_fock(fock_gxc_ao_imim, molecule, basis,
                                          perturbed_dm_ao_list_imim, zero_dm_ao_list,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")
                xc_drv.integrate_kxc_fock(fock_gxc_ao_reim, molecule, basis,
                                          perturbed_dm_ao_list_reim, zero_dm_ao_list,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")
                xc_drv.integrate_kxc_fock(fock_gxc_ao_imre, molecule, basis,
                                          perturbed_dm_ao_list_imre, zero_dm_ao_list,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")

                for idx in range(len(fock_gxc_ao_rere)):
                    fock_gxc_ao_rere[idx] = self.comm.reduce(fock_gxc_ao_rere[idx], root=mpi_master())
                for idx in range(len(fock_gxc_ao_imim)):
                    fock_gxc_ao_imim[idx] = self.comm.reduce(fock_gxc_ao_imim[idx], root=mpi_master())
                for idx in range(len(fock_gxc_ao_reim)):
                    fock_gxc_ao_reim[idx] = self.comm.reduce(fock_gxc_ao_reim[idx], root=mpi_master())
                for idx in range(len(fock_gxc_ao_imre)):
                    fock_gxc_ao_imre[idx] = self.comm.reduce(fock_gxc_ao_imre[idx], root=mpi_master())

            fock_ao_rhs_real = self._comp_lr_fock(dm_ao_rhs_real_list, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)
            fock_ao_rhs_imag = self._comp_lr_fock(dm_ao_rhs_imag_list, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)

            # calculate the RHS
            if self.rank == mpi_master():
                # extract the 1PDM contributions
                fock_ao_rhs_1pdm = np.zeros((dof**2, nao, nao), dtype=np.dtype('complex128'))
                fock_ao_rhs_1pdm_real = np.zeros((dof**2, nao, nao))
                fock_ao_rhs_1pdm_imag = np.zeros((dof**2, nao, nao))
                for i in range(dof**2):
                    fock_ao_rhs_1pdm_real[i] = fock_ao_rhs_real[i]
                    fock_ao_rhs_1pdm_imag[i] = fock_ao_rhs_imag[i]
                # combine to complex array
                fock_ao_rhs_1pdm = fock_ao_rhs_1pdm_real + 1j * fock_ao_rhs_1pdm_imag

                # transform to MO basis: mi,xmn,na->xia
                fock_mo_rhs_1pdm = np.array([
                    np.linalg.multi_dot([mo_occ.T, fock_ao_rhs_1pdm[x], mo_vir])
                    for x in range(dof**2)
                ])

                # extract the x_plus_y and x_minus_y contributions
                fock_ao_rhs_x_plus_y_real = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y_real = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_plus_y_imag = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y_imag = np.zeros((dof, nao, nao))
                for i in range(dof):
                    fock_ao_rhs_x_plus_y_real[i] = fock_ao_rhs_real[ dof**2 + i]
                    fock_ao_rhs_x_minus_y_real[i] = fock_ao_rhs_real[ dof**2 + dof + i]
                    fock_ao_rhs_x_plus_y_imag[i] = fock_ao_rhs_imag[ dof**2 + i]
                    fock_ao_rhs_x_minus_y_imag[i] = fock_ao_rhs_imag[ dof**2 + dof + i]
                # combine to complex
                fock_ao_rhs_x_plus_y = (fock_ao_rhs_x_plus_y_real +
                                        1j * fock_ao_rhs_x_plus_y_imag)
                fock_ao_rhs_x_minus_y = (fock_ao_rhs_x_minus_y_real +
                                         1j * fock_ao_rhs_x_minus_y_imag)
               
                # calculate 2-particle density matrix contribution
                fock_mo_rhs_2pdm = self.calculate_rhs_2pdm_contrib(molecule, scf_tensors,
                    x_plus_y_ao, x_minus_y_ao, fock_ao_rhs_x_plus_y, fock_ao_rhs_x_minus_y)

                # Calculate dipole contribution
                rhs_dipole_contrib = self.calculate_rhs_dipole_contrib(
                    molecule, basis, scf_tensors, x_minus_y)

                # sum RHS contributions
                rhs_mo = fock_mo_rhs_1pdm + fock_mo_rhs_2pdm + rhs_dipole_contrib

                # add DFT E[3] contribution to the RHS
                if self._dft:
                    gxc_ao = np.zeros((dof**2, nao, nao), dtype=np.dtype('complex128'))

                    for i in range(dof**2):
                        gxc_ao[i] = (fock_gxc_ao_rere[2 * i]
                                     - fock_gxc_ao_imim[2 * i]
                                     + 1j * (fock_gxc_ao_reim[2 * i]
                                     + fock_gxc_ao_imre[2 * i]
                                     ))

                    # transform to MO basis: mi,xmn,na->xia
                    gxc_mo = np.array([
                        np.linalg.multi_dot([mo_occ.T, gxc_ao[x], mo_vir])
                        for x in range(dof**2)
                    ])
                    # different factor compared to TDDFT orbital response
                    # because here vectors are scaled by 1/sqrt(2)
                    rhs_mo += 0.5 * (gxc_mo)

            self.profiler.stop_timer('RHS')

            if self.rank == mpi_master():
                # reduce dimensions of RHS to unique operator component combinations
                rhs_tmp = rhs_mo.reshape(dof, dof, nocc, nvir).copy()
                rhs_red = []

                for x, y in xy_pairs:
                    rhs_red.append(rhs_tmp[x, y].copy())
                rhs_red = np.array(rhs_red)

                # TODO purge the dict
                orbrsp_rhs[(w)] = {
                    'cphf_rhs': rhs_red,
                    'dm_oo': dm_oo,
                    'dm_vv': dm_vv,
                    'x_plus_y_ao': x_plus_y_ao,
                    'x_minus_y_ao': x_minus_y_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'fock_ao_rhs_real': fock_ao_rhs_real,
                    'fock_ao_rhs_imag': fock_ao_rhs_imag,
                    'fock_gxc_ao_rere': fock_gxc_ao_rere,  # None if not DFT
                    'fock_gxc_ao_imim': fock_gxc_ao_imim,  # None if not DFT
                    'fock_gxc_ao_reim': fock_gxc_ao_reim,  # None if not DFT
                    'fock_gxc_ao_imre': fock_gxc_ao_imre,  # None if not DFT
                }
                if (f == 0):
                    tot_rhs_mo = np.concatenate((rhs_red.real, rhs_red.imag))
                else:
                    tot_rhs_mo = np.concatenate(
                        (tot_rhs_mo, rhs_red.real, rhs_red.imag))

            # save RHS in distributed array
            dist_cphf_rhs_re = []
            dist_cphf_rhs_im = []
            for k in range(dof_red):
                if self.rank == mpi_master():
                    cphf_rhs_k_re = rhs_red[k].real.reshape(nocc * nvir)
                    cphf_rhs_k_im = rhs_red[k].imag.reshape(nocc * nvir)
                else:
                    cphf_rhs_k_re = None
                    cphf_rhs_k_im = None
                dist_cphf_rhs_re.append(DistributedArray(cphf_rhs_k_re, self.comm, root=mpi_master()))
                dist_cphf_rhs_im.append(DistributedArray(cphf_rhs_k_im, self.comm, root=mpi_master()))

            dist_cphf_rhs.extend(dist_cphf_rhs_re + dist_cphf_rhs_im)


        if self.rank == mpi_master():
            orbrsp_rhs['cphf_rhs'] = tot_rhs_mo
            orbrsp_rhs['dist_cphf_rhs'] = dist_cphf_rhs
            return orbrsp_rhs
        else:
            return {
                'dist_cphf_rhs': dist_cphf_rhs,
            }

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing the orbrsp RHS for '
            valstr += '{} frequencies: '.format(len(self.frequencies))
            valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_rhs_real(self, molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict, lr_results):
        """
        Computes the right-hand side (RHS) of the Real polarizability
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

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Run the solver => in parent class

        loop_start_time = tm.time()

        nocc = molecule.number_of_alpha_electrons()
        nao = basis.get_dimensions_of_basis()

        # number of vector components
        dof = len(self.vector_components)

        if self.rank == mpi_master():
            # check if response vectors exist for desired frequency of gradient
            polgrad_sanity_check(self, self.flag, lr_results)

            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']  # only alpha part

            # MO coefficients
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            # reduce dimensions of RHS to unique operator component combinations
            xy_pairs = [(x,y) for x in range(3) for y in range(x,3)]
            dof_red = len(xy_pairs)
        else:
            density = None
            mo = None
            dof_red = None

        mo = self.comm.bcast(mo, root=mpi_master())
        density = self.comm.bcast(density, root=mpi_master())
        dof_red = self.comm.bcast(dof_red, root=mpi_master())

        # TODO: double check dft_dict['gs_density']
        # TODO: maybe rename to mol_grid
        molgrid = dft_dict['molgrid']
        gs_density = [density]

        orbrsp_rhs = {}

        dist_cphf_rhs = []

        for f, w in enumerate(self.frequencies):

            if self.rank == mpi_master():
                self.ostream.print_info(
                    'Building RHS for w = {:4.3f}'.format(w))
                self.ostream.flush()

            full_vec = [
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ]

            if self.rank == mpi_master():

                # Note: polorbitalresponse uses r instead of mu for dipole operator
                for idx in range(len(full_vec)):
                    full_vec[idx] *= -1.0

                # extract the excitation and de-excitation components
                # from the full solution vector.
                sqrt2 = np.sqrt(2.0)
                exc_vec = (1.0 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1.0 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                # construct plus/minus combinations of excitation and
                # de-excitation part
                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # transform to AO basis: mi,xia,na->xmn
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(dof)
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(dof)
                ])

                # turn them into a list for AODensityMatrix
                xpmy_ao_list = list(x_plus_y_ao) + list(x_minus_y_ao)

                # calculate symmetrized unrelaxed one-particle density matrix
                unrel_dm_ao, dm_oo, dm_vv = self.calculate_unrel_dm(molecule, scf_tensors,
                                                      x_plus_y, x_minus_y)
                # create lists
                dm_ao_list = list(unrel_dm_ao.reshape(dof**2, nao, nao))

                # density matrix for RHS
                dm_ao_rhs_list = dm_ao_list + xpmy_ao_list

                if self._dft:
                    # construct density matrices for E[3] term:
                    perturbed_dm_ao_list, zero_dm_ao_list = self.construct_dft_e3_dm_real(x_minus_y_ao)
            else:
                #dof = None
                dm_ao_rhs_list = None

                if self._dft:
                    perturbed_dm_ao_list = None
                    zero_dm_ao_list = None

            dm_ao_rhs_list = self.comm.bcast(dm_ao_rhs_list, root=mpi_master())

            if self._dft:
                perturbed_dm_ao_list = self.comm.bcast(perturbed_dm_ao_list, root=mpi_master())
                zero_dm_ao_list = self.comm.bcast(zero_dm_ao_list, root=mpi_master())

            if self._dft:
                fock_gxc_ao = []
                for i_mat in range(len(zero_dm_ao_list)):
                    fock_gxc_ao.append(zero_dm_ao_list[i_mat].copy())

                xc_drv = XCIntegrator()
                xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis,
                                          perturbed_dm_ao_list, zero_dm_ao_list,
                                          gs_density, molgrid,
                                          self.xcfun.get_func_label(), "qrf")

                for idx in range(len(fock_gxc_ao)):
                    fock_gxc_ao[idx] = self.comm.reduce(fock_gxc_ao[idx], root=mpi_master())
            else:
                fock_gxc_ao = None

            # vector-related components to general Fock matrix
            # (not 1PDM part)
            fock_ao_rhs = self._comp_lr_fock(dm_ao_rhs_list, molecule, basis,
                               eri_dict, dft_dict, pe_dict, self.profiler)

            # calculate the RHS
            if self.rank == mpi_master():
                # extract the 1PDM contributions
                fock_ao_rhs_1pdm = np.zeros((dof**2, nao, nao))
                for i in range(dof**2):
                    fock_ao_rhs_1pdm[i] = fock_ao_rhs[i]

                # transform to MO basis: mi,xmn,na->xia
                fock_mo_rhs_1pdm = np.array([
                    np.linalg.multi_dot([mo_occ.T, fock_ao_rhs_1pdm[x], mo_vir])
                    for x in range(dof**2)
                ])

                # extract the x_plus_y and x_minus_y contributions
                fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
                fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
                for i in range(dof):
                    fock_ao_rhs_x_plus_y[i] = fock_ao_rhs[ dof**2 + i]
                    fock_ao_rhs_x_minus_y[i] = fock_ao_rhs[ dof**2 + dof + i]

                # calculate 2-particle density matrix contribution
                fock_mo_rhs_2pdm = self.calculate_rhs_2pdm_contrib(molecule, scf_tensors,
                    x_plus_y_ao, x_minus_y_ao, fock_ao_rhs_x_plus_y, fock_ao_rhs_x_minus_y)

                # calculate dipole contribution
                rhs_dipole_contrib = self.calculate_rhs_dipole_contrib(
                    molecule, basis, scf_tensors, x_minus_y)

                # sum RHS contributions
                rhs_mo = fock_mo_rhs_1pdm + fock_mo_rhs_2pdm + rhs_dipole_contrib

                # add DFT E[3] contribution to the RHS
                if self._dft:
                    gxc_ao = np.zeros((dof**2, nao, nao))

                    for i in range(dof**2):
                        gxc_ao[i] = fock_gxc_ao[2 * i]

                    # mi,xmn,na->xia
                    gxc_mo = np.array([
                        np.linalg.multi_dot([mo_occ.T, gxc_ao[x], mo_vir])
                        for x in range(dof**2)
                    ])
                    # different factor compared to TDDFT orbital response
                    # because here vectors are scaled by 1/sqrt(2)
                    rhs_mo += 0.5 * gxc_mo

            self.profiler.stop_timer('RHS')

            if self.rank == mpi_master():
                # reduce dimensions of RHS to unique operator component combinations
                rhs_red = []
                rhs_tmp = rhs_mo.reshape(dof, dof, nocc, nvir).copy()

                for x, y in xy_pairs:
                    rhs_red.append(rhs_tmp[x, y].copy())
                rhs_red = np.array(rhs_red)

                # TODO purge dict
                orbrsp_rhs[(w)] = {
                    #'cphf_rhs': rhs_red,
                    'dm_oo': dm_oo,
                    'dm_vv': dm_vv,
                    #'x_plus_y_ao': x_plus_y_ao,
                    #'x_minus_y_ao': x_minus_y_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'fock_ao_rhs': fock_ao_rhs,
                    'fock_gxc_ao': fock_gxc_ao,  # None if not DFT
                }
                if (f == 0):
                    tot_rhs_mo = rhs_red
                else:
                    tot_rhs_mo = np.append(tot_rhs_mo, rhs_red, axis=0)

            # save RHS in distributed array
            for k in range(dof_red):
                if self.rank == mpi_master():
                    cphf_rhs_k = rhs_red[k].reshape(nocc * nvir)
                else:
                    cphf_rhs_k = None
                dist_cphf_rhs.append(DistributedArray(cphf_rhs_k, self.comm, root=mpi_master()))

        if self.rank == mpi_master():
            orbrsp_rhs['dist_cphf_rhs'] = dist_cphf_rhs
            if not self.use_subspace_solver:
                orbrsp_rhs['cphf_rhs'] = tot_rhs_mo
            return orbrsp_rhs
        else:
            return {
                'dist_cphf_rhs': dist_cphf_rhs,
            }

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing the orbrsp RHS for '
            valstr += '{} frequencies: '.format(len(self.frequencies))
            valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def calculate_unrel_dm(self, molecule, scf_tensors, x_plus_y, x_minus_y):
        """
        Calculates the symmetrized unrelaxed one-particle density matrix
        in AO basis.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param x_plus_y_ao:
            The X+Y response vectors in AO basis.
        :param x_minus_y_ao:
            The X-Y response vectors in AO basis.

        :return dm_oo:
            Occ/Occ block of unrelaxed one-particle density matrix in MO basis.
        :return dm_vv:
            Vir/vir block of unrelaxed one-particle density matrix in MO basis.
        :return unrel_dm_ao:
            Unrelaxed one-particle density matrix in AO basis.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        # number of AOs
        nao = mo.shape[0]

        # determine data type of RHS
        if self.is_complex:
            rhs_dt = np.dtype('complex128')
        else:
            rhs_dt = np.dtype('float64')

        # calculate the symmetrized unrelaxed one-particle density matrix
        # in MO basis
        dm_oo = np.zeros((dof, dof, nocc, nocc), dtype = rhs_dt)
        dm_vv = np.zeros((dof, dof, nvir, nvir), dtype = rhs_dt)

        for x in range(dof):
            for y in range(x, dof):
                dm_vv[x, y] = 0.25 * (
                    # xib,yia->xyab
                    np.linalg.multi_dot([x_plus_y[y].T, x_plus_y[x]])
                    # xib,yia->xyab
                    + np.linalg.multi_dot([x_minus_y[y].T, x_minus_y[x]])
                    # yib,xia->xyab
                    + np.linalg.multi_dot([x_plus_y[x].T, x_plus_y[y]])
                    # yib,xia->xyab
                    + np.linalg.multi_dot([x_minus_y[x].T, x_minus_y[y]]))

                dm_oo[x, y] = -0.25 * (
                    # xja,yia->xyij
                    np.linalg.multi_dot([x_plus_y[x], x_plus_y[y].T])
                    # xja,yia->xyij
                    + np.linalg.multi_dot([x_minus_y[x], x_minus_y[y].T])
                    # yja,xia->xyij
                    + np.linalg.multi_dot([x_plus_y[y], x_plus_y[x].T])
                    # yja,xia->xyij
                    + np.linalg.multi_dot([x_minus_y[y], x_minus_y[x].T]))

                if (y != x):
                    dm_vv[y,x] = dm_vv[x,y]
                    dm_oo[y,x] = dm_oo[x,y]

        # transform to AO basis: mi,xia,na->xmn
        unrel_dm_ao = np.zeros((dof, dof, nao, nao), dtype = rhs_dt)
        for x in range(dof):
            for y in range(x, dof):
                unrel_dm_ao[x, y] = (
                        # mi,xyij,nj->xymn
                        np.linalg.multi_dot([mo_occ, dm_oo[x, y], mo_occ.T])
                        # ma,xyab,nb->xymn
                        + np.linalg.multi_dot([mo_vir, dm_vv[x, y], mo_vir.T]))

                if (y != x):
                    unrel_dm_ao[y, x] = unrel_dm_ao[x, y]

        return unrel_dm_ao, dm_oo, dm_vv

    def calculate_rhs_2pdm_contrib(self, molecule, scf_tensors, x_plus_y_ao, x_minus_y_ao,
                                  fock_ao_rhs_x_plus_y, fock_ao_rhs_x_minus_y):
        """
        Calculates the 2-particle density matrix contribution to the RHS.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param x_plus_y_ao:
            The X+Y response vectors in AO basis.
        :param x_minus_y_ao:
            The X-Y response vectors in AO basis.
        :param fock_ao_rhs_x_plus_y:
            Fock matrix from unrel. DM and X+Y response vectors in AO basis.
        :param fock_ao_rhs_x_minus_y:
            Fock matrix from unrel. DM and X-Y response vectors in AO basis.

        :return fock_mo_rhs_2pdm:
            The Fock 2-particle DM contribution to RHS in MO basis.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # overlap
        ovlp = scf_tensors['S']

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        # determine data type of RHS
        if self.is_complex:
            rhs_dt = np.dtype('complex128')
        else:
            rhs_dt = np.dtype('float64')

        fock_mo_rhs_2pdm = np.zeros((dof, dof, nocc, nvir), dtype = rhs_dt)

        for x in range(dof):
            for y in range(x, dof):
                # xmt,ymc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[x], x_plus_y_ao[y]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # xmt,ymc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[x], x_minus_y_ao[y]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # ymt,xmc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[y], x_plus_y_ao[x]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # ymt,xmc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[y], x_minus_y_ao[x]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # xtm,ymc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[x].T, x_plus_y_ao[y]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # xtm,ymc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[x].T, x_minus_y_ao[y]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # ytm,xmc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[y].T, x_plus_y_ao[x]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # ytm,xmc->xytc
                tmp = np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[y].T, x_minus_y_ao[x]])
                # cl,ti,la,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_occ.T, tmp, ovlp, mo_vir])
                # xmt,ycm->xyct
                tmp = np.linalg.multi_dot(
                    [x_plus_y_ao[y], fock_ao_rhs_x_plus_y[x].T])
                # cl,li,ta,xyct->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, ovlp.T, tmp, mo_vir])
                # xmt,ycm->xyct
                tmp = np.linalg.multi_dot(
                    [x_minus_y_ao[y], fock_ao_rhs_x_minus_y[x].T])
                # cl,li,ta,xyct->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, ovlp.T, tmp, mo_vir])
                # ymt,xcm->xyct
                tmp = np.linalg.multi_dot(
                    [x_plus_y_ao[x], fock_ao_rhs_x_plus_y[y].T])
                # cl,li,ta,xyct->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, ovlp.T, tmp, mo_vir])
                # ymt,xcm->xyct
                tmp = np.linalg.multi_dot(
                    [x_minus_y_ao[x], fock_ao_rhs_x_minus_y[y].T])
                # cl,li,ta,xyct->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_occ.T, ovlp.T, tmp, mo_vir])
                # xtm,ycm->xytc
                tmp = np.linalg.multi_dot(
                    [x_plus_y_ao[y], fock_ao_rhs_x_plus_y[x]]).T
                # cl,li,ta,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_vir.T, tmp, ovlp, mo_occ]).T
                # xtm,ymc->xytc
                tmp = np.linalg.multi_dot(
                    [x_minus_y_ao[y], fock_ao_rhs_x_minus_y[x]]).T
                # cl,li,ta,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_vir.T, tmp, ovlp, mo_occ]).T
                # ytm,xcm->xytc
                tmp = np.linalg.multi_dot(
                    [x_plus_y_ao[x], fock_ao_rhs_x_plus_y[y]]).T
                # cl,li,ta,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += -1.0 * np.linalg.multi_dot(
                    [mo_vir.T, tmp, ovlp, mo_occ]).T
                # ytm,xcm->xytc
                tmp = np.linalg.multi_dot(
                    [x_minus_y_ao[x], fock_ao_rhs_x_minus_y[y]]).T
                # cl,li,ta,xytc->xyia
                fock_mo_rhs_2pdm[x, y] += np.linalg.multi_dot(
                    [mo_vir.T, tmp, ovlp, mo_occ]).T

                # save value in lower off-diagonal block
                if (y != x):
                    fock_mo_rhs_2pdm[y, x] = fock_mo_rhs_2pdm[x, y]

        fock_mo_rhs_2pdm = 0.25 * fock_mo_rhs_2pdm.reshape(
            dof**2, nocc, nvir)

        return fock_mo_rhs_2pdm

    def calculate_rhs_dipole_contrib(self, molecule, basis, scf_tensors,
                                     x_minus_y):
        """
        Calculates the dipole contribution to the RHS.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param x_minus_y:
            The X-Y response vectors.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        # number of AOs
        nao = mo.shape[0]

        # determine data type of RHS
        if self.is_complex:
            rhs_dt = np.dtype('complex128')
        else:
            rhs_dt = np.dtype('float64')

        # get the dipole integrals in AO basis
        dipole_mats = compute_electric_dipole_integrals(molecule, basis, [0.0, 0.0, 0.0])

        dipole_ints_ao = np.zeros((dof, nao, nao))
        k = 0
        if 'x' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[0]
            k += 1
        if 'y' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[1]
            k += 1
        if 'z' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[2]

        # transform to MO basis (oo and vv blocks only)
        dipole_ints_oo = np.array([
            np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
            for x in range(dof)
        ])
        dipole_ints_vv = np.array([
            np.linalg.multi_dot([mo_vir.T, dipole_ints_ao[x], mo_vir])
            for x in range(dof)
        ])

        rhs_dipole_contrib = np.zeros((dof, dof, nocc, nvir), dtype=rhs_dt)

        # calculate dipole contributions to the RHS
        for x in range(dof):
            for y in range(x, dof):
                rhs_dipole_contrib[x, y] = (
                    0.5 * (np.linalg.multi_dot( # xja,yji->xyia
                    [dipole_ints_oo[y].T, x_minus_y[x]])
                    + np.linalg.multi_dot( # yja,xji->xyia
                    [dipole_ints_oo[x], x_minus_y[y]]))
                    - 0.5 * (np.linalg.multi_dot( # xib,yab->xyia
                    [x_minus_y[x], dipole_ints_vv[y]])
                    + np.linalg.multi_dot( # yib,xab->xyia
                    [x_minus_y[y], dipole_ints_vv[x].T])))
                
                if (y != x):
                    rhs_dipole_contrib[y, x] = rhs_dipole_contrib[x, y]

        rhs_dipole_contrib = rhs_dipole_contrib.reshape(
                    dof**2, nocc, nvir)

        return rhs_dipole_contrib

    def construct_dft_e3_dm_real(self, x_minus_y_ao):
        """
        Constructs the density matrices for E[3] term
        for the real RHS.

        :param x_minus_y_ao:
            the X-Y response vectors in AO basis.

        :return perturbed_dm_ao:
            Perturbed density matrix as an AODensityMatrix.
        :return zero_dm_ao:
            Empty matrix same size as perturbed density matrix as an AODensityMatrix.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # sanity check: should not be carried out in parallel
        if self.rank != mpi_master():
            return

        # XCIntegrator expects a DM with real and imaginary part,
        # so we set the imaginary part to zero.

        perturbed_dm_ao_list = []
        zero_dm_ao_list = []

        for x in range(dof):
            for y in range(dof):
                perturbed_dm_ao_list.extend([
                    x_minus_y_ao[x], 0 * x_minus_y_ao[x],
                    x_minus_y_ao[y], 0 * x_minus_y_ao[y]
                ])
                zero_dm_ao_list.extend(
                    [0 * x_minus_y_ao[x], 0 * x_minus_y_ao[y]])

        return perturbed_dm_ao_list, zero_dm_ao_list

    def construct_dft_e3_dm_complex(self, x_minus_y_ao):
        """
        Constructs the density matrices for E[3] term
        for the complex RHS.

        :param x_minus_y_ao:
            the X-Y response vectors in AO basis.

        :return perturbed_dm_ao_rere:
            The perturbed density matrix from Re/Re parts of the X-Y response vectors as an AODensityMatrix.
        :return perturbed_dm_ao_imim:
            The perturbed density matrix from Im/Im parts of the X-Y response vectors as an AODensityMatrix.
        :return perturbed_dm_ao_reim:
            The perturbed density matrix from Re/Im parts of the X-Y response vectors as an AODensityMatrix.
        :return perturbed_dm_ao_imre:
            The perturbed density matrix from Im/Re parts of the X-Y response vectors as an AODensityMatrix.
        :return zero_dm_ao:
            Empty matrix same size as perturbed density matrices as an AODensityMatrix.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # sanity check: should not be carried out in parallel
        if self.rank != mpi_master():
            return

        perturbed_dm_ao_list_rere = []
        perturbed_dm_ao_list_imim = []
        perturbed_dm_ao_list_reim = []
        perturbed_dm_ao_list_imre = []
        zero_dm_ao_list = []

        for x in range(dof):
            for y in range(dof):
                perturbed_dm_ao_list_rere.extend([
                    np.array(x_minus_y_ao[x].real),
                    np.array(0 * x_minus_y_ao[x].real),
                    np.array(x_minus_y_ao[y].real),
                    np.array(0 * x_minus_y_ao[y].real)])
                perturbed_dm_ao_list_imim.extend([
                    np.array(x_minus_y_ao[x].imag),
                    np.array(0 * x_minus_y_ao[x].imag),
                    np.array(x_minus_y_ao[y].imag),
                    np.array(0 * x_minus_y_ao[y].imag)])
                # complex cross-terms
                perturbed_dm_ao_list_reim.extend([
                    np.array(x_minus_y_ao[x].real),
                    np.array(0 * x_minus_y_ao[x].real),
                    np.array(x_minus_y_ao[y].imag),
                    np.array(0 * x_minus_y_ao[y].imag)])
                perturbed_dm_ao_list_imre.extend([
                    np.array(x_minus_y_ao[x].imag),
                    np.array(0 * x_minus_y_ao[x].imag),
                    np.array(x_minus_y_ao[y].real),
                    np.array(0 * x_minus_y_ao[y].real)])

                zero_dm_ao_list.extend([
                    np.array(0 * x_minus_y_ao[x].real),
                    np.array(0 * x_minus_y_ao[y].real)])

        perturbed_dm_ao_rere = AODensityMatrix(perturbed_dm_ao_list_rere,
                                          denmat.rest)
        perturbed_dm_ao_imim = AODensityMatrix(perturbed_dm_ao_list_imim,
                                          denmat.rest)
        perturbed_dm_ao_reim = AODensityMatrix(perturbed_dm_ao_list_reim,
                                          denmat.rest)
        perturbed_dm_ao_imre = AODensityMatrix(perturbed_dm_ao_list_imre,
                                          denmat.rest)

        # corresponds to rho^{omega_b,omega_c} in quadratic response,
        # which is zero for orbital response
        zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)

        return (perturbed_dm_ao_rere, perturbed_dm_ao_imim,
               perturbed_dm_ao_reim, perturbed_dm_ao_imre, zero_dm_ao)

    def set_dft_fmat_factor_and_type_real(self, fock_ao_rhs, zero_dm_ao):#fock_gxc_ao):
        """
        Sets of scale factor and fock type for Fock matrices to compute DFT E[3] term
        for real RHS.

        :param fock_ao_rhs:
            The general Fock matrix in AO basis.
        :param zero_dm_ao:
            Empty AODensityMatrix in same dimensions as the Fock matrices.
        :param fock_gxc_ao:
            Fock matrix for g^xc term in AO basis.

        :return fock_ao_rhs:
            Input param modified with set scale factor and fock type.
        :return fock_gxc_ao:
            Input param modified with set scale factor and fock type.
        """

        # degrees of freedom
        dof = len(self.vector_components)

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

        return fock_ao_rhs, fock_gxc_ao

    def set_dft_fmat_factor_and_type_complex(self, fock_ao_rhs_real, fock_ao_rhs_imag,
                                          zero_dm_ao):
        """
        Sets of scale factor and fock type for Fock matrices to compute DFT E[3] term
        for complex RHS.

        :param fock_ao_rhs_real/imag:
            The general Fock matrix in AO basis.
        :param zero_dm_ao:
            Empty AODensityMatrix in same dimensions as the Fock matrices.
        :param fock_gxc_ao:
            Fock matrix for g^xc term in AO basis.

        :return fock_ao_rhs_real/imag:
            Input param modified with set scale factor and fock type.
        :return fock_gxc_ao_rere/imim/reim/imre:
            Input param modified with set scale factor and fock type.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # Fock matrix for computing the DFT E[3] term g^xc
        fock_gxc_ao_rere = AOFockMatrix(zero_dm_ao)
        fock_gxc_ao_imim = AOFockMatrix(zero_dm_ao)
        fock_gxc_ao_reim = AOFockMatrix(zero_dm_ao)
        fock_gxc_ao_imre = AOFockMatrix(zero_dm_ao)

        if self.xcfun.is_hybrid():
            fact_xc = self.xcfun.get_frac_exact_exchange()
            for ifock in range(fock_ao_rhs_real.number_of_fock_matrices()):
                fock_ao_rhs_real.set_scale_factor(fact_xc, ifock)
                fock_ao_rhs_imag.set_scale_factor(fact_xc, ifock)
            for ifock in range(fock_gxc_ao_rere.number_of_fock_matrices()):
                fock_gxc_ao_rere.set_scale_factor(fact_xc, ifock)
                fock_gxc_ao_imim.set_scale_factor(fact_xc, ifock)
                fock_gxc_ao_reim.set_scale_factor(fact_xc, ifock)
                fock_gxc_ao_imre.set_scale_factor(fact_xc, ifock)
                fock_gxc_ao_rere.set_fock_type(fockmat.rgenjkx, ifock)
                fock_gxc_ao_imim.set_fock_type(fockmat.rgenjkx, ifock)
                fock_gxc_ao_reim.set_fock_type(fockmat.rgenjkx, ifock)
                fock_gxc_ao_imre.set_fock_type(fockmat.rgenjkx, ifock)
            for ifock in range(dof**2):
                fock_ao_rhs_real.set_fock_type(fockmat.restjkx, ifock)
                fock_ao_rhs_imag.set_fock_type(fockmat.restjkx, ifock)
            for ifock in range(dof**2, dof**2 + 2 * dof):
                fock_ao_rhs_real.set_fock_type(fockmat.rgenjkx, ifock)
                fock_ao_rhs_imag.set_fock_type(fockmat.rgenjkx, ifock)
        else:
            for ifock in range(dof**2):
                fock_ao_rhs_real.set_fock_type(fockmat.restj, ifock)
                fock_ao_rhs_imag.set_fock_type(fockmat.restj, ifock)
            for ifock in range(dof**2, dof**2 + 2 * dof):
                fock_ao_rhs_real.set_fock_type(fockmat.rgenj, ifock)
                fock_ao_rhs_imag.set_fock_type(fockmat.rgenj, ifock)
            for ifock in range(fock_gxc_ao_rere.number_of_fock_matrices()):
                fock_gxc_ao_rere.set_fock_type(fockmat.rgenj, ifock)
                fock_gxc_ao_imim.set_fock_type(fockmat.rgenj, ifock)
                fock_gxc_ao_reim.set_fock_type(fockmat.rgenj, ifock)
                fock_gxc_ao_imre.set_fock_type(fockmat.rgenj, ifock)

        return (fock_ao_rhs_real, fock_ao_rhs_imag, fock_gxc_ao_rere,
                fock_gxc_ao_imim, fock_gxc_ao_reim, fock_gxc_ao_imre)

    def integrate_gxc_real(self, molecule, basis, molgrid, gs_density,
                           zero_dm_ao, perturbed_dm_ao, fock_gxc_ao):
        """
        Drives the integration of the real XC Fock thingything.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param molgrid:
            The grid level for the DFT calculation.
        :param gs_density:
            The ground-state density.
        :param zero_dm_ao:
            AODensityMatrix with zeros.
        :param perturbed_dm_ao:
            Perturbed density matrix in AO basis.
        :param fock_gxc_ao:
            The g^xc Fock matrix in AO basis.

        :return fock_gxc_ao:
            The g^xc Fock matrix after integration.
        """

        xc_drv = XCIntegrator(self.comm)
        xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis,
                                  perturbed_dm_ao, zero_dm_ao,
                                  gs_density, molgrid,
                                  self.xcfun.get_func_label(), "qrf")
        return fock_gxc_ao

    def integrate_gxc_complex(self, molecule, basis, molgrid, gs_density, zero_dm_ao,
                              perturbed_dm_ao_list, fock_gxc_ao_list):
        """
        Drives the integration of the complex XC Fock thingything.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param molegrid:
            The grid level for the DFT calculation.
        :param gs_density:
            The ground-state density.
        :param zero_dm_ao:
            AODensityMatrix with zeros.
        :param perturbed_dm_ao_list:
            List with the rere/imim/reim/imre parts of the perturbed
            density matrix in AO basis.
        :param fock_gxc_ao_list:
            List with the rere/imim/reim/imre parts of the g^xc Fock matrix in AO basis.

        :return fock_gxc_ao_rere/imim/reim/imre:
            The g^xc Fock matrix after integration.
        """

        # unpack perturbed dm list
        perturbed_dm_ao_rere = perturbed_dm_ao_list[0]
        perturbed_dm_ao_imim = perturbed_dm_ao_list[1]
        perturbed_dm_ao_reim = perturbed_dm_ao_list[2]
        perturbed_dm_ao_imre = perturbed_dm_ao_list[3]

        # unpack XC Fock matrix list
        fock_gxc_ao_rere = fock_gxc_ao_list[0]
        fock_gxc_ao_imim = fock_gxc_ao_list[1]
        fock_gxc_ao_reim = fock_gxc_ao_list[2]
        fock_gxc_ao_imre = fock_gxc_ao_list[3]

        if not self.xcfun.is_hybrid():
            for ifock in range(fock_gxc_ao_rere.number_of_fock_matrices()):
                fock_gxc_ao_rere.scale(2.0, ifock)
                fock_gxc_ao_imim.scale(2.0, ifock)
                fock_gxc_ao_reim.scale(2.0, ifock)
                fock_gxc_ao_imre.scale(2.0, ifock)
               
        xc_drv = XCIntegrator(self.comm)
        xc_drv.integrate_kxc_fock(fock_gxc_ao_rere, molecule, basis,
                                  perturbed_dm_ao_rere, zero_dm_ao,
                                  gs_density, molgrid,
                                  self.xcfun.get_func_label(), "qrf")
        xc_drv.integrate_kxc_fock(fock_gxc_ao_imim, molecule, basis,
                                  perturbed_dm_ao_imim, zero_dm_ao,
                                  gs_density, molgrid,
                                  self.xcfun.get_func_label(), "qrf")
        xc_drv.integrate_kxc_fock(fock_gxc_ao_reim, molecule, basis,
                                  perturbed_dm_ao_reim, zero_dm_ao,
                                  gs_density, molgrid,
                                  self.xcfun.get_func_label(), "qrf")
        xc_drv.integrate_kxc_fock(fock_gxc_ao_imre, molecule, basis,
                                  perturbed_dm_ao_imre, zero_dm_ao,
                                  gs_density, molgrid,
                                  self.xcfun.get_func_label(), "qrf")

        return fock_gxc_ao_rere, fock_gxc_ao_imim, fock_gxc_ao_reim, fock_gxc_ao_imre

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

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results from the linear response calculation.
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

        if self.rank == mpi_master():
            # MO coefficients
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]

            # number of atomic orbitals
            nao = basis.get_dimensions_of_basis()

            # number of vector components
            dof = len(self.vector_components)
            xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]
            dof_red = len(xy_pairs)

            # number of frequencies
            n_freqs = len(self.frequencies)

            if not self.use_subspace_solver:
                all_cphf_red = self.cphf_results['cphf_ov']
                all_cphf_red = all_cphf_red.reshape(n_freqs, dof_red, nocc * nvir)
        else:
            xy_pairs = None
            dof_red = None
            dof = None
            nocc = None
            nvir = None

        xy_pairs = self.comm.bcast(xy_pairs, root=mpi_master())
        dof_red = self.comm.bcast(dof_red, root=mpi_master())
        dof = self.comm.bcast(dof, root=mpi_master())
        nocc = self.comm.bcast(nocc , root=mpi_master())
        nvir = self.comm.bcast(nvir , root=mpi_master())

        # timings
        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):

            full_vec = [
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ]

            # cphf subspace solver returns the solution as distributed array
            if self.use_subspace_solver:
                if self.rank == mpi_master():
                    cphf_ov = np.zeros((dof, dof, nocc * nvir))
                else:
                    cphf_ov = None

                for idx, xy in enumerate(xy_pairs):
                    # get lambda multipliers from distributed array
                    tmp_cphf_ov = self.cphf_results['dist_cphf_ov'][dof_red * f + idx].get_full_vector()

                    if self.rank == mpi_master():
                        x = xy[0]
                        y = xy[1]
                        cphf_ov[x, y] += tmp_cphf_ov

                        if (y != x):
                            cphf_ov[y, x] += cphf_ov[x, y]

            if self.rank == mpi_master():

                self.ostream.print_info('Building omega for w = {:4.3f}'.format(w))
                self.ostream.flush()

                # Note: polorbitalresponse uses r instead of mu for dipole operator
                for idx in range(len(full_vec)):
                    full_vec[idx] *= -1.0

                # get fock matrices from cphf_results
                fock_ao_rhs = self.cphf_results[w]['fock_ao_rhs']
                fock_gxc_ao = self.cphf_results[w]['fock_gxc_ao']
                dm_oo = self.cphf_results[w]['dm_oo']
                dm_vv = self.cphf_results[w]['dm_vv']

                # note: conjugate gradient solver returns array in dictionary
                # and the coefficients are therefore imported differently
                # from subspace solver
                if not self.use_subspace_solver:
                    cphf_ov_red = all_cphf_red[f]
                    cphf_ov = np.zeros((dof, dof, nocc * nvir))

                    for idx, xy in enumerate(xy_pairs):
                        x = xy[0]
                        y = xy[1]

                        cphf_ov[x, y] = cphf_ov_red[idx]

                        if (y != x):
                            cphf_ov[y, x] += cphf_ov[x, y]
                           
                cphf_ov = cphf_ov.reshape(dof**2, nocc, nvir)

                # create response vectors in MO basis
                sqrt2 = np.sqrt(2.0)
                exc_vec = (1.0 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1.0 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # transform to AO basis: mi,xia,na->xmn
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(x_plus_y.shape[0])
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(x_minus_y.shape[0])
                ])

                # calculate dipole contribution to omega
                omega_dipole_contrib_ao = self.calculate_omega_dipole_contrib(
                    molecule, basis, scf_tensors, x_minus_y)

                # calculate the density matrices, alpha block only
                D_occ = np.matmul(mo_occ, mo_occ.T)

                # construct fock_lambda (or fock_cphf)
                # transform to AO basis: mi,xia,na->xmn
                cphf_ao = np.array([
                    np.linalg.multi_dot([mo_occ, cphf_ov[xy], mo_vir.T])
                    for xy in range(dof**2)
                ])
                cphf_ao_list = [cphf_ao[x] for x in range(dof**2)]
            else:
                cphf_ao_list = None

            cphf_ao_list = self.comm.bcast(cphf_ao_list, root=mpi_master())

            # TODO: what has to be on MPI master and what not?
            fock_cphf = self._comp_lr_fock(cphf_ao_list, molecule, basis,
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
                omega = np.zeros((dof, dof, nao, nao))

                # construct epsilon density matrix
                epsilon_dm_ao = self.calculate_epsilon_dm(molecule, scf_tensors,
                                                          dm_oo, dm_vv, cphf_ov)

                for m in range(dof):
                    for n in range(m, dof):
                        # TODO: move outside for-loop when all Fock matrices can be
                        # extracted into a numpy array at the same time.

                        # since the excitation vector is not symmetric,
                        # we need both the matrix (OO block in omega, and probably VO)
                        # and its transpose (VV, OV blocks)
                        # this comes from the transformation of the 2PDM contribution
                        # from MO to AO basis
                        fock_ao_rhs_1_m = fock_ao_rhs[ dof**2 + m]  # x_plus_y
                        fock_ao_rhs_2_m = fock_ao_rhs[ dof**2 + dof + m]  # x_minus_y

                        fock_ao_rhs_1_n = fock_ao_rhs[ dof**2 + n]  # x_plus_y
                        fock_ao_rhs_2_n = fock_ao_rhs[ dof**2 + dof + n]  # x_minus_y

                        fmat = (fock_cphf[m * dof + n] +
                                fock_cphf[m * dof + n].T +
                                fock_ao_rhs[m * dof + n])

                        # dof=3  (0,0), (0,1), (0,2); (1,0), (1,1), (1,2),
                        #        (2,0), (2,1), (2,2) * dof

                        # calculate the contributions from 2PDM and relaxed 1PDM
                        omega_1pdm_2pdm_contrib_mn = self.calculate_omega_1pdm_2pdm_contrib(
                            molecule, scf_tensors, x_plus_y_ao[m], x_plus_y_ao[n],
                            x_minus_y_ao[m], x_minus_y_ao[n], fock_ao_rhs_1_m,
                            fock_ao_rhs_2_m, fock_ao_rhs_1_n, fock_ao_rhs_2_n, fmat)

                        # sum contributions to omega
                        omega[m, n] = (epsilon_dm_ao[m, n] +
                                              omega_1pdm_2pdm_contrib_mn +
                                              omega_dipole_contrib_ao[m, n])
                        if self._dft:
                            omega_gxc_contrib = self.calculate_omega_gxc_contrib_real(
                                scf_tensors, fock_gxc_ao[2 * (m * dof + n)],
                                D_occ)

                            omega[m, n] += omega_gxc_contrib

                        # set values in lower off-diagonal
                        if (n != m):
                            omega[n, m] = omega[m, n]

                # TODO make distributed
                omega = omega.reshape(dof * dof, nao, nao)

                # save omega multipliers in cphf_results dictionary
                self.cphf_results[(w)]['omega_ao'] = omega
                # TODO remove
                #self.cphf_results[(w)]['lambda_ao'] = cphf_ao

        if self.rank == mpi_master():
            self.ostream.print_blank()
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

        if self.rank == mpi_master():
            # MO coefficients
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]

            # number of atomic orbitals
            nao = basis.get_dimensions_of_basis()

            # number of vector components
            dof = len(self.vector_components)
            xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]
            dof_red = len(xy_pairs)

            # number of frequencies
            n_freqs = len(self.frequencies)

            # get CPHF results from conj. gradient solver
            if not self.use_subspace_solver:
                all_cphf_red = self.cphf_results['cphf_ov']
                all_cphf_red = all_cphf_red.reshape(n_freqs, 2 * dof_red, nocc * nvir)
        else:
            xy_pairs = None
            dof_red = None
            dof = None
            nocc = None
            nvir = None

        xy_pairs = self.comm.bcast(xy_pairs, root=mpi_master())
        dof_red = self.comm.bcast(dof_red, root=mpi_master())
        dof = self.comm.bcast(dof, root=mpi_master())
        nocc = self.comm.bcast(nocc , root=mpi_master())
        nvir = self.comm.bcast(nvir , root=mpi_master())

        # timings
        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):

            full_vec = [
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ]

            if self.use_subspace_solver:
                if self.rank == mpi_master():
                    cphf_ov = np.zeros((dof, dof, nocc * nvir), dtype = np.dtype('complex128'))

                for idx, xy in enumerate(xy_pairs):
                    tmp_cphf_re = self.cphf_results['dist_cphf_ov'][
                        2 * dof_red * f + idx].get_full_vector()
                    tmp_cphf_im = self.cphf_results['dist_cphf_ov'][
                        2 * dof_red * f + dof_red + idx].get_full_vector()

                    if self.rank == mpi_master():
                        x = xy[0]
                        y = xy[1]
                        cphf_ov[x, y] += tmp_cphf_re + 1j * tmp_cphf_im
                        if (y != x):
                            cphf_ov[y, x] += cphf_ov[x, y]
                del tmp_cphf_re
                del tmp_cphf_im

            if self.rank == mpi_master():

                # Note: polorbitalresponse uses r instead of mu for dipole operator
                for idx in range(len(full_vec)):
                    full_vec[idx] *= -1.0

                self.ostream.print_info('Building omega for w = {:4.3f}'.format(w))
                self.ostream.flush()
               
                # get fock matrices from cphf_results
                fock_ao_rhs_real = self.cphf_results[w]['fock_ao_rhs_real']
                fock_ao_rhs_imag = self.cphf_results[w]['fock_ao_rhs_imag']
                fock_gxc_ao_rere = self.cphf_results[w]['fock_gxc_ao_rere']
                fock_gxc_ao_imim = self.cphf_results[w]['fock_gxc_ao_imim']
                fock_gxc_ao_reim = self.cphf_results[w]['fock_gxc_ao_reim']
                fock_gxc_ao_imre = self.cphf_results[w]['fock_gxc_ao_imre']
                dm_oo = self.cphf_results[w]['dm_oo']  # complex
                dm_vv = self.cphf_results[w]['dm_vv']  # complex

                # get Lagrangian lambda multipliers
                if not self.use_subspace_solver:
                    cphf_ov_red = all_cphf_red[f]
                    cphf_ov = np.zeros((dof, dof, nocc * nvir), dtype = np.dtype('complex128'))

                    tmp_cphf_ov = cphf_ov_red[:dof_red] + 1j * cphf_ov_red[dof_red:]

                    for idx, xy in enumerate(xy_pairs):
                        x = xy[0]
                        y = xy[1]

                        cphf_ov[x, y] = tmp_cphf_ov[idx]

                        if (y != x):
                            cphf_ov[y, x] += cphf_ov[x, y]

                cphf_ov = cphf_ov.reshape(dof**2, nocc, nvir)

                # extract the excitation and de-excitation components
                # from the full solution vector.
                sqrt2 = np.sqrt(2.0)
                exc_vec = (1.0 / sqrt2 *
                           np.array(full_vec)[:, :nocc * nvir].reshape(
                               dof, nocc, nvir))
                deexc_vec = (1.0 / sqrt2 *
                             np.array(full_vec)[:, nocc * nvir:].reshape(
                                 dof, nocc, nvir))

                x_plus_y = exc_vec + deexc_vec
                x_minus_y = exc_vec - deexc_vec

                # transform to AO basis: mi,xia,na->xmn
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y[x], mo_vir.T])
                    for x in range(x_plus_y.shape[0])
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y[x], mo_vir.T])
                    for x in range(x_minus_y.shape[0])
                ])

                # calculate dipole contribution to omega
                omega_dipole_contrib_ao = self.calculate_omega_dipole_contrib(
                    molecule, basis, scf_tensors, x_minus_y)

                # calculate the density matrices, alpha block only
                D_occ = np.matmul(mo_occ, mo_occ.T)

                # construct fock_lambda (or fock_cphf)
                # transform to AO basis: mi,xia,na->xmn
                cphf_ao = np.array([
                    np.linalg.multi_dot([mo_occ, cphf_ov[x], mo_vir.T])
                    for x in range(dof**2)
                ])

                cphf_ao_list_real = list(
                    np.array([cphf_ao[x].real for x in range(dof**2)]))
                cphf_ao_list_imag = list(
                    np.array([cphf_ao[x].imag for x in range(dof**2)]))
            else:
                cphf_ao_list_real = None
                cphf_ao_list_imag = None

            dof = self.comm.bcast(dof, root=mpi_master())

            cphf_ao_list_real = self.comm.bcast(cphf_ao_list_real, root=mpi_master())
            cphf_ao_list_imag = self.comm.bcast(cphf_ao_list_imag, root=mpi_master())

            # TODO: what has to be on MPI master and what not?
            fock_cphf_real = self._comp_lr_fock(cphf_ao_list_real, molecule,
                               basis, eri_dict, dft_dict, pe_dict,
                               self.profiler)
            fock_cphf_imag = self._comp_lr_fock(cphf_ao_list_imag, molecule,
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
                omega = np.zeros((dof, dof, nao, nao), dtype=np.dtype('complex128'))

                # construct epsilon density matrix
                epsilon_dm_ao = self.calculate_epsilon_dm(molecule, scf_tensors,
                                                          dm_oo, dm_vv, cphf_ov)

                for m in range(dof):
                    for n in range(m, dof):
                        # TODO: move outside for-loop when all Fock matrices can be
                        # extracted into a numpy array at the same time.

                        # Because the excitation vector is not symmetric,
                        # we need both the matrix (OO block in omega, and probably VO)
                        # and its transpose (VV, OV blocks)
                        # this comes from the transformation of the 2PDM contribution
                        # from MO to AO basis
                        fock_ao_rhs_1_m = (
                            fock_ao_rhs_real[dof**2 + m] +
                            1j * fock_ao_rhs_imag[dof**2 + m]
                        )  # x_plus_y
                        fock_ao_rhs_2_m = (
                            fock_ao_rhs_real[dof**2 + dof + m] +
                            1j *
                            fock_ao_rhs_imag[dof**2 + dof + m]
                        )  # x_minus_y

                        fock_ao_rhs_1_n = (
                            fock_ao_rhs_real[dof**2 + n] +
                            1j * fock_ao_rhs_imag[dof**2 + n]
                        )  # x_plus_y
                        fock_ao_rhs_2_n = (
                            fock_ao_rhs_real[dof**2 + dof + n] +
                            1j *
                            fock_ao_rhs_imag[dof**2 + dof + n]
                        )  # x_minus_y

                        fmat = ((fock_cphf_real[m * dof + n] +
                                 fock_cphf_real[m * dof + n].T +
                                 fock_ao_rhs_real[m * dof + n]) +
                                1j *
                                (fock_cphf_imag[m * dof + n] +
                                 fock_cphf_imag[m * dof + n].T +
                                 fock_ao_rhs_imag[m * dof + n]))

                        # dof=3  (0,0), (0,1), (0,2); (1,0), (1,1), (1,2),
                        #        (2,0), (2,1), (2,2) * dof

                        # compute the contributions from 2PDM and relaxed 1PDM
                        omega_1pdm_2pdm_contrib = self.calculate_omega_1pdm_2pdm_contrib(
                            molecule, scf_tensors, x_plus_y_ao[m], x_plus_y_ao[n],
                            x_minus_y_ao[m], x_minus_y_ao[n], fock_ao_rhs_1_m,
                            fock_ao_rhs_2_m, fock_ao_rhs_1_n, fock_ao_rhs_2_n, fmat)

                        # sum contributions to omega
                        omega[m, n] = (epsilon_dm_ao[m, n] +
                                              omega_1pdm_2pdm_contrib +
                                              omega_dipole_contrib_ao[m, n])

                        if self._dft:
                            fock_gxc_ao_mn_list = [
                                fock_gxc_ao_rere[2 * (m * dof + n)],
                                fock_gxc_ao_imim[2 * (m * dof + n)],
                                fock_gxc_ao_reim[2 * (m * dof + n)],
                                fock_gxc_ao_imre[2 * (m * dof + n)],
                            ]
                            omega_gxc_contrib = self.calculate_omega_gxc_contrib_complex(
                                fock_gxc_ao_mn_list, D_occ
                            )
                            omega[m, n] += omega_gxc_contrib

                        # set values in lower off-diagonal
                        if (n != m):
                            omega[n, m] = omega[m, n]

                # TODO distributed
                omega = omega.reshape(dof * dof, nao, nao)

                # save omega multipliers in cphf_results dictionary as list of dist. arrays
                self.cphf_results[(w)]['omega_ao'] = omega
                # TODO remove
                # self.cphf_results[(w)]['lambda_ao'] = cphf_ao

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing omega multipliers '
            valstr += 'for {} frequencies: '.format(n_freqs)
            valstr += '{:.6f} sec **'.format(tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def calculate_omega_dipole_contrib(self, molecule, basis, scf_tensors,
                                     x_minus_y):
        """
        Calculates the dipole contribution to the omega multipliers.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param x_minus_y:
            The X-Y response vectors.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()

        # number of AOs
        nao = mo.shape[0]

        # determine data type of the omega multipliers
        if self.is_complex:
            omega_dt = np.dtype('complex128')
        else:
            omega_dt = np.dtype('float64')

        # get dipole moment integrals
        dipole_mats = compute_electric_dipole_integrals(molecule, basis, [0.0, 0.0, 0.0])

        dipole_ints_ao = np.zeros((dof, nao, nao))
        k = 0
        if 'x' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[0]
            k += 1
        if 'y' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[1]
            k += 1
        if 'z' in self.vector_components:
            # Note: polorbitalresponse uses r instead of mu for dipole operator
            dipole_ints_ao[k] = -1.0*dipole_mats[2]

        # transform to MO basis (oo and ov blocks only)
        dipole_ints_oo = np.array([
            np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ])
            for x in range(dof)
        ])
        dipole_ints_ov = np.array([
            np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_vir])
            for x in range(dof)
        ])

        omega_dipole_contrib = np.zeros((dof, dof, nao, nao), dtype = omega_dt)
        for x in range(dof):
            for y in range(x, dof):
                tmp_oo = 0.5 * (np.linalg.multi_dot([ # xjc,yic->xyij
                    dipole_ints_ov[y], x_minus_y[x].T])
                    + np.linalg.multi_dot( # yjc,xic->xyij
                    [dipole_ints_ov[x], x_minus_y[y].T]))
                tmp_ov = 0.5 * (np.linalg.multi_dot([ # xka,yki->xyia
                    dipole_ints_oo[y].T, x_minus_y[x]])
                    + np.linalg.multi_dot( # yka,xki->xyia
                    [dipole_ints_oo[x].T, x_minus_y[y]]))
                tmp_vv = 0.5 * (np.linalg.multi_dot([ # xkb,yka->xyab
                    dipole_ints_ov[y].T, x_minus_y[x]])
                    + np.linalg.multi_dot( # ykb,xka->xyab
                    [dipole_ints_ov[x].T, x_minus_y[y]]))
                omega_dipole_contrib[x, y] = (
                    # mi,xyij,nj->xymn
                    np.linalg.multi_dot([mo_occ, tmp_oo, mo_occ.T]) +
                    # mi,xyia,na->xymn
                    np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]) +
                    # mi,xyia,na->xymn
                    np.linalg.multi_dot([mo_occ, tmp_ov, mo_vir.T]).T +
                    # ma,xyab,nb->xymn
                    np.linalg.multi_dot([mo_vir, tmp_vv, mo_vir.T]))
                if (y != x):
                    omega_dipole_contrib[y, x] = omega_dipole_contrib[x, y]
               
        return omega_dipole_contrib

    def calculate_epsilon_dm(self, molecule, scf_tensors, dm_oo, dm_vv, lambda_ov):
        """
        Calculates the epsilon density matrix for the omega multipliers

        :param molecule:
            The molecule.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param dm_oo:
            Occ/occ block of the one-particle density matrix.
        :param dm_vv:
            Vir/vir block of the one-particle density matrix.
        :param lambda_ov:
            The lambda multipliers in MO basis.

        :return epsilon_dm:
            The epsilon density matrix in AO basis.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        # number of AOs
        nao = mo.shape[0]

        # orbital energies
        mo_energies = scf_tensors['E_alpha']
        eocc = mo_energies[:nocc]
        evir = mo_energies[nocc:]
        eo_diag = np.diag(eocc)
        ev_diag = np.diag(evir)

        # determine data type of the density matrix
        if self.is_complex:
            epsilon_dt = np.dtype('complex128')
        else:
            epsilon_dt = np.dtype('float64')

        # construct epsilon density matrix
        epsilon_dm = np.zeros((dof, dof, nao, nao), dtype = epsilon_dt)
        epsilon_lambda= np.zeros((dof, dof, nao, nao), dtype = epsilon_dt)
        for x in range(dof):
            for y in range(x, dof):
                # mi,ii,xyij,nj->xymn
                epsilon_dm[x, y] = -1.0 * np.linalg.multi_dot(
                    [mo_occ, eo_diag, dm_oo[x, y], mo_occ.T])
                # ma,aa,xyab,nb->xymn
                epsilon_dm[x, y] -= np.linalg.multi_dot(
                    [mo_vir, ev_diag, dm_vv[x, y], mo_vir.T])
                # mi,ii,xyia,na->xymn
                epsilon_lambda[x, y] = np.linalg.multi_dot([
                    mo_occ, eo_diag,
                    lambda_ov.reshape(dof, dof, nocc, nvir)[x, y],
                    mo_vir.T
                ])
                if (y != x):
                    epsilon_dm[y, x] = epsilon_dm[x, y]

        # symmetrize (OV + VO)
        epsilon_dm -= (epsilon_lambda +
                          epsilon_lambda.transpose(0, 1, 3, 2))

        return epsilon_dm

    def calculate_omega_1pdm_2pdm_contrib(self, molecule, scf_tensors, x_plus_y_ao_m, x_plus_y_ao_n,
                                          x_minus_y_ao_m, x_minus_y_ao_n, fock_ao_rhs_1_m,
                                          fock_ao_rhs_2_m, fock_ao_rhs_1_n, fock_ao_rhs_2_n,
                                          fmat):
        """
        Calculates the one-particle and two-particle density matrix contributions to the
        omega multipliers.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The SCF tensors.
        :param x_plus_y_ao_m/n:
            The m/n component of the X+Y response vectors in AO basis.
        :param x_minus_y_ao:
            The m/n component of the X-Y response vectors in AO basis.
        :param fock_ao_rhs_1_m/n:
            The m/n component of the RHS fock/X+Y matrix in AO basis.
        :param fock_ao_rhs_2_m/n:
            The m/n component of the RHS fock/X-Y matrix in AO basis.
        :param fmat:
            The contributions from the two-particle and relaxed one-particle
            density matrix in AO basis.

        :return omega_1pdm_2pdm_contrib:
            The one-particle and two-particle density matric contributions to
            the omega multipliers.
        """

        # MO coefficients
        mo = scf_tensors['C_alpha']  # only alpha part
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()

        # overlap
        ovlp = scf_tensors['S']

        # calculate the ground state density matrices
        D_occ = np.matmul(mo_occ, mo_occ.T)
        D_vir = np.matmul(mo_vir, mo_vir.T)

        # contract Fock RHS matrices with response vectors
        Fp1_vv = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_1_m.T, x_plus_y_ao_n, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_1_n.T, x_plus_y_ao_m, ovlp.T]))

        Fm1_vv = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_2_m.T, x_minus_y_ao_n, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_2_n.T, x_minus_y_ao_m, ovlp.T]))

        Fp2_vv = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_1_m, x_plus_y_ao_n, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_1_n, x_plus_y_ao_m, ovlp.T]))

        Fm2_vv = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_2_m, x_minus_y_ao_n, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_2_n, x_minus_y_ao_m, ovlp.T]))

        Fp1_oo = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_1_m, x_plus_y_ao_n.T, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_1_n, x_plus_y_ao_m.T, ovlp.T]))

        Fm1_oo = 0.25 * (np.linalg.multi_dot([
            fock_ao_rhs_2_m, x_minus_y_ao_n.T, ovlp.T
        ]) + np.linalg.multi_dot(
            [fock_ao_rhs_2_n, x_minus_y_ao_m.T, ovlp.T]))

        Fp2_oo = 0.25 * (np.linalg.multi_dot([
            ovlp, x_plus_y_ao_n, fock_ao_rhs_1_m]).T
            + np.linalg.multi_dot([
            ovlp, x_plus_y_ao_m, fock_ao_rhs_1_n]).T)

        Fm2_oo = 0.25 * (np.linalg.multi_dot([
            ovlp, x_minus_y_ao_n, fock_ao_rhs_2_m]).T
            + np.linalg.multi_dot([
            ovlp, x_minus_y_ao_m, fock_ao_rhs_2_n]).T)
        # We see that:
        # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
        # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

        omega_1pdm_2pdm_contrib = -1.0 * (np.linalg.multi_dot([
            D_vir, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
        ]) + np.linalg.multi_dot([
            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
        ]) + np.linalg.multi_dot([
            D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir
        ]).T + np.linalg.multi_dot([
            D_occ, Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ
        ]) + np.linalg.multi_dot([D_occ, fmat, D_occ]))

        return omega_1pdm_2pdm_contrib

    def calculate_omega_gxc_contrib_real(self, scf_tensors, fock_gxc_ao_mn, D_occ):
        """
        Calculates the contribution to the real omega multipliers from the
        DFT E[3] g^xc term.

        :param scf_tensors:
            The tensors from the SCF calculation.
        :param fock_gxc_ao_mn:
            The mn component of the integrated g^xc Fock matrix
        :param D_occ:
            The occ/occ MO density matrix.

        :return omega_gxc_contrib:
            The E[3] g^xc contribution to the omega multipliers.
        """

        # degrees of freedom
        dof = len(self.vector_components)

        factor = -0.5
        omega_gxc_contrib = factor * np.linalg.multi_dot([
            D_occ, fock_gxc_ao_mn, D_occ])

        return omega_gxc_contrib

    def calculate_omega_gxc_contrib_complex(self, fock_gxc_ao_mn_list, D_occ):
        """
        Calculates the contribution to the complex omega multipliers from the
        DFT E[3] g^xc term.

        :param fock_gxc_ao_mn_list:
            List with the rere/imim/reim/imre parts of the mn component
            of the integrated g^xc Fock matrix
        :param D_occ:
            The occ/occ MO density matrix.

        :return omega_gxc_contrib:
            The E[3] g^xc contribution to the omega multipliers.
        """

        # unpack list with Fock matrices
        fock_gxc_ao_rere_mn = fock_gxc_ao_mn_list[0]
        fock_gxc_ao_imim_mn = fock_gxc_ao_mn_list[1]
        fock_gxc_ao_reim_mn = fock_gxc_ao_mn_list[2]
        fock_gxc_ao_imre_mn = fock_gxc_ao_mn_list[3]

        factor = -0.5
        omega_gxc_contrib = factor * (
                np.linalg.multi_dot([D_occ, fock_gxc_ao_rere_mn, D_occ])
                - np.linalg.multi_dot([D_occ, fock_gxc_ao_imim_mn, D_occ])
                + 1j * (np.linalg.multi_dot([D_occ, fock_gxc_ao_reim_mn, D_occ])
                + np.linalg.multi_dot([D_occ, fock_gxc_ao_imre_mn, D_occ]))
        )

        return omega_gxc_contrib
   
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

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

import time as tm
import math
import sys
from mpi4py import MPI
import numpy as np

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCMolecularGradient
from .veloxchemlib import T4CScreener
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver, NuclearPotentialGeom010Driver,
                           FockGeom1000Driver, ElectricDipoleMomentGeom100Driver)
from .veloxchemlib import mpi_master, hartree_in_wavenumber, denmat
from .veloxchemlib import partition_atoms, make_matrix, mat_t

from .polorbitalresponse import PolOrbitalResponse
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .molecule import Molecule
from .outputstream import OutputStream
from .matrices import Matrices
from .griddriver import GridDriver
from .inputparser import parse_input
from .profiler import Profiler
from .dftutils import get_default_grid_level
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, polgrad_sanity_check)


class PolarizabilityGradient:
    """
    Implements the dipole polarizability gradient.

    Instance variables
        - polgradient: The polarizability gradient
        - delta_h: the numerical perturbation
        - is_complex: Complex polarizability
        - damping: Damping factor for complex polarizability gradient
        - numerical: Numerical differentiation
        - do_four_point: Four-point numerical differentiation
        - frequencies: The frequencies
        - vector_components: Cartesian components of the tensor
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes polarizability gradient driver to default setup.

        :param scf_drv:
            The SCF driver
        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        # timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        self._scf_drv = scf_drv

        self.method_dict = {}
        self.orbrsp_dict = {}

        self.polgradient = None
        self.delta_h = 0.001

        self.is_complex = False
        self.grad_dt = np.dtype('float64')  # data type for pol. gradient (real/complex)
        self.damping = 1000.0 / hartree_in_wavenumber()

        self.numerical = False
        self.do_four_point = False
        self.do_print_polgrad = False

        self._dft = False
        self.grid_level = None
        self.xcfun = None
        self._xcfun_ldstaging = 1024

        self.flag = 'Polarizability Gradient Driver'
        self.frequencies = (0,)
        self.vector_components = 'xyz'

        self._input_keywords = {
            'polarizabilitygradient': {
                'vector_components': ('str_lower', 'Cartesian components of operator'),
                'frequencies': ('seq_range', 'frequencies'),
                'numerical': ('bool', 'do numerical integration'),
                'do_four_point': ('bool', 'do four-point numerical integration'),
                'delta_h': ('float', 'the displacement for finite difference'),
                'is_complex': ('bool', 'whether the polarizability is complex'),
                'damping': ('float', 'damping parameter for complex numerical (a.u.)'),
                'do_print_polgrad': ('bool', 'whether to print the pol. gradient'),
                'timing': ('bool', 'print timing information'),
                'profiling': ('bool', 'print profiling information'),
                'memory_profiling': ('bool', 'print memory usage'),
                'memory_tracing': ('bool', 'trace memory allocation'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
            }
        }

    def update_settings(self, grad_dict, orbrsp_dict=None, method_dict=None):
        """
        Updates response and method settings in polarizability gradient
        computation driver.

        :param grad_dict:
            The input dictionary of gradient input.
        :param orbrsp_dict:
            The dictionary of orbital response (CPHF) input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}
        if orbrsp_dict is None:
            orbrsp_dict = {}

        grad_keywords = {
            key: val[0] for key, val in
            self._input_keywords['polarizabilitygradient'].items()
        }

        parse_input(self, grad_keywords, grad_dict)

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

        self.frequencies = list(self.frequencies)

        self.method_dict = dict(method_dict)
        self.orbrsp_dict = dict(orbrsp_dict)

    def compute(self, molecule, basis, scf_tensors, lr_results=None):
        """
        Calls the correct function to perform the calculation of
        the polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.
        """

        if self.rank == mpi_master():
            self.print_header()

        # set data type of pol. gradient for use in compute_analytical()
        if self.is_complex:
            self.grad_dt = np.dtype('complex128')

        # sanity checks
        molecule_sanity_check(molecule)
        scf_results_sanity_check(self, self._scf_drv.scf_tensors)
        dft_sanity_check(self, 'compute')

        start_time = tm.time()

        if self.numerical:
            # compute
            self.polgradient = self.compute_numerical(molecule, basis, self._scf_drv)
        else:
            # sanity checks linear response input
            if (self.rank == mpi_master()) and (lr_results is None):
                error_message = 'PolarizabilityGradient missing input: LR results'
                error_message += 'for analytical gradient'
                raise ValueError(error_message)
            if self.rank == mpi_master():
                polgrad_sanity_check(self, self.flag, lr_results)
                self._check_real_or_complex_input(lr_results)
            # compute
            self.polgradient = self.compute_analytical(molecule, basis,
                                                       scf_tensors, lr_results)

        if self.rank == mpi_master():

            if self.do_print_polgrad:
                self.print_polarizability_gradient(molecule)

            valstr = '*** Time spent in polarizability gradient driver: '
            valstr += f'{(tm.time() - start_time):.6f} sec ***'
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        return self.polgradient

    def compute_analytical(self, molecule, basis, scf_tensors, lr_results):
        """
        Performs calculation of the analytical polarizability gradient
        using the veloxchem integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.
        """

        # profiling
        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # get orbital response results
        orbrsp_results = self.compute_orbital_response(
            molecule, basis, scf_tensors, lr_results)

        # dictionary for polarizability gradient
        polgrad_results = {}

        # timings
        loop_start_time = tm.time()

        # partition atoms for parallelisation
        natm = molecule.number_of_atoms()

        local_atoms = partition_atoms(natm, self.rank, self.nodes)

        # operator components and permutation pairs
        dof = len(self.vector_components)
        xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]
        dof_red = len(xy_pairs)

        # number of atomic orbitals
        nao = basis.get_dimensions_of_basis()

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']  # only alpha part

            # MO coefficients
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            # number of input frequencies
            n_freqs = len(self.frequencies)

            # if Lagr. multipliers not available as dist. array,
            # read it from dict into variable on master
            if 'dist_cphf_ov' not in orbrsp_results.keys():
                all_cphf_red = orbrsp_results['cphf_ov']
                if self.is_complex:
                    all_cphf_red = all_cphf_red.reshape(
                        n_freqs, 2 * dof_red, nocc * nvir)
                else:
                    all_cphf_red = all_cphf_red.reshape(n_freqs, dof_red, nocc * nvir)
        else:
            nocc = None
            nvir = None

        nocc, nvir = self.comm.bcast((nocc, nvir), root=mpi_master())

        for f, w in enumerate(self.frequencies):

            profiler.set_timing_key(f"Polgrad w = {w}")

            full_vec = [
                self.get_full_solution_vector(lr_results['solutions'][x, w])
                for x in self.vector_components
            ]

            if 'dist_cphf_ov' in orbrsp_results.keys():
                # get lambda multipliers from distributed arrays
                cphf_ov = self.get_lambda_response_vector(
                    molecule, scf_tensors, orbrsp_results['dist_cphf_ov'], f)

            else:
                # TODO get conj. gradient solver to also return dist. array
                # get lambda multipliers from array in orbrsp dict
                cphf_ov_red = all_cphf_red[f]
                cphf_ov = np.zeros((dof, dof, nocc * nvir),
                                   dtype=self.grad_dt)

                if self.is_complex:
                    tmp_cphf_ov = cphf_ov_red[:dof_red] + 1j * cphf_ov_red[dof_red:]

                    for idx, xy in enumerate(xy_pairs):
                        x = xy[0]
                        y = xy[1]

                        cphf_ov[x, y] = tmp_cphf_ov[idx]

                        if y != x:
                            cphf_ov[y, x] += cphf_ov[x, y]
                else:
                    for idx, xy in enumerate(xy_pairs):
                        x = xy[0]
                        y = xy[1]

                        cphf_ov[x, y] = cphf_ov_red[idx]

                        if y != x:
                            cphf_ov[y, x] += cphf_ov[x, y]

            omega_ao = self.get_omega_response_vector(
                basis, orbrsp_results['dist_omega_ao'], f
            )

            if self.rank == mpi_master():
                info_msg = f'Building gradient for frequency = {w:4.3f}'
                self.ostream.print_info(info_msg)
                self.ostream.print_blank()
                self.ostream.flush()

                # Note: polorbitalresponse uses r instead of mu for dipole operator
                # for idx in range(len(full_vec)):
                for idx, vec in enumerate(full_vec):
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
                x_plus_y_mo = exc_vec + deexc_vec
                x_minus_y_mo = exc_vec - deexc_vec

                del exc_vec, deexc_vec

                # transform to AO basis: mi,xia,na->xmn
                x_plus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_plus_y_mo[x], mo_vir.T])
                    for x in range(dof)
                ])
                x_minus_y_ao = np.array([
                    np.linalg.multi_dot([mo_occ, x_minus_y_mo[x], mo_vir.T])
                    for x in range(dof)
                ])

                gs_dm = scf_tensors['D_alpha']  # only alpha part

                cphf_ov = cphf_ov.reshape((dof**2, nocc, nvir))

                lambda_ao = np.array([
                    np.linalg.multi_dot([mo_occ, cphf_ov[xy], mo_vir.T])
                    for xy in range(dof**2)
                ])

                lambda_ao = lambda_ao.reshape((dof, dof, nao, nao))
                lambda_ao += lambda_ao.transpose(0, 1, 3, 2)

                # calculate symmetrized unrelaxed density matrix
                unrel_dm_ao = self.calculate_unrel_dm(molecule, scf_tensors,
                                                      x_plus_y_mo, x_minus_y_mo)

                # calculate relaxed density matrix
                rel_dm_ao = unrel_dm_ao + lambda_ao

                del x_plus_y_mo, x_minus_y_mo, cphf_ov, unrel_dm_ao, lambda_ao
            else:
                gs_dm = None
                x_plus_y_ao = None
                x_minus_y_ao = None
                omega_ao = None
                rel_dm_ao = None

            gs_dm = self.comm.bcast(gs_dm, root=mpi_master())
            x_plus_y_ao = self.comm.bcast(x_plus_y_ao, root=mpi_master())
            x_minus_y_ao = self.comm.bcast(x_minus_y_ao, root=mpi_master())
            omega_ao = self.comm.bcast(omega_ao, root=mpi_master())
            rel_dm_ao = self.comm.bcast(rel_dm_ao, root=mpi_master())

            # initiate polarizability gradient variable with data type set in init()
            pol_gradient = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

            # kinetic energy gradient driver
            kin_grad_drv = KineticEnergyGeom100Driver()

            # nuclear potential gradienr drivers
            npot_grad_100_drv = NuclearPotentialGeom100Driver()
            npot_grad_010_drv = NuclearPotentialGeom010Driver()

            profiler.start_timer("1-elec contribs")

            # loop over atoms and contract integral derivatives with density matrices
            for iatom in local_atoms:
                # core Hamiltonian contribution to gradients
                gmats_kin = kin_grad_drv.compute(molecule, basis, iatom)
                gmats_npot_100 = npot_grad_100_drv.compute(molecule, basis, iatom)
                gmats_npot_010 = npot_grad_010_drv.compute(molecule, basis, iatom)

                for icoord, label in enumerate(['X', 'Y', 'Z']):
                    gmat_kin = gmats_kin.matrix_to_numpy(label)
                    gmat_npot_100 = gmats_npot_100.matrix_to_numpy(label)
                    gmat_npot_010 = gmats_npot_010.matrix_to_numpy(label)
                    gmat_hcore = gmat_kin + gmat_kin.T
                    gmat_hcore -= gmat_npot_100 + gmat_npot_100.T + gmat_npot_010

                    # loop over operator components
                    for x, y in xy_pairs:
                        pol_gradient[x, y, iatom, icoord] += (
                            np.linalg.multi_dot([  # xymn,amn->xya
                                2.0 * rel_dm_ao[x, y].reshape(nao**2),
                                gmat_hcore.reshape(nao**2)])
                        )

                gmats_kin = Matrices()
                gmats_npot_100 = Matrices()
                gmats_npot_010 = Matrices()
                gmat_hcore = Matrices()

                # overlap contribution to gradient
                ovlp_grad_drv = OverlapGeom100Driver()
                gmats_ovlp = ovlp_grad_drv.compute(molecule, basis, iatom)

                for icoord, label in enumerate(['X', 'Y', 'Z']):
                    gmat_ovlp = gmats_ovlp.matrix_to_numpy(label)
                    gmat_ovlp += gmat_ovlp.T
                    # loop over operator components
                    for x, y in xy_pairs:
                        pol_gradient[x, y, iatom, icoord] += (
                            1.0 * np.linalg.multi_dot([  # xymn,amn->xya
                                2.0 * omega_ao[x, y],
                                gmat_ovlp.reshape(nao**2)])
                        )

                gmats_ovlp = Matrices()

                # dipole contribution to gradient
                dip_grad_drv = ElectricDipoleMomentGeom100Driver()
                gmats_dip = dip_grad_drv.compute(molecule, basis,
                                                 [0.0, 0.0, 0.0], iatom)

                # the keys of the dipole gmat
                gmats_dip_components = (['X_X', 'X_Y', 'X_Z', 'Y_X']
                                        + ['Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'])

                # dictionary to convert from string idx to integer idx
                comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

                d_dipole = np.zeros((dof, 3, nao, nao))

                for label in gmats_dip_components:
                    gmat_dip = gmats_dip.matrix_to_numpy(label)
                    gmat_dip += gmat_dip.T

                    icoord = comp_to_idx[label[0]]  # atom coordinate component
                    icomp = comp_to_idx[label[-1]]  # dipole operator component

                    # reorder indices to first is operator comp, second is coord
                    d_dipole[icomp, icoord] += gmat_dip

                for icoord in range(3):
                    # loop over operator components
                    for x, y in xy_pairs:
                        pol_gradient[x, y, iatom, icoord] += (
                            - 2.0 * np.linalg.multi_dot([  # xmn,yamn->xya
                                x_minus_y_ao[x].reshape(nao**2),
                                d_dipole[y, icoord].reshape(nao**2)
                            ]) - 2.0 * np.linalg.multi_dot([  # xmn,yamn->yxa
                                d_dipole[x, icoord].reshape(nao**2),
                                x_minus_y_ao[y].reshape(nao**2)
                            ]))

                gmats_dip = Matrices()

            profiler.stop_timer("1-elec contribs")
            profiler.check_memory_usage(f"1-elec grad w = {w}")

            # ERI contribution
            profiler.start_timer("ERI contrib")
            eri_contrib = self.compute_eri_contrib(molecule, basis, gs_dm,
                                                   rel_dm_ao, x_plus_y_ao, x_minus_y_ao,
                                                   local_atoms)
            profiler.stop_timer("ERI contrib")
            profiler.check_memory_usage(f"ERI grad w = {w}")

            pol_gradient += eri_contrib
            pol_gradient = self.comm.reduce(pol_gradient, root=mpi_master())

            del eri_contrib

            if self.rank == mpi_master():
                for x in range(dof):
                    for y in range(x + 1, dof):
                        pol_gradient[y, x] += pol_gradient[x, y]

            if self._dft:
                profiler.start_timer("XC contrib")

                # compute the XC contribution
                polgrad_xc_contrib = self.compute_polgrad_xc_contrib(
                    molecule, basis, gs_dm, rel_dm_ao, x_minus_y_ao, profiler)

                # add contribution to the SCF polarizability gradient
                if self.rank == mpi_master():
                    pol_gradient += polgrad_xc_contrib

                profiler.stop_timer("XC contrib")
                profiler.check_memory_usage(f'XC contrib w = {w}')

                del polgrad_xc_contrib

            if self.rank == mpi_master():
                polgrad_results[w] = pol_gradient.reshape(dof, dof, 3 * natm)

        profiler.print_timing(self.ostream)
        profiler.print_memory_usage(self.ostream)

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing the analytical gradient for '
            valstr += f'{n_freqs:d} frequencies: '
            valstr += f'{(tm.time() - loop_start_time):.6f} sec **'
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        # NOTE distributed?
        return polgrad_results

    def compute_eri_contrib(self, molecule, basis, gs_dm, rel_dm_ao,
                            x_plus_y_ao, x_minus_y_ao, local_atoms):
        """
        Directs the computation of the contribution from ERI derivative integrals
        to the polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param gs_dm:
            The ground state density.
        :param rel_dm_ao:
            The relaxed density matrix.
        :param x_plus_y_ao:
            The X+Y response vectors.
        :param x_minus_y_ao:
            The X-Y response vectors.
        :param local_atoms:
            The atom partition for the MPI node.

        :return eri_deriv_contrib:
            The ERI derivative integral contribution to the polarizability gradient.
        """

        if self.is_complex:
            eri_deriv_contrib = self.compute_eri_contrib_complex(molecule, basis, gs_dm,
                                                                 rel_dm_ao, x_plus_y_ao,
                                                                 x_minus_y_ao,
                                                                 local_atoms)
        else:
            eri_deriv_contrib = self.compute_eri_contrib_real(molecule, basis, gs_dm,
                                                              rel_dm_ao, x_plus_y_ao,
                                                              x_minus_y_ao, local_atoms)

        return eri_deriv_contrib

    def compute_eri_contrib_real(self, molecule, basis, gs_dm, rel_dm_ao,
                                 x_plus_y_ao, x_minus_y_ao, local_atoms):
        """
        Computes the contribution from ERI derivative integrals
        to the real polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param gs_dm:
            The ground state density.
        :param rel_dm_ao:
            The relaxed density matrix.
        :param x_plus_y_ao:
            The X+Y response vectors.
        :param x_minus_y_ao:
            The X-Y response vectors.
        :param local_atoms:
            The atom partition for the MPI node.

        :return eri_deriv_contrib:
            The ERI derivative integral contribution to the polarizability gradient.
        """

        fock_type, exchange_scaling_factor = self.get_fock_type_and_x_frac()

        # scaling of ERI gradient for non-hybrid functionals
        factor = 2.0 if fock_type == 'j' else 1.0

        # ERI threshold
        thresh_int = int(-math.log10(self._scf_drv.eri_thresh))

        # Fock gradient driver
        fock_grad_drv = FockGeom1000Driver()

        # screening
        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        # contraction with ground state density
        den_mat_for_fock_gs = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock_gs.set_values(gs_dm)

        # number of atomic orbitals
        natm = molecule.number_of_atoms()

        # degrees of freedom
        dof = len(self.vector_components)

        # operator component combinations
        xy_pairs = [(x,y) for x in range(dof) for y in range(x,dof)]

        eri_deriv_contrib = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

        # ERI gradient
        for iatom in local_atoms:
            # screening
            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            # contraction with density matrices
            for x, y in xy_pairs:
                # relaxed DM
                den_mat_for_fock_rel = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel.set_values(2.0 * rel_dm_ao[x,y])
                # (X+Y)_x
                den_mat_for_fock_xpy_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x.set_values(x_plus_y_ao[x])
                # (X+Y)_y
                den_mat_for_fock_xpy_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y.set_values(x_plus_y_ao[y])
                # (X+Y)_x - (X+Y)_x
                den_mat_for_fock_xpy_m_xpyT_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x.set_values(x_plus_y_ao[x]
                                                         - x_plus_y_ao[x].T)
                # (X+Y)_y - (X+Y)_y
                den_mat_for_fock_xpy_m_xpyT_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y.set_values(x_plus_y_ao[y]
                                                         - x_plus_y_ao[y].T)
                # (X-Y)_x
                den_mat_for_fock_xmy_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x.set_values(x_minus_y_ao[x])
                # (X-Y)_y
                den_mat_for_fock_xmy_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y.set_values(x_minus_y_ao[y])
                # (X-Y)_x + (X-Y)_x
                den_mat_for_fock_xmy_p_xmyT_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x.set_values(x_minus_y_ao[x]
                                                         + x_minus_y_ao[x].T)
                # (X-Y)_y + (X-Y)_y
                den_mat_for_fock_xmy_p_xmyT_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y.set_values(x_minus_y_ao[y]
                                                         + x_minus_y_ao[y].T)
                # contraction of integrals and DMs
                erigrad_rel = fock_grad_drv.compute(basis, screener_atom, screener,
                                                    den_mat_for_fock_gs,
                                                    den_mat_for_fock_rel,
                                                    iatom, fock_type,
                                                    exchange_scaling_factor,
                                                    0.0, thresh_int)
                erigrad_xpy_xy = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_xpy_x,
                                                       den_mat_for_fock_xpy_m_xpyT_y,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                erigrad_xpy_yx = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_xpy_m_xpyT_x,
                                                       den_mat_for_fock_xpy_y,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                erigrad_xmy_xy = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_xmy_x,
                                                       den_mat_for_fock_xmy_p_xmyT_y,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                erigrad_xmy_yx = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_xmy_p_xmyT_x,
                                                       den_mat_for_fock_xmy_y,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                eri_deriv_contrib[x, y, iatom] += np.array(erigrad_rel)
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xpy_xy)
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xpy_yx)
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xmy_xy)
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xmy_yx)
                eri_deriv_contrib[x, y, iatom] *= factor

        return eri_deriv_contrib

    def compute_eri_contrib_complex(self, molecule, basis, gs_dm, rel_dm_ao,
                                    x_plus_y_ao, x_minus_y_ao, local_atoms):
        """
        Computes the contribution from ERI derivative integrals
        to the complex polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param gs_dm:
            The ground state density.
        :param rel_dm_ao:
            The relaxed density matrix.
        :param x_plus_y_ao:
            The X+Y response vectors.
        :param x_minus_y_ao:
            The X-Y response vectors.
        :param local_atoms:
            The atom partition for the MPI node.

        :return eri_deriv_contrib:
            The ERI derivative integral contribution to the polarizability gradient.
        """

        # Fock type and sclaing of exact exhange
        fock_type, exchange_scaling_factor = self.get_fock_type_and_x_frac()

        # scaling of ERI gradient for non-hybrid functionals
        factor = 2.0 if fock_type == 'j' else 1.0

        # ERI threshold
        thresh_int = int(-math.log10(self._scf_drv.eri_thresh))

        # Fock gradient driver
        fock_grad_drv = FockGeom1000Driver()

        # screening
        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        # contraction with ground state density
        den_mat_for_fock_gs = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock_gs.set_values(gs_dm)

        # number of atomic orbitals
        natm = molecule.number_of_atoms()

        # degrees of freedom
        dof = len(self.vector_components)

        # operator component combinations
        xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]

        eri_deriv_contrib = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

        for iatom in local_atoms:
            # screening
            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            # contraction with density matrices
            for x, y in xy_pairs:
                # relaxed DM
                den_mat_for_fock_rel_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel_real.set_values(2.0 * rel_dm_ao[x,y].real)
                den_mat_for_fock_rel_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel_imag.set_values(2.0 * rel_dm_ao[x,y].imag)
                # (X+Y)_x
                den_mat_for_fock_xpy_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x_real.set_values(x_plus_y_ao[x].real)
                den_mat_for_fock_xpy_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x_imag.set_values(x_plus_y_ao[x].imag)
                # (X+Y)_y
                den_mat_for_fock_xpy_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y_real.set_values(x_plus_y_ao[y].real)
                den_mat_for_fock_xpy_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y_imag.set_values(x_plus_y_ao[y].imag)
                # (X+Y)_x - (X+Y)_x
                den_mat_for_fock_xpy_m_xpyT_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x_real.set_values((x_plus_y_ao[x]
                                                               - x_plus_y_ao[x].T).real)
                den_mat_for_fock_xpy_m_xpyT_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x_imag.set_values((x_plus_y_ao[x]
                                                               - x_plus_y_ao[x].T).imag)
                # (X+Y)_y - (X+Y)_y
                den_mat_for_fock_xpy_m_xpyT_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y_real.set_values((x_plus_y_ao[y]
                                                               - x_plus_y_ao[y].T).real)
                den_mat_for_fock_xpy_m_xpyT_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y_imag.set_values((x_plus_y_ao[y]
                                                               - x_plus_y_ao[y].T).imag)
                # (X-Y)_x
                den_mat_for_fock_xmy_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x_real.set_values(x_minus_y_ao[x].real)
                den_mat_for_fock_xmy_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x_imag.set_values(x_minus_y_ao[x].imag)
                # (X-Y)_y
                den_mat_for_fock_xmy_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y_real.set_values(x_minus_y_ao[y].real)
                den_mat_for_fock_xmy_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y_imag.set_values(x_minus_y_ao[y].imag)
                # (X-Y)_x + (X-Y)_x
                den_mat_for_fock_xmy_p_xmyT_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x_real.set_values(
                    (x_minus_y_ao[x] + x_minus_y_ao[x].T).real)
                den_mat_for_fock_xmy_p_xmyT_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x_imag.set_values((x_minus_y_ao[x]
                                                               + x_minus_y_ao[x].T).imag)
                # (X-Y)_y + (X-Y)_y
                den_mat_for_fock_xmy_p_xmyT_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y_real.set_values((x_minus_y_ao[y]
                                                               + x_minus_y_ao[y].T).real)
                den_mat_for_fock_xmy_p_xmyT_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y_imag.set_values((x_minus_y_ao[y]
                                                               + x_minus_y_ao[y].T).imag)
                # contraction of integrals and DMs
                # Re
                erigrad_rel_re = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_gs,
                                                       den_mat_for_fock_rel_real,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                # Im
                erigrad_rel_im = fock_grad_drv.compute(basis, screener_atom, screener,
                                                       den_mat_for_fock_gs,
                                                       den_mat_for_fock_rel_imag,
                                                       iatom, fock_type,
                                                       exchange_scaling_factor,
                                                       0.0, thresh_int)
                # ReRe
                erigrad_xpy_xy_rere = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_x_real,
                                                            den_mat_for_fock_xpy_m_xpyT_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImIm
                erigrad_xpy_xy_imim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_x_imag,
                                                            den_mat_for_fock_xpy_m_xpyT_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReIm
                erigrad_xpy_xy_reim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_x_real,
                                                            den_mat_for_fock_xpy_m_xpyT_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImRe
                erigrad_xpy_xy_imre = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_x_imag,
                                                            den_mat_for_fock_xpy_m_xpyT_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReRe
                erigrad_xpy_yx_rere = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_m_xpyT_x_real,
                                                            den_mat_for_fock_xpy_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImIm
                erigrad_xpy_yx_imim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_m_xpyT_x_imag,
                                                            den_mat_for_fock_xpy_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReIm
                erigrad_xpy_yx_reim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_m_xpyT_x_real,
                                                            den_mat_for_fock_xpy_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImRe
                erigrad_xpy_yx_imre = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xpy_m_xpyT_x_imag,
                                                            den_mat_for_fock_xpy_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReRe
                erigrad_xmy_xy_rere = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_x_real,
                                                            den_mat_for_fock_xmy_p_xmyT_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImIm
                erigrad_xmy_xy_imim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_x_imag,
                                                            den_mat_for_fock_xmy_p_xmyT_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReIm
                erigrad_xmy_xy_reim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_x_real,
                                                            den_mat_for_fock_xmy_p_xmyT_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImRe
                erigrad_xmy_xy_imre = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_x_imag,
                                                            den_mat_for_fock_xmy_p_xmyT_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReRe
                erigrad_xmy_yx_rere = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_p_xmyT_x_real,
                                                            den_mat_for_fock_xmy_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImIm
                erigrad_xmy_yx_imim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_p_xmyT_x_imag,
                                                            den_mat_for_fock_xmy_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ReIm
                erigrad_xmy_yx_reim = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_p_xmyT_x_real,
                                                            den_mat_for_fock_xmy_y_imag,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)
                # ImRe
                erigrad_xmy_yx_imre = fock_grad_drv.compute(basis, screener_atom,
                                                            screener,
                                                            den_mat_for_fock_xmy_p_xmyT_x_imag,
                                                            den_mat_for_fock_xmy_y_real,
                                                            iatom, fock_type,
                                                            exchange_scaling_factor,
                                                            0.0, thresh_int)

                # Real fmat contribution to gradient
                erigrad_real = np.array(erigrad_rel_re)  # real DM
                erigrad_real += 0.5 * np.array(erigrad_xpy_xy_rere)  # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xpy_xy_imim)  # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xpy_yx_rere)  # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xpy_yx_imim)  # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xmy_xy_rere)  # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xmy_xy_imim)  # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xmy_yx_rere)  # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xmy_yx_imim)  # ImIm
                erigrad_real *= factor

                # Imaginary fmat contribution to gradient
                erigrad_imag = np.array(erigrad_rel_im) # imag DM
                erigrad_imag += 0.5 * np.array(erigrad_xpy_xy_reim)  # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xpy_xy_imre)  # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xpy_yx_reim)  # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xpy_yx_imre)  # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xmy_xy_reim)  # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xmy_xy_imre)  # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xmy_yx_reim)  # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xmy_yx_imre)  # ImRe
                erigrad_imag *= factor

                # add to complex variable
                eri_deriv_contrib[x, y, iatom] += erigrad_real + 1j * erigrad_imag

        return eri_deriv_contrib

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

    def get_lambda_response_vector(self, molecule, scf_tensors, lambda_list, fdx):
        """
        Gets the full lambda multipliers vector from distributed array.

        :param lambda_list:
            The list with distributed arrays.
        :param fdx:
            The frequency index.

        :return:
            The orbital response lambda vector.
        """

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']  # only alpha part
            nocc = molecule.number_of_alpha_electrons()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            del mo, mo_vir

            dof = len(self.vector_components)
            xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]
            dof_red = len(xy_pairs)

            lambda_ov = np.zeros((dof, dof, nocc * nvir),
                                 dtype=self.grad_dt)
        else:
            dof_red = None
            xy_pairs = None
            lambda_ov = None

        dof_red, xy_pairs = self.comm.bcast((dof_red, xy_pairs), root=mpi_master())

        if self.is_complex:

            for idx, xy in enumerate(xy_pairs):
                tmp_lambda_re = lambda_list[
                    2 * dof_red * fdx + idx].get_full_vector()
                tmp_lambda_im = lambda_list[
                    2 * dof_red * fdx + dof_red + idx].get_full_vector()

                if self.rank == mpi_master():
                    x = xy[0]
                    y = xy[1]

                    lambda_ov[x, y] += tmp_lambda_re + 1j * tmp_lambda_im

                    if y != x:
                        lambda_ov[y, x] += lambda_ov[x, y]
        else:
            for idx, xy in enumerate(xy_pairs):
                tmp_lambda_ov = lambda_list[dof_red * fdx + idx].get_full_vector()

                if self.rank == mpi_master():
                    x = xy[0]
                    y = xy[1]

                    lambda_ov[x, y] += tmp_lambda_ov

                    if y != x:
                        lambda_ov[y, x] += lambda_ov[x, y]

        return lambda_ov

    def get_omega_response_vector(self, basis, omega_list, fdx):
        """
        Gets the full omega multipliers vector from distributed array.

        :param omega_list:
            The list with distributed arrays.
        :param fdx:
            The frequency index.

        :return:
            The orbital response omega vector.
        """

        if self.rank == mpi_master():
            dof = len(self.vector_components)
            xy_pairs = [(x, y) for x in range(dof) for y in range(x, dof)]
            dof_red = len(xy_pairs)
            nao = basis.get_dimensions_of_basis()

            omega_ao = np.zeros((dof, dof, nao * nao), dtype=self.grad_dt)
        else:
            dof_red = None
            xy_pairs = None
            omega_ao = None

        dof_red, xy_pairs = self.comm.bcast((dof_red, xy_pairs), root=mpi_master())

        if self.is_complex:
            for idx, xy in enumerate(xy_pairs):
                tmp_omega_re = omega_list[
                    2 * dof_red * fdx + idx].get_full_vector()
                tmp_omega_im = omega_list[
                    2 * dof_red * fdx + dof_red + idx].get_full_vector()

                if self.rank == mpi_master():
                    x = xy[0]
                    y = xy[1]

                    omega_ao[x, y] += tmp_omega_re + 1j * tmp_omega_im

                    if y != x:
                        omega_ao[y, x] += omega_ao[x, y]
        else:
            for idx, xy in enumerate(xy_pairs):
                tmp_omega_ao = omega_list[dof_red * fdx + idx].get_full_vector()

                if self.rank == mpi_master():
                    x = xy[0]
                    y = xy[1]
                    omega_ao[x, y] += tmp_omega_ao

                    if y != x:
                        omega_ao[y, x] += omega_ao[x, y]

        return omega_ao

    def compute_numerical(self, molecule, ao_basis, scf_drv):
        """
        Performs calculation of numerical nuclear gradient
        of the electric dipole polarizability.

        Moved from tddftgradientdriver.py 18/10/2023

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        """

        # set linear response driver
        lr_drv = self.set_lr_driver()

        # mute the outputstreams
        scf_drv.ostream.mute()
        lr_drv.ostream.mute()

        # construct polarizability gradient
        num_polgradient = self.construct_numerical_gradient(molecule, ao_basis, scf_drv, lr_drv)

        # unmute the output streams
        scf_drv.ostream.unmute()
        lr_drv.ostream.unmute()

        return num_polgradient

    def set_lr_driver(self):
        """
        Sets the linear response driver for numerical calculations.

        :return lr_drv:
            The linear response: LR or CPP
        """

        if self.is_complex:
            lr_drv = ComplexResponse(self.comm, self.ostream)
            lr_drv.frequencies = self.frequencies
            lr_drv.damping = self.damping
        else:
            lr_drv = LinearResponseSolver(self.comm, self.ostream)
            lr_drv.frequencies = self.frequencies

        return lr_drv

    def compute_orbital_response(self, molecule, basis, scf_tensors, lr_results):
        """
        Directs the calculation of the orbital response.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_tensors:
            The SCF tensors.
        :param lr_results:
            The results of the linear response calculation.
        """

        # timings
        orbrsp_start_time = tm.time()

        # setup orbital response driver
        orbrsp_drv = PolOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)

        # enforce the frequencies to be the same as in current driver
        orbrsp_drv.frequencies = self.frequencies
        orbrsp_drv.is_complex = self.is_complex

        # compute orbital response
        orbrsp_drv.compute(molecule, basis, scf_tensors, lr_results)
        orbrsp_drv.compute_omega(molecule, basis, scf_tensors, lr_results)

        if self.rank == mpi_master():
            valstr = f'** Time spent on orbital response for {len(self.frequencies)} frequencies: '
            valstr += f'{(tm.time() - orbrsp_start_time):.6f} sec **'
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        return orbrsp_drv.cphf_results

    def compute_polgrad_xc_contrib(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                   x_minus_y_ao, profiler):
        """
        Directs the calculation of the exchange-correlation contribution to the DFT
        polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param gs_dm:
            The ground state density matrix.
        :param rel_dm_ao:
            The relaxed density matric in AO basis.
        :param x_minus_y_ao:
            The X-Y response vector.

        :return xc_contrib:
            The XC-contribution to the polarizability gradient
        """

        xcfun_label = self.xcfun.get_func_label()

        if self.is_complex:
            xc_contrib = self.compute_polgrad_xc_contrib_complex(
                molecule, ao_basis, gs_dm, rel_dm_ao, x_minus_y_ao, xcfun_label, profiler)
        else:
            xc_contrib = self.compute_polgrad_xc_contrib_real(
                molecule, ao_basis, gs_dm, rel_dm_ao, x_minus_y_ao, xcfun_label, profiler)

        return xc_contrib

    def compute_polgrad_xc_contrib_real(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                        x_minus_y_ao, xcfun_label, profiler):
        """
        Calculates the exchange-correlation contribution to the real
        polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param gs_dm:
            The ground state density matrix.
        :param rel_dm_ao:
            The relaxed density matric in AO basis.
        :param x_minus_y_ao:
            The X-Y response vector.
        :param xcfun_label:
            The label for the XC functional

        :return xc_pol_gradient:
            The XC-contribution to the polarizability gradient
        """

        if self.rank == mpi_master():
            natm = molecule.number_of_atoms()
            dof = len(self.vector_components)
            xc_pol_gradient = np.zeros((dof, dof, natm, 3))
        else:
            dof = None
            xc_pol_gradient = None

        dof = self.comm.bcast(dof, root=mpi_master())
        xc_pol_gradient = self.comm.bcast(xc_pol_gradient, root=mpi_master())

        for m in range(dof):
            for n in range(m, dof):
                if self.rank == mpi_master():
                    rhow_dm = 1.0 * rel_dm_ao[m, n]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                    del rhow_dm

                    # symmetrize
                    x_minus_y_sym_m = np.sqrt(2) * 0.5 * (x_minus_y_ao[m] +
                                                          x_minus_y_ao[m].T)
                    x_minus_y_sym_n = np.sqrt(2) * 0.5 * (x_minus_y_ao[n] +
                                                          x_minus_y_ao[n].T)

                else:
                    rhow_dm_sym = None
                    x_minus_y_sym_m = None
                    x_minus_y_sym_n = None

                rhow_dm_sym = self.comm.bcast(rhow_dm_sym, root=mpi_master())
                x_minus_y_sym_m = self.comm.bcast(x_minus_y_sym_m, root=mpi_master())
                x_minus_y_sym_n = self.comm.bcast(x_minus_y_sym_n, root=mpi_master())

                polgrad_xcgrad = self.calculate_xc_mn_contrib_real(
                    molecule, ao_basis, [rhow_dm_sym],
                    [x_minus_y_sym_m], [x_minus_y_sym_n],
                    [gs_dm], xcfun_label, profiler)

                #if self.rank == mpi_master():
                xc_pol_gradient[m, n] += polgrad_xcgrad

                if n != m:
                    xc_pol_gradient[n, m] = xc_pol_gradient[m, n]

        xc_pol_gradient = self.comm.reduce(xc_pol_gradient, root=mpi_master())

        return xc_pol_gradient

    def compute_polgrad_xc_contrib_complex(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                           x_minus_y_ao, xcfun_label, profiler):
        """
        Calculates the exchange-correlation contribution to the complex
        polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param gs_dm:
            The ground state density matrix.
        :param rel_dm_ao:
            The relaxed density matric in AO basis.
        :param x_minus_y_ao:
            The X-Y response vector.
        :param xcfun_label:
            The label for the XC functional

        :return xc_pol_gradient:
            The XC-contribution to the polarizability gradient
        """

        if self.rank == mpi_master():
            natm = molecule.number_of_atoms()
            dof = len(self.vector_components)
            xc_pol_gradient = np.zeros((dof, dof, natm, 3),
                                       dtype=np.dtype('complex128'))
        else:
            dof = None
            xc_pol_gradient = None

        dof = self.comm.bcast(dof, root=mpi_master())
        xc_pol_gradient = self.comm.bcast(xc_pol_gradient, root=mpi_master())

        # FIXME change "m,n" to "x,y" for consistency
        for m in range(dof):
            for n in range(m, dof):
                if self.rank == mpi_master():
                    rhow_dm = 1.0 * rel_dm_ao[m, n]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)

                    rhow_dm_sym_list_real = [np.array(rhow_dm_sym.real)]
                    rhow_dm_sym_list_imag = [np.array(rhow_dm_sym.imag)]

                    # symmetrize
                    x_minus_y_sym_m = np.sqrt(2) * 0.5 * (x_minus_y_ao[m] +
                                                          x_minus_y_ao[m].T)
                    x_minus_y_sym_m_list_real = [np.array(x_minus_y_sym_m.real)]
                    x_minus_y_sym_m_list_imag = [np.array(x_minus_y_sym_m.imag)]

                    x_minus_y_sym_n = np.sqrt(2) * 0.5 * (x_minus_y_ao[n] +
                                                          x_minus_y_ao[n].T)
                    x_minus_y_sym_n_list_real = [np.array(x_minus_y_sym_n.real)]
                    x_minus_y_sym_n_list_imag = [np.array(x_minus_y_sym_n.imag)]

                else:
                    rhow_dm_sym_list_real = None
                    rhow_dm_sym_list_imag = None
                    x_minus_y_sym_m_list_real = None
                    x_minus_y_sym_m_list_imag = None
                    x_minus_y_sym_n_list_real = None
                    x_minus_y_sym_n_list_imag = None

                rhow_dm_sym_list_real = self.comm.bcast(rhow_dm_sym_list_real,
                                                        root=mpi_master())
                rhow_dm_sym_list_imag = self.comm.bcast(rhow_dm_sym_list_imag,
                                                        root=mpi_master())
                x_minus_y_sym_m_list_real = self.comm.bcast(x_minus_y_sym_m_list_real,
                                                            root=mpi_master())
                x_minus_y_sym_m_list_imag = self.comm.bcast(x_minus_y_sym_m_list_imag,
                                                            root=mpi_master())
                x_minus_y_sym_n_list_real = self.comm.bcast(x_minus_y_sym_n_list_real,
                                                            root=mpi_master())
                x_minus_y_sym_n_list_imag = self.comm.bcast(x_minus_y_sym_n_list_imag,
                                                            root=mpi_master())

                polgrad_xcgrad = self.calculate_xc_mn_contrib_complex(
                    molecule, ao_basis, rhow_dm_sym_list_real,
                    rhow_dm_sym_list_imag,
                    x_minus_y_sym_m_list_real,
                    x_minus_y_sym_n_list_real,
                    x_minus_y_sym_m_list_imag,
                    x_minus_y_sym_n_list_imag,
                    [gs_dm], xcfun_label,
                    profiler)

                #if self.rank == mpi_master():
                xc_pol_gradient[m, n] += polgrad_xcgrad

                if n != m:
                    xc_pol_gradient[n, m] = xc_pol_gradient[m, n]

        xc_pol_gradient = self.comm.reduce(xc_pol_gradient, root=mpi_master())

        return xc_pol_gradient

    def calculate_xc_mn_contrib_real(self, molecule, ao_basis, rhow_den,
                                     x_minus_y_den_m, x_minus_y_den_n,
                                     gs_density, xcfun_label, profiler):
        """
        Calculates exchange-correlation contribution to the (m,n) component of
        the real polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den:
            The perturbed density.
        :param x_minus_y_den(_m/n):
            The X-Y density.
        :param gs_density:
            The ground state density.
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The exchange-correlation contribution to polarizability gradient.
        """

        mol_grid = self._scf_drv._mol_grid

        xcgrad_drv = XCMolecularGradient()

        polgrad_xcgrad = xcgrad_drv.integrate_vxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, mol_grid, xcfun_label)

        polgrad_xcgrad += xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, gs_density, mol_grid,
            xcfun_label)

        polgrad_xcgrad += 0.5 * xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, x_minus_y_den_m, x_minus_y_den_n, gs_density,
            mol_grid, xcfun_label)

        polgrad_xcgrad += 0.5 * xcgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, x_minus_y_den_m, x_minus_y_den_n, gs_density,
            mol_grid, xcfun_label)

        polgrad_xcgrad += 0.5 * xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, x_minus_y_den_n, x_minus_y_den_m, gs_density,
            mol_grid, xcfun_label)

        polgrad_xcgrad += 0.5 * xcgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, x_minus_y_den_n, x_minus_y_den_m, gs_density,
            mol_grid, xcfun_label)

        #polgrad_xcgrad = self.comm.reduce(polgrad_xcgrad, root=mpi_master())

        return polgrad_xcgrad

    def calculate_xc_mn_contrib_complex(self, molecule, ao_basis, rhow_den_real,
                                        rhow_den_imag, x_minus_y_den_real_m,
                                        x_minus_y_den_real_n, x_minus_y_den_imag_m,
                                        x_minus_y_den_imag_n, gs_density, xcfun_label,
                                        profiler):
        """
        Calculates exchange-correlation contribution to the (m,n) component of
        the complex polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            ehe AO basis set.
        :param rhow_den(_real/_imag):
            The (Real/Imaginary) perturbed density.
        :param x_minus_y_den(_real/_imag):
            The (Real/Imaginary) X-Y density.
        :param gs_density:
            The ground state density.
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The exchange-correlation contribution to complex polarizability gradient.
        """

        mol_grid = self._scf_drv._mol_grid

        xcgrad_drv = XCMolecularGradient()

        # real contribution
        polgrad_xcgrad_real = xcgrad_drv.integrate_vxc_gradient(  # Re DM
            molecule, ao_basis, rhow_den_real, gs_density, mol_grid,
            xcfun_label)
        polgrad_xcgrad_real += xcgrad_drv.integrate_fxc_gradient(  # Re DM
            molecule, ao_basis, rhow_den_real, gs_density, gs_density, mol_grid,
            xcfun_label)

        polgrad_xcgrad_real += 0.5 * xcgrad_drv.integrate_fxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5 * xcgrad_drv.integrate_fxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5 * xcgrad_drv.integrate_kxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5 * xcgrad_drv.integrate_kxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)

        polgrad_xcgrad_real -= 0.5 * xcgrad_drv.integrate_fxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5 * xcgrad_drv.integrate_fxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5 * xcgrad_drv.integrate_kxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5 * xcgrad_drv.integrate_kxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)

        #polgrad_xcgrad_real = self.comm.reduce(polgrad_xcgrad_real,
        #                                       root=mpi_master())

        # imaginary contribution
        polgrad_xcgrad_imag = xcgrad_drv.integrate_vxc_gradient(  # Im DM
            molecule, ao_basis, rhow_den_imag, gs_density, mol_grid,
            xcfun_label)
        polgrad_xcgrad_imag += xcgrad_drv.integrate_fxc_gradient(  # Im DM
            molecule, ao_basis, rhow_den_imag, gs_density, gs_density, mol_grid,
            xcfun_label)

        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ReIm
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ReIm
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ReIm
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ReIm
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)

        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ImRe
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ImRe
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ImRe
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_imag += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ImRe
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)

        #polgrad_xcgrad_imag = self.comm.reduce(polgrad_xcgrad_imag,
        #                                       root=mpi_master())

        #if self.rank == mpi_master():
        #    return polgrad_xcgrad_real + 1j * polgrad_xcgrad_imag
        #else:
        #    return None
        return polgrad_xcgrad_real + 1j * polgrad_xcgrad_imag

    def get_fock_type_and_x_frac(self):
        """
        Determines fock type and fraction of exact exchange.

        :return fock_type:
            Fock type for fock build.
        :return exchange_scaling_factor:
            The fraction of exact exchange.
        """

        if self._dft:
            # TODO: range-separated Fock
            if self.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            exchange_scaling_factor = 1.0
            fock_type = "2jk"

        return fock_type, exchange_scaling_factor

    def calculate_unrel_dm(self, molecule, scf_tensors, x_plus_y_mo, x_minus_y_mo):
        """
        Calculates the symmetrized unrelaxed one-particle density matrix
        in AO basis.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param x_plus_y_mo:
            The X+Y response vectors in MO basis.
        :param x_minus_y_mo:
            The X-Y response vectors in MO basis.

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

        # calculate the symmetrized unrelaxed one-particle density matrix
        # in MO basis
        dm_oo = np.zeros((dof, dof, nocc, nocc), dtype=self.grad_dt)
        dm_vv = np.zeros((dof, dof, nvir, nvir), dtype=self.grad_dt)

        for x in range(dof):
            for y in range(x, dof):
                dm_vv[x, y] = 0.25 * (
                    # xib,yia->xyab
                    np.linalg.multi_dot([x_plus_y_mo[y].T, x_plus_y_mo[x]])
                    # xib,yia->xyab
                    + np.linalg.multi_dot([x_minus_y_mo[y].T, x_minus_y_mo[x]])
                    # yib,xia->xyab
                    + np.linalg.multi_dot([x_plus_y_mo[x].T, x_plus_y_mo[y]])
                    # yib,xia->xyab
                    + np.linalg.multi_dot([x_minus_y_mo[x].T, x_minus_y_mo[y]]))

                dm_oo[x, y] = -0.25 * (
                    # xja,yia->xyij
                    np.linalg.multi_dot([x_plus_y_mo[x], x_plus_y_mo[y].T])
                    # xja,yia->xyij
                    + np.linalg.multi_dot([x_minus_y_mo[x], x_minus_y_mo[y].T])
                    # yja,xia->xyij
                    + np.linalg.multi_dot([x_plus_y_mo[y], x_plus_y_mo[x].T])
                    # yja,xia->xyij
                    + np.linalg.multi_dot([x_minus_y_mo[y], x_minus_y_mo[x].T]))

                if y != x:
                    dm_vv[y,x] = dm_vv[x,y]
                    dm_oo[y,x] = dm_oo[x,y]

        # transform to AO basis: mi,xia,na->xmn
        unrel_dm_ao = np.zeros((dof, dof, nao, nao), dtype=self.grad_dt)
        for x in range(dof):
            for y in range(x, dof):
                unrel_dm_ao[x, y] = (
                    # mi,xyij,nj->xymn
                    np.linalg.multi_dot([mo_occ, dm_oo[x, y], mo_occ.T])
                    # ma,xyab,nb->xymn
                    + np.linalg.multi_dot([mo_vir, dm_vv[x, y], mo_vir.T]))

                if y != x:
                    unrel_dm_ao[y, x] = unrel_dm_ao[x, y]

        return unrel_dm_ao

    def construct_numerical_gradient(self, molecule, ao_basis, scf_drv, lr_drv):
        """
        Constructs the numerical polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_driver:
            The SCF driver.
        :param lr_drv:
            The linear response/CPP driver.

        :return num_polgradient:
            The numerical polarizability gradient.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()

        # number of frequencies
        n_freqs = len(self.frequencies)

        # array: numerical polarizability gradient
        num_polgradient = np.zeros((n_freqs, 3, 3, 3 * natm), dtype=self.grad_dt)

        # timings
        loop_start_time = tm.time()

        # dictionary for polarizability gradients
        polgrad_results = {}

        for i in range(natm):
            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                scf_drv.compute(new_mol, ao_basis)
                lr_drv._is_converged = False
                lr_results_p1 = lr_drv.compute(new_mol, ao_basis,
                                               scf_drv.scf_tensors)

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                scf_drv.compute(new_mol, ao_basis)
                lr_drv._is_converged = False
                lr_results_m1 = lr_drv.compute(new_mol, ao_basis,
                                              scf_drv.scf_tensors)
                # reset coordinates
                coords[i, d] += self.delta_h

                if self.do_four_point:
                    coords[i, d] += 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_p2 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] -= 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_m2 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)
                    # reset coordinates
                    coords[i, d] += 2.0 * self.delta_h

                    # construct gradient
                    for f, w in enumerate(self.frequencies):
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                # f'(x) ~ [ f(x - 2h) - 8 f(x - h)
                                # + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = ((
                                        lr_results_m2['response_functions'][key]
                                        - 8.0 *
                                        lr_results_m1['response_functions'][key]
                                        + 8.0 *
                                        lr_results_p1['response_functions'][key]
                                        -
                                        lr_results_p2['response_functions'][key]
                                    ) / (12.0 * self.delta_h))
                else:
                    for f, w in enumerate(self.frequencies):
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = ((
                                        lr_results_p1['response_functions'][key]
                                        -
                                        lr_results_m1['response_functions'][key])
                                        / (2.0 * self.delta_h))

        if self.rank == mpi_master():
            # convert array to dict
            for f, w in enumerate(self.frequencies):
                polgrad_results[w] = num_polgradient[f]

            valstr = '** Time spent on constructing the numerical gradient for '
            valstr += f'{n_freqs:d} frequencies: {(tm.time() - loop_start_time):.6f} sec **'
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        return polgrad_results

    def print_header(self):
        """
        Prints polarizability gradient calculation setup details to output stream.
        """

        str_width = 70

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.flush()

        cur_str = 'Polarizability gradient type    : '
        if self.is_complex:
            cur_str += 'Complex '
            cur_str2 = 'Damping value                   : '
            cur_str2 += f'{self.damping:.4f} a.u.'
        else:
            cur_str += 'Real '
        if self.numerical:
            cur_str += 'Numerical'
            cur_str3 = 'Numerical Method                : '
            if self.do_four_point:
                cur_str3 += 'Five-Point Stencil'
            else:
                cur_str3 += 'Symmetric Difference Quotient'
            cur_str4 = 'Finite Difference Step Size     : '
            cur_str4 += f'{self.delta_h:.4f} a.u.'
        else:
            cur_str += 'Analytical'

        self.ostream.print_blank()
        self.ostream.print_header(cur_str.ljust(str_width))
        if self.is_complex:
            self.ostream.print_header(cur_str2.ljust(str_width))

        if self.numerical:
            self.ostream.print_header(cur_str3.ljust(str_width))
            self.ostream.print_header(cur_str4.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = f'Molecular Grid Level            : {grid_level}'
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_polarizability_gradient(self, molecule):
        """
        Prints the polarizability gradient.

        :param molecule:
            The molecule.
        """

        labels = molecule.get_labels()

        natm = molecule.number_of_atoms()

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        # title for result output
        header = title + 'Polarizability Gradients'
        line = '=' * 40
        info_on_threshold = 'Printing only components with absolute values > 1.0e-6'

        self.ostream.print_blank()
        self.ostream.print_header(header)
        self.ostream.print_header(line)
        self.ostream.print_blank()
        self.ostream.print_info(info_on_threshold)
        self.ostream.print_blank()
        self.ostream.flush()

        for i in range(natm):
            atom_info = f'** Atom #{i + 1}: {labels[i]} **'
            self.ostream.print_header(atom_info)
            self.ostream.print_blank()

            if self.is_complex:
                column_headers = '{:<14s} {:>14s} {:>15s} {:>18s}'.format(
                    'Component', 'Frequency', 'Real', 'Imaginary')
                column_headers += '\n' + '-' * len(column_headers) + '\n'
                gradient_block = column_headers
                for w in self.frequencies:
                    current_gradient = self.polgradient[w].reshape(3, 3, natm, 3)
                    for aop, acomp in enumerate('xyz'):
                        for bop, bcomp in enumerate('xyz'):
                            for cop, ccomp in enumerate('xyz'):
                                row = ''
                                grad_element = current_gradient[aop, bop, i, cop]
                                if (abs(grad_element.real) < 1.0e-6) and (abs(grad_element.imag) < 1.0e-6):
                                    continue
                                grad_str = 'd<<{:>1s};{:<1s}>>/d{:<3s} {:>10.4f} + i{:<8} '.format(
                                        acomp.lower(), bcomp.lower(), ccomp.lower(), w, self.damping)
                                result = '{:>12.6f} {:>14.6f}'.format(round(grad_element.real,6), grad_element.imag)
                                row += grad_str + result + '\n'
                                gradient_block += row
            else:
                column_headers = '{:<14s} {:>12s} {:>15s}'.format(
                    'Component', 'Frequency', 'Value')
                column_headers += '\n' + '-' * len(column_headers) + '\n'
                gradient_block = column_headers
                for w in self.frequencies:
                    current_gradient = self.polgradient[w].reshape(3, 3, natm, 3)
                    for aop, acomp in enumerate('xyz'):
                        for bop, bcomp in enumerate('xyz'):
                            for cop, ccomp in enumerate('xyz'):
                                row = ''
                                grad_element = current_gradient[aop, bop, i, cop]
                                if (abs(grad_element.real) < 1.0e-6) and (abs(grad_element.imag) < 1.0e-6):
                                    continue
                                grad_str = 'd<<{:>1s};{:<1s}>>/d{:<3s} {:>10.4f} '.format(
                                        acomp.lower(), bcomp.lower(), ccomp.lower(), w)
                                result = '{:>18.6f}'.format(round(grad_element,6))
                                row += grad_str + result + '\n'
                                gradient_block += row

            self.ostream.print_block(gradient_block)
            self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.flush()

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def _check_real_or_complex_input(self, lr_results):
        """
        Checks if the input LR results and polgrad settings are
        both real or both complex.

        :param lr_results:
            The results of the linear response calculation.
        """

        response_functions = lr_results.get('response_functions', None)
        keys = list(response_functions.keys())
        is_complex_response = response_functions[keys[0]].dtype is np.dtype('complex128')

        if is_complex_response != self.is_complex:
            error_text = 'Mismatch between LR results and polgrad settings!'
            error_text += 'One is complex, the other is not.'
            raise ValueError(error_text)

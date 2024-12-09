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
import sys
import math
from mpi4py import MPI

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCMolecularGradient
from .veloxchemlib import mpi_master, hartree_in_wavenumber, denmat

from .polorbitalresponse import PolOrbitalResponse
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .molecule import Molecule
from .outputstream import OutputStream
from .matrices import Matrices
from .griddriver import GridDriver
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check, polgrad_sanity_check
from .dftutils import get_default_grid_level
from .profiler import Profiler

# VLX integrals
from .veloxchemlib import ElectricDipoleMomentGeom100Driver
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver, NuclearPotentialGeom010Driver,
                           FockGeom1000Driver, ElectricDipoleMomentGeom100Driver)
from .veloxchemlib import partition_atoms, make_matrix, mat_t
from .veloxchemlib import T4CScreener

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import dipole_deriv


class PolarizabilityGradient():
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

        self._scf_drv = scf_drv

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

        if self._dft and (self.grid_level is None):
            self.grid_level = get_default_grid_level(self.xcfun) 

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

        # sanity check
        dft_sanity_check(self, 'compute')

        start_time = tm.time()

        if self.numerical:
            # compute
            self.compute_numerical(molecule, basis, self._scf_drv)
        else:
            # sanity checks linear response input
            if (self.rank == mpi_master()) and (lr_results is None):
                error_message = 'PolarizabilityGradient missing input: LR results'
                error_message += 'for analytical gradient'
                raise ValueError(error_message)
            if self.rank == mpi_master():
                polgrad_sanity_check(self, self.flag, lr_results)
                self.check_real_or_complex_input(lr_results)
            # compute
            self.compute_analytical(molecule, basis, scf_tensors, lr_results)

        if self.rank == mpi_master():
            self.print_geometry(molecule)
            if self.do_print_polgrad:
                self.print_polarizability_gradient(molecule)

            valstr = '*** Time spent in polarizability gradient driver: '
            valstr += '{:.6f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_analytical_pyscf(self, molecule, basis, scf_tensors, lr_results):
        """
        Performs calculation of the both real and complex analytical
        polarizability gradient using pySCF integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.
        """

        # get orbital response results
        all_orbrsp_results = self.compute_orbital_response(
            molecule, basis, scf_tensors, lr_results)

        # number of frequencies
        n_freqs = len(self.frequencies)

        # dictionary for polarizability gradient
        polgrad_results = {}

        # timings
        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):
            info_msg = 'Building gradient for frequency = {:4.3f}'.format(w)
            self.ostream.print_info(info_msg)
            self.ostream.print_blank()
            self.ostream.flush()

            if self.rank == mpi_master():
                orbrsp_results = all_orbrsp_results[w]
                mo = scf_tensors['C_alpha']  # only alpha part
                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]
                natm = molecule.number_of_atoms()
                gs_dm = scf_tensors['D_alpha']  # only alpha part

                x_plus_y = orbrsp_results['x_plus_y_ao']
                x_minus_y = orbrsp_results['x_minus_y_ao']

                dof = len(self.vector_components)

                # Lagrange multipliers
                omega_ao = orbrsp_results['omega_ao'].reshape(
                    dof, dof, nao, nao)
                lambda_ao = orbrsp_results['lambda_ao'].reshape(dof, dof, nao, nao)
                lambda_ao += lambda_ao.transpose(0, 1, 3, 2)  # vir-occ

                # calculate relaxed density matrix
                rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao

                # initiate polarizability gradient variable with data type set in init()
                pol_gradient = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

                # loop over atoms and contract integral derivatives with density matrices
                for i in range(natm):
                    # import integrals
                    d_hcore, d_ovlp, d_eri, d_dipole = self.import_integrals(
                        molecule, basis, i)

                    # construct polarizability gradient
                    pol_gradient[:,:,i,:] = self.construct_scf_polgrad(
                        gs_dm, rel_dm_ao, omega_ao, x_plus_y, x_minus_y,
                        d_hcore, d_ovlp, d_eri, d_dipole, nao, i)
            else:
                gs_dm = None
                rel_dm_ao = None
                x_minus_y = None
                pol_gradient = None

            gs_dm = self.comm.bcast(gs_dm, root=mpi_master())

            if self._dft:
                xcfun_label = self.xcfun.get_func_label()
                # compute the XC contribution
                polgrad_xc_contrib = self.compute_polgrad_xc_contrib(
                    molecule, basis, gs_dm, rel_dm_ao, x_minus_y, xcfun_label)
                # add contribution to the SCF polarizability gradient
                if self.rank == mpi_master():
                    pol_gradient += polgrad_xc_contrib

            if self.rank == mpi_master():
                polgrad_results[w] = pol_gradient.reshape(dof, dof, 3 * natm)

        if self.rank == mpi_master():
            self.polgradient = dict(polgrad_results)

            valstr = '** Time spent on constructing the analytical gradient for '
            valstr += '{:d} frequencies: {:.6f} sec **'.format(
                n_freqs, tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()
           
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

        # get orbital response results
        all_orbrsp_results = self.compute_orbital_response(
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

        # number of atomic orbitals
        nao = basis.get_dimensions_of_basis()

        for f, w in enumerate(self.frequencies):
            info_msg = 'Building gradient for frequency = {:4.3f}'.format(w)
            self.ostream.print_info(info_msg)
            self.ostream.print_blank()
            self.ostream.flush()

            if self.rank == mpi_master():
                orbrsp_results = all_orbrsp_results[w]
                #mo = scf_tensors['C_alpha']  # only alpha part
                #nao = mo.shape[0]
                #nocc = molecule.number_of_alpha_electrons()
                #mo_occ = mo[:, :nocc].copy()
                #mo_vir = mo[:, nocc:].copy()
                #nvir = mo_vir.shape[1]
                gs_dm = scf_tensors['D_alpha']  # only alpha part

                x_plus_y = orbrsp_results['x_plus_y_ao']
                x_minus_y = orbrsp_results['x_minus_y_ao']

                # Lagrange multipliers
                omega_ao = orbrsp_results['omega_ao'].reshape(
                    dof, dof, nao, nao)
                lambda_ao = orbrsp_results['lambda_ao'].reshape(dof, dof, nao, nao)
                lambda_ao += lambda_ao.transpose(0, 1, 3, 2)  # vir-occ

                # calculate relaxed density matrix
                rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao
            else:
                orbrsp_results = None
                gs_dm = None
                x_plus_y = None
                x_minus_y = None
                omega_ao = None
                rel_dm_ao = None

            orbrsp_results = self.comm.bcast(orbrsp_results, root=mpi_master())
            gs_dm = self.comm.bcast(gs_dm, root=mpi_master())
            x_plus_y= self.comm.bcast(x_plus_y, root=mpi_master())
            x_minus_y = self.comm.bcast(x_minus_y, root=mpi_master())
            omega_ao = self.comm.bcast(omega_ao, root=mpi_master())
            rel_dm_ao = self.comm.bcast(rel_dm_ao, root=mpi_master())

            # initiate polarizability gradient variable with data type set in init()
            pol_gradient = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

            # WIP
            #tmp_scf_gradient = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

            # kinetic energy gradient driver
            kin_grad_drv = KineticEnergyGeom100Driver()

            # nuclear potential gradienr drivers
            npot_grad_100_drv = NuclearPotentialGeom100Driver()
            npot_grad_010_drv = NuclearPotentialGeom010Driver()

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
                        #tmp_scf_gradient[x, y, iatom, icoord] += (
                        pol_gradient[x, y, iatom, icoord] += (
                                np.linalg.multi_dot([ # xymn,amn->xya
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
                        #tmp_scf_gradient[x, y, iatom, icoord] += (
                        pol_gradient[x, y, iatom, icoord] += (
                                1.0 * np.linalg.multi_dot([ # xymn,amn->xya
                                2.0 * omega_ao[x, y].reshape(nao**2),
                                      gmat_ovlp.reshape(nao**2)]))

                gmats_ovlp = Matrices()

                # dipole contribution to gradient
                dip_grad_drv = ElectricDipoleMomentGeom100Driver()
                gmats_dip = dip_grad_drv.compute(molecule, basis, [0.0, 0.0, 0.0], iatom)

                # the keys of the dipole gmat
                gmats_dip_components = (['X_X', 'X_Y', 'X_Z', 'Y_X' ]
                                     + [ 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'])

                # dictionary to convert from string idx to integer idx
                comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

                d_dipole = np.zeros((dof, 3, nao, nao))

                for i, label in enumerate(gmats_dip_components):
                    gmat_dip = gmats_dip.matrix_to_numpy(label)
                    gmat_dip += gmat_dip.T

                    icoord = comp_to_idx[label[0]] # atom coordinate component
                    icomp = comp_to_idx[label[-1]] # dipole operator component

                    # reorder indices to first is operator comp, second is coord
                    d_dipole[icomp, icoord] += gmat_dip

                for icoord in range(3):
                    # loop over operator components
                    for x, y in xy_pairs:
                        #tmp_scf_gradient[x, y, iatom, icoord] += (
                        pol_gradient[x, y, iatom, icoord] += (
                                - 2.0 * np.linalg.multi_dot([ # xmn,yamn->xya
                                    x_minus_y[x].reshape(nao**2),
                                    d_dipole[y, icoord].reshape(nao**2)
                                ]) - 2.0 * np.linalg.multi_dot([ # xmn,yamn->yxa
                                    d_dipole[x, icoord].reshape(nao**2),
                                    x_minus_y[y].reshape(nao**2)
                                ]))

                        # symmetric
                        #if (y != x):
                        #    tmp_scf_gradient[y, x, iatom, icoord] += tmp_scf_gradient[x, y, iatom, icoord]

                gmats_dip = Matrices()

            #tmp_scf_gradient = self.comm.reduce(tmp_scf_gradient, root=mpi_master())
            # ERI contribution
            eri_contrib = self.compute_eri_contrib(molecule, basis, gs_dm,
                                    rel_dm_ao, x_plus_y, x_minus_y, local_atoms)

            pol_gradient += eri_contrib
            pol_gradient = self.comm.reduce(pol_gradient, root=mpi_master())

            if self.rank == mpi_master():
                #pol_gradient = tmp_scf_gradient + eri_contrib
                for x in range(dof):
                    for y in range(x + 1, dof):
                        pol_gradient[y, x] += pol_gradient[x, y]

                #pol_gradient = tmp_scf_gradient + eri_contrib
            #else:
            #    gs_dm = None
            #    rel_dm_ao = None
            #    x_minus_y = None
            #    pol_gradient = None

            #gs_dm = self.comm.bcast(gs_dm, root=mpi_master())

            if self._dft:
                xcfun_label = self.xcfun.get_func_label()
                # compute the XC contribution
                polgrad_xc_contrib = self.compute_polgrad_xc_contrib(
                    molecule, basis, gs_dm, rel_dm_ao, x_minus_y, xcfun_label)

                # add contribution to the SCF polarizability gradient
                if self.rank == mpi_master():
                    pol_gradient += polgrad_xc_contrib

            if self.rank == mpi_master():
                polgrad_results[w] = pol_gradient.reshape(dof, dof, 3 * natm)

        if self.rank == mpi_master():
            self.polgradient = dict(polgrad_results)

            n_freqs = len(self.frequencies)

            valstr = '** Time spent on constructing the analytical gradient for '
            valstr += '{:d} frequencies: {:.6f} sec **'.format(
                n_freqs, tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_eri_contrib(self, molecule, basis, gs_dm, rel_dm_ao,
                                         x_plus_y, x_minus_y, local_atoms):
        """
        oDirects the computation of the contribution from ERI derivative integrals
        to the polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param gs_dm:
            The ground state density.
        :param rel_dm_ao:
            The relaxed density matrix.
        :param x_plus_y:
            The X+Y response vectors.
        :param x_minus_y:
            The X-Y response vectors.
        :param local_atoms:
            The atom partition for the MPI node.

        :return eri_deriv_contrib:
            The ERI derivative integral contribution to the polarizability gradient.
        """

        if self.is_complex:
            eri_deriv_contrib = self.compute_eri_contrib_complex(molecule, basis, gs_dm,
                                            rel_dm_ao, x_plus_y, x_minus_y, local_atoms)
        else:
            eri_deriv_contrib = self.compute_eri_contrib_real(molecule, basis, gs_dm,
                                            rel_dm_ao, x_plus_y, x_minus_y, local_atoms)

        return eri_deriv_contrib

    def compute_eri_contrib_real(self, molecule, basis, gs_dm, rel_dm_ao,
                                         x_plus_y, x_minus_y, local_atoms):
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
        :param x_plus_y:
            The X+Y response vectors.
        :param x_minus_y:
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
                den_mat_for_fock_rel.set_values(2.0*rel_dm_ao[x,y])
                # (X+Y)_x
                den_mat_for_fock_xpy_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x.set_values(x_plus_y[x])
                # (X+Y)_y
                den_mat_for_fock_xpy_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y.set_values(x_plus_y[y])
                # (X+Y)_x - (X+Y)_x
                den_mat_for_fock_xpy_m_xpyT_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x.set_values( x_plus_y[x]
                                                        - x_plus_y[x].T)
                # (X+Y)_y - (X+Y)_y
                den_mat_for_fock_xpy_m_xpyT_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y.set_values( x_plus_y[y]
                                                        - x_plus_y[y].T)
                # (X-Y)_x
                den_mat_for_fock_xmy_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x.set_values(x_minus_y[x])
                # (X-Y)_y
                den_mat_for_fock_xmy_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y.set_values(x_minus_y[y])
                # (X-Y)_x + (X-Y)_x
                den_mat_for_fock_xmy_p_xmyT_x = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x.set_values( x_minus_y[x]
                                                        + x_minus_y[x].T)
                # (X-Y)_y + (X-Y)_y
                den_mat_for_fock_xmy_p_xmyT_y = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y.set_values( x_minus_y[y]
                                                        + x_minus_y[y].T)
                # contraction of integrals and DMs
                erigrad_rel = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_gs,
                                             den_mat_for_fock_rel, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                erigrad_xpy_xy = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_x,
                                             den_mat_for_fock_xpy_m_xpyT_y, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                erigrad_xpy_yx = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_m_xpyT_x,
                                             den_mat_for_fock_xpy_y, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                erigrad_xmy_xy = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_x,
                                             den_mat_for_fock_xmy_p_xmyT_y, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                erigrad_xmy_yx = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_p_xmyT_x,
                                             den_mat_for_fock_xmy_y, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                eri_deriv_contrib[x, y, iatom] += np.array(erigrad_rel) * factor
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xpy_xy) * factor
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xpy_yx) * factor
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xmy_xy) * factor
                eri_deriv_contrib[x, y, iatom] += 0.5 * np.array(erigrad_xmy_yx) * factor

                # symmetric
                #if (y != x):
                #    eri_deriv_contrib[y, x, iatom] += eri_deriv_contrib[x, y, iatom]

        eri_deriv_contrib = self.comm.reduce(eri_deriv_contrib, root=mpi_master())

        return eri_deriv_contrib

    def compute_eri_contrib_complex(self, molecule, basis, gs_dm, rel_dm_ao,
                                         x_plus_y, x_minus_y, local_atoms):
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
        :param x_plus_y:
            The X+Y response vectors.
        :param x_minus_y:
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
        xy_pairs = [(x,y) for x in range(dof) for y in range(x,dof)]

        eri_deriv_contrib = np.zeros((dof, dof, natm, 3), dtype=self.grad_dt)

        for iatom in local_atoms:
            # screening
            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            # contraction with density matrices
            for x, y in xy_pairs:
                # relaxed DM
                den_mat_for_fock_rel_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel_real.set_values(2.0*rel_dm_ao[x,y].real)
                den_mat_for_fock_rel_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel_imag.set_values(2.0*rel_dm_ao[x,y].imag)
                # (X+Y)_x
                den_mat_for_fock_xpy_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x_real.set_values(x_plus_y[x].real)
                den_mat_for_fock_xpy_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_x_imag.set_values(x_plus_y[x].imag)
                # (X+Y)_y
                den_mat_for_fock_xpy_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y_real.set_values(x_plus_y[y].real)
                den_mat_for_fock_xpy_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_y_imag.set_values(x_plus_y[y].imag)
                # (X+Y)_x - (X+Y)_x
                den_mat_for_fock_xpy_m_xpyT_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x_real.set_values( (x_plus_y[x]
                                                        - x_plus_y[x].T).real )
                den_mat_for_fock_xpy_m_xpyT_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_x_imag.set_values( (x_plus_y[x]
                                                        - x_plus_y[x].T).imag )
                # (X+Y)_y - (X+Y)_y
                den_mat_for_fock_xpy_m_xpyT_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y_real.set_values( (x_plus_y[y]
                                                        - x_plus_y[y].T).real )
                den_mat_for_fock_xpy_m_xpyT_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xpy_m_xpyT_y_imag.set_values( (x_plus_y[y]
                                                        - x_plus_y[y].T).imag )
                # (X-Y)_x
                den_mat_for_fock_xmy_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x_real.set_values(x_minus_y[x].real)
                den_mat_for_fock_xmy_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_x_imag.set_values(x_minus_y[x].imag)
                # (X-Y)_y
                den_mat_for_fock_xmy_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y_real.set_values(x_minus_y[y].real)
                den_mat_for_fock_xmy_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_y_imag.set_values(x_minus_y[y].imag)
                # (X-Y)_x + (X-Y)_x
                den_mat_for_fock_xmy_p_xmyT_x_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x_real.set_values( (x_minus_y[x]
                                                        + x_minus_y[x].T).real)
                den_mat_for_fock_xmy_p_xmyT_x_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_x_imag.set_values( (x_minus_y[x]
                                                        + x_minus_y[x].T).imag)
                # (X-Y)_y + (X-Y)_y
                den_mat_for_fock_xmy_p_xmyT_y_real = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y_real.set_values( (x_minus_y[y]
                                                        + x_minus_y[y].T).real)
                den_mat_for_fock_xmy_p_xmyT_y_imag = make_matrix(basis, mat_t.general)
                den_mat_for_fock_xmy_p_xmyT_y_imag.set_values( (x_minus_y[y]
                                                        + x_minus_y[y].T).imag)
                # contraction of integrals and DMs
                # Re
                erigrad_rel_re = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_gs,
                                             den_mat_for_fock_rel_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # Im
                erigrad_rel_im = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_gs,
                                             den_mat_for_fock_rel_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReRe
                erigrad_xpy_xy_rere = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_x_real,
                                             den_mat_for_fock_xpy_m_xpyT_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImIm
                erigrad_xpy_xy_imim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_x_imag,
                                             den_mat_for_fock_xpy_m_xpyT_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReIm
                erigrad_xpy_xy_reim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_x_real,
                                             den_mat_for_fock_xpy_m_xpyT_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImRe
                erigrad_xpy_xy_imre = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_x_imag,
                                             den_mat_for_fock_xpy_m_xpyT_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReRe
                erigrad_xpy_yx_rere = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_m_xpyT_x_real,
                                             den_mat_for_fock_xpy_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImIm
                erigrad_xpy_yx_imim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_m_xpyT_x_imag,
                                             den_mat_for_fock_xpy_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReIm
                erigrad_xpy_yx_reim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_m_xpyT_x_real,
                                             den_mat_for_fock_xpy_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImRe
                erigrad_xpy_yx_imre = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xpy_m_xpyT_x_imag,
                                             den_mat_for_fock_xpy_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReRe
                erigrad_xmy_xy_rere = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_x_real,
                                             den_mat_for_fock_xmy_p_xmyT_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImIm
                erigrad_xmy_xy_imim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_x_imag,
                                             den_mat_for_fock_xmy_p_xmyT_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReIm
                erigrad_xmy_xy_reim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_x_real,
                                             den_mat_for_fock_xmy_p_xmyT_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImRe
                erigrad_xmy_xy_imre = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_x_imag,
                                             den_mat_for_fock_xmy_p_xmyT_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReRe
                erigrad_xmy_yx_rere = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_p_xmyT_x_real,
                                             den_mat_for_fock_xmy_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImIm
                erigrad_xmy_yx_imim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_p_xmyT_x_imag,
                                             den_mat_for_fock_xmy_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ReIm
                erigrad_xmy_yx_reim = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_p_xmyT_x_real,
                                             den_mat_for_fock_xmy_y_imag, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                # ImRe
                erigrad_xmy_yx_imre = fock_grad_drv.compute(basis, screener_atom,
                                             screener,
                                             den_mat_for_fock_xmy_p_xmyT_x_imag,
                                             den_mat_for_fock_xmy_y_real, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)

                # Real fmat contribution to gradient
                erigrad_real = np.array(erigrad_rel_re) # real DM
                erigrad_real += 0.5 * np.array(erigrad_xpy_xy_rere) # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xpy_xy_imim) # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xpy_yx_rere) # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xpy_yx_imim) # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xmy_xy_rere) # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xmy_xy_imim) # ImIm
                erigrad_real += 0.5 * np.array(erigrad_xmy_yx_rere) # ReRe
                erigrad_real -= 0.5 * np.array(erigrad_xmy_yx_imim) # ImIm
                erigrad_real *= factor

                # Imaginary fmat contribution to gradient
                erigrad_imag = np.array(erigrad_rel_im) # imag DM
                erigrad_imag += 0.5 * np.array(erigrad_xpy_xy_reim) # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xpy_xy_imre) # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xpy_yx_reim) # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xpy_yx_imre) # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xmy_xy_reim) # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xmy_xy_imre) # ImRe
                erigrad_imag += 0.5 * np.array(erigrad_xmy_yx_reim) # ReIm
                erigrad_imag += 0.5 * np.array(erigrad_xmy_yx_imre) # ImRe
                erigrad_imag *= factor

                # add to complex variable
                eri_deriv_contrib[x, y, iatom] += erigrad_real + 1j * erigrad_imag

                # symmetric
                #if (y != x):
                #    eri_deriv_contrib[y, x, iatom] += eri_deriv_contrib[x, y, iatom]

        eri_deriv_contrib = self.comm.reduce(eri_deriv_contrib, root=mpi_master())

        return eri_deriv_contrib

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

        # dictionary for results
        polgrad_results = {}

        # construct polarizability gradient
        num_polgradient = self.construct_numerical_gradient(molecule, ao_basis, scf_drv, lr_drv)

        if self.rank == mpi_master():
            for f, w in enumerate(self.frequencies):
                polgrad_results[w] = num_polgradient[f]
            self.polgradient = dict(polgrad_results)

        # unmute the output streams
        scf_drv.ostream.unmute()
        lr_drv.ostream.unmute()

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
            valstr = '** Time spent on orbital response for {} frequencies: '.format(
                len(self.frequencies))
            valstr += '{:.6f} sec **'.format(tm.time() - orbrsp_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        return orbrsp_drv.cphf_results

    def import_integrals(self, molecule, basis, idx):
        """
        Imports integrals from PySCF for analytical polarizability
        gradient.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param idx:
            The atom index

        :return d_hcore:
            The derivative H_core integrals.
        :return d_ovlp:
            The derivative overlap integrals.
        :return d_eri:
            The derivative electron-repulsion integrals.
        :return d_dipole:
            The derivative dipole moment integrals.
        """

        # timings
        integral_start_time = tm.time()

        # importing integral derivatives from pyscf
        d_hcore = hcore_deriv(molecule, basis, idx)
        d_ovlp = overlap_deriv(molecule, basis, idx)
        d_eri = eri_deriv(molecule, basis, idx)
        d_dipole = dipole_deriv(molecule, basis, idx)

        valstr = ' * Time spent importing integrals for atom #{}: '.format(
            idx + 1)
        valstr += '{:.6f} sec * '.format(tm.time() - integral_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        return d_hcore, d_ovlp, d_eri, d_dipole

    def construct_scf_polgrad(self, gs_dm, rel_dm_ao, omega_ao, x_plus_y, x_minus_y,
                              d_hcore, d_ovlp, d_eri, d_dipole, nao, idx):
        """
        Constructs the SCF polarizability gradient

        :param gs_dm:
            The ground state density matrix.
        :param rel_dm_ao:
            The relaxed one-particle density matrix in AO basis.
        :param omega_ao:
            The omega Lagrangian multipliers.
        :param x_plus_y:
            The X+Y response vectors.
        :param x_minus_y:
            The X-Y response vectors.
        :param d_hcore:
            The derivative H_core integrals.
        :param d_ovlp:
            The derivative overlap integrals.
        :param d_eri:
            The derivative electron-repulsion integrals.
        :param d_dipole:
            The derivative dipole moment integrals.
        :param nao:
            The number of atomic orbitals.
        :param idx:
            The index of the atom.

        :return scf_polgrad:
            The SCF polarizability gradient.
        """

        dof = len(self.vector_components)

        #frac_K = self.get_k_fraction()
        dummy, frac_K = self.get_fock_type_and_x_frac()

        scf_polgrad = np.zeros((dof, dof, 3), dtype=self.grad_dt)

        gradient_start_time = tm.time()

        # construct the analytic polarizability gradient
        for x in range(dof):
            for y in range(x, dof):
                for a in range(3): # nuclear cartesian coordinate components
                    scf_polgrad[x, y, a] += (
                        np.linalg.multi_dot([ # xymn,amn->xya
                            2.0 * rel_dm_ao[x, y].reshape(nao**2),
                            d_hcore[a].reshape(nao**2)
                        ]) + 1.0 * np.linalg.multi_dot([ # xymn,amn->xya
                            2.0 * omega_ao[x, y].reshape(nao**2),
                            d_ovlp[a].reshape(nao * nao)
                        ]) + 2.0 * np.linalg.multi_dot([ # mt,xynp,amtnp->xya
                            gs_dm.reshape(nao**2), d_eri[a].reshape(
                                nao**2, nao**2),
                            2.0 * rel_dm_ao[x, y].reshape(nao**2)
                        ]) - 1.0 * frac_K * np.linalg.multi_dot([ # mt,xynp,amnpt->xya
                            gs_dm.reshape(nao**2),
                            d_eri[a].transpose(0, 3, 1, 2).reshape(
                                nao**2, nao**2),
                            2.0 * rel_dm_ao[x, y].reshape(nao**2)
                        ]) + 1.0 * np.linalg.multi_dot([ # xmn,ypt,atpmn->xya
                            x_plus_y[x].reshape(nao**2),
                            d_eri[a].transpose(2, 3, 1, 0).reshape(
                                nao**2, nao**2),
                            (x_plus_y[y] - x_plus_y[y].T).reshape(
                                nao**2)
                        ]) - 0.5 * frac_K * np.linalg.multi_dot([ # xmn,ypt,atnmp->xya
                            x_plus_y[x].reshape(nao**2),
                            d_eri[a].transpose(2, 1, 3, 0).reshape(
                                nao**2, nao**2),
                            (x_plus_y[y] - x_plus_y[y].T).reshape(
                                nao**2)
                        ]) + 1.0 * np.linalg.multi_dot([ # xmn,ypt,atpmn->xya
                            x_minus_y[x].reshape(nao**2),
                            d_eri[a].transpose(2, 3, 1, 0).reshape(
                                nao**2, nao**2),
                            (x_minus_y[y] + x_minus_y[y].T).reshape(
                                nao**2)
                        ]) - 0.5 * frac_K * np.linalg.multi_dot([ # xmn,ypt,atnmp->xya
                            x_minus_y[x].reshape(nao**2),
                            d_eri[a].transpose(2, 1, 3, 0).reshape(
                                nao**2, nao**2),
                            (x_minus_y[y] + x_minus_y[y].T).reshape(
                                nao**2)
                        ]) - 2.0 * np.linalg.multi_dot([ # xmn,yamn->xya
                            x_minus_y[x].reshape(nao**2),
                            d_dipole[y, a].reshape(nao**2)
                        ])
                        + 1.0 * np.linalg.multi_dot([ # xmn,ypt,atpmn->yxa
                            (x_plus_y[x] - x_plus_y[x].T
                            ).reshape(nao**2), d_eri[a].transpose(
                                1, 0, 2, 3).reshape(nao**2, nao**2),
                            x_plus_y[y].reshape(nao**2)
                        ]) - 0.5 * frac_K * np.linalg.multi_dot([ # xmn,ypt,atnmp->yxa
                            (x_plus_y[x] - x_plus_y[x].T
                            ).reshape(nao**2), d_eri[a].transpose(
                                3, 0, 2, 1).reshape(nao**2, nao**2),
                            x_plus_y[y].reshape(nao**2)
                        ]) + 1.0 * np.linalg.multi_dot([ # xmn,ypt,atpmn->yxa
                            (x_minus_y[x] + x_minus_y[x].T
                            ).reshape(nao**2), d_eri[a].transpose(
                                1, 0, 2, 3).reshape(nao**2, nao**2),
                            x_minus_y[y].reshape(nao**2)
                        ]) - 0.5 * frac_K * np.linalg.multi_dot([ # xmn,ypt,atnmp->yxa
                            (x_minus_y[x] + x_minus_y[x].T
                            ).reshape(nao**2), d_eri[a].transpose(
                                3, 0, 2, 1).reshape(nao**2, nao**2),
                            x_minus_y[y].reshape(nao**2)
                        ]) - 2.0 * np.linalg.multi_dot([ # xmn,yamn->yxa
                            d_dipole[x, a].reshape(nao**2),
                            x_minus_y[y].reshape(nao**2)
                        ]))

                    if (y != x):
                        scf_polgrad[y, x, a] += scf_polgrad[x, y, a]

        valstr = ' * Time spent constructing pol. gradient for '
        valstr += 'atom #{:d}: {:.6f} sec * '.format(
            (idx + 1),
            tm.time() - gradient_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        return scf_polgrad

    def compute_polgrad_xc_contrib(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                    x_minus_y, xcfun_label):
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
        :param x_minus_y:
            The X-Y response vector.
        :param xcfun_label:
            The label for the XC functional

        :return xc_contrib:
            The XC-contribution to the polarizability gradient
        """

        if self.is_complex:
            xc_contrib = self.compute_polgrad_xc_contrib_complex(
                molecule, ao_basis, gs_dm, rel_dm_ao, x_minus_y, xcfun_label)
        else:
            xc_contrib = self.compute_polgrad_xc_contrib_real(
                molecule, ao_basis, gs_dm, rel_dm_ao, x_minus_y, xcfun_label)

        return xc_contrib

    def compute_polgrad_xc_contrib_real(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                    x_minus_y, xcfun_label):
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
        :param x_minus_y:
            The X-Y response vector.
        :param xcfun_label:
            The label for the XC functional

        :return xc_pol_gradient:
            The XC-contribution to the polarizability gradient
        """

        natm = molecule.number_of_atoms()
        dof = len(self.vector_components)
        xc_pol_gradient = np.zeros((dof, dof, natm, 3))

        for m in range(dof):
            for n in range(m, dof):
                if self.rank == mpi_master():
                    rhow_dm = 1.0 * rel_dm_ao[m, n]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)

                    # symmetrize
                    x_minus_y_sym_m = np.sqrt(2) * 0.5 * (x_minus_y[m] +
                                                          x_minus_y[m].T)
                    x_minus_y_sym_n = np.sqrt(2) * 0.5 * (x_minus_y[n] +
                                                          x_minus_y[n].T)

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
                    [gs_dm], xcfun_label)

                if self.rank == mpi_master():
                    xc_pol_gradient[m, n] += polgrad_xcgrad

                if (n != m):
                    xc_pol_gradient[n, m] = xc_pol_gradient[m, n]


        return xc_pol_gradient

    def compute_polgrad_xc_contrib_complex(self, molecule, ao_basis, gs_dm, rel_dm_ao,
                                    x_minus_y, xcfun_label):
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
        :param x_minus_y:
            The X-Y response vector.
        :param xcfun_label:
            The label for the XC functional

        :return xc_pol_gradient:
            The XC-contribution to the polarizability gradient
        """

        natm = molecule.number_of_atoms()
        dof = len(self.vector_components)
        xc_pol_gradient = np.zeros((dof, dof, natm, 3), dtype=np.dtype('complex128'))

        # TODO change "m,n" to "x,y" for consistency
        for m in range(dof):
            for n in range(m, dof):
                if self.rank == mpi_master():
                    rhow_dm = 1.0 * rel_dm_ao[m, n]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)

                    rhow_dm_sym_list_real = [np.array(rhow_dm_sym.real)]
                    rhow_dm_sym_list_imag = [np.array(rhow_dm_sym.imag)]

                    # symmetrize
                    x_minus_y_sym_m = np.sqrt(2) * 0.5 * (x_minus_y[m] +
                                                          x_minus_y[m].T)
                    x_minus_y_sym_m_list_real = [ np.array(x_minus_y_sym_m.real) ]
                    x_minus_y_sym_m_list_imag = [ np.array(x_minus_y_sym_m.imag) ]

                    x_minus_y_sym_n = np.sqrt(2) * 0.5 * (x_minus_y[n] +
                                                          x_minus_y[n].T)
                    x_minus_y_sym_n_list_real = [np.array(x_minus_y_sym_n.real)]
                    x_minus_y_sym_n_list_imag = [np.array(x_minus_y_sym_n.imag)]

                else:
                    rhow_dm_sym_list_real = None
                    rhow_dm_sym_list_imag = None
                    x_minus_y_sym_m_list_real = None
                    x_minus_y_sym_m_list_imag = None
                    x_minus_y_sym_n_list_real = None
                    x_minus_y_sym_n_list_imag = None

                rhow_dm_sym_list_real = self.comm.bcast(rhow_dm_sym_list_real, root=mpi_master())
                rhow_dm_sym_list_imag = self.comm.bcast(rhow_dm_sym_list_imag, root=mpi_master())
                x_minus_y_sym_m_list_real = self.comm.bcast(x_minus_y_sym_m_list_real, root=mpi_master())
                x_minus_y_sym_m_list_imag = self.comm.bcast(x_minus_y_sym_m_list_imag, root=mpi_master())
                x_minus_y_sym_n_list_real = self.comm.bcast(x_minus_y_sym_n_list_real, root=mpi_master())
                x_minus_y_sym_n_list_imag = self.comm.bcast(x_minus_y_sym_n_list_imag, root=mpi_master())

                polgrad_xcgrad = self.calculate_xc_mn_contrib_complex(
                    molecule, ao_basis, rhow_dm_sym_list_real,
                    rhow_dm_sym_list_imag,
                    x_minus_y_sym_m_list_real,
                    x_minus_y_sym_n_list_real,
                    x_minus_y_sym_m_list_imag,
                    x_minus_y_sym_n_list_imag,
                    [gs_dm], xcfun_label)

                if self.rank == mpi_master():
                    xc_pol_gradient[m, n] += polgrad_xcgrad

                if (n != m):
                    xc_pol_gradient[n, m] = xc_pol_gradient[m, n]

        return xc_pol_gradient

    def calculate_xc_mn_contrib_real(self, molecule, ao_basis, rhow_den, x_minus_y_den_m,
                                     x_minus_y_den_n, gs_density, xcfun_label):
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

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

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

        polgrad_xcgrad = self.comm.reduce(polgrad_xcgrad, root=mpi_master())

        return polgrad_xcgrad

    def calculate_xc_mn_contrib_complex(self, molecule, ao_basis, rhow_den_real,
                                        rhow_den_imag, x_minus_y_den_real_m,
                                        x_minus_y_den_real_n, x_minus_y_den_imag_m,
                                        x_minus_y_den_imag_n, gs_density, xcfun_label):
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

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xcgrad_drv = XCMolecularGradient()

        # real contribution
        polgrad_xcgrad_real = xcgrad_drv.integrate_vxc_gradient(  # Re DM
            molecule, ao_basis, rhow_den_real, gs_density, mol_grid,
            xcfun_label)
        polgrad_xcgrad_real += xcgrad_drv.integrate_fxc_gradient(  # Re DM
            molecule, ao_basis, rhow_den_real, gs_density, gs_density, mol_grid,
            xcfun_label)

        polgrad_xcgrad_real += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5*xcgrad_drv.integrate_fxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_m, x_minus_y_den_real_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real += 0.5*xcgrad_drv.integrate_kxc_gradient(  # ReRe
            molecule, ao_basis, x_minus_y_den_real_n, x_minus_y_den_real_m,
            gs_density, mol_grid, xcfun_label)

        polgrad_xcgrad_real -= 0.5*xcgrad_drv.integrate_fxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5*xcgrad_drv.integrate_fxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5*xcgrad_drv.integrate_kxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_m, x_minus_y_den_imag_n,
            gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad_real -= 0.5*xcgrad_drv.integrate_kxc_gradient(  # ImIm
            molecule, ao_basis, x_minus_y_den_imag_n, x_minus_y_den_imag_m,
            gs_density, mol_grid, xcfun_label)

        polgrad_xcgrad_real = self.comm.reduce(polgrad_xcgrad_real,
                                               root=mpi_master())

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

        polgrad_xcgrad_imag = self.comm.reduce(polgrad_xcgrad_imag,
                                               root=mpi_master())

        if self.rank == mpi_master():
            return polgrad_xcgrad_real + 1j * polgrad_xcgrad_imag
        else:
            return None

    def get_k_fraction(self):
        """
        Determines fraction of exact exchange, prefactor for K.

        :return frac_k:
            The fraction
        """
        if self._dft:
            if self.xcfun.is_hybrid():
                frac_k = self.xcfun.get_frac_exact_exchange()
            else:
                frac_k = 0.0
        else:
            frac_k = 1.0

        return frac_k

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
                                    num_polgradient[f, aop, bop, 3 * i + d] = (
                                        (lr_results_p1['response_functions'][key]
                                         -
                                         lr_results_m1['response_functions'][key]
                                         ) / (2.0 * self.delta_h))

        if self.rank == mpi_master():
            valstr = '** Time spent on constructing the analytical gradient for '
            valstr += '{:d} frequencies: {:.6f} sec **'.format(
                n_freqs, tm.time() - loop_start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        return num_polgradient

    def _init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

        # generate integration grid
        if self._dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule)

            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                gs_density = AODensityMatrix([scf_tensors['D_alpha']],
                                             denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        return {
            'molgrid': molgrid,
            'gs_density': gs_density,
            'dft_func_label': dft_func_label,
        }

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
            cur_str2 += '{:.4f} a.u.'.format(self.damping)
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
            cur_str4 += '{:.4f} a.u.'.format(self.delta_h)
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
            cur_str = 'Molecular Grid Level            : {}'.format(grid_level)
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
            atom_info = '** Atom #{0}: {1} **'.format(i + 1, labels[i])
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

    def check_real_or_complex_input(self, lr_results):
        """
        Checks if the input LR results and polgrad settings are
        both real or both complex.

        :param lr_results:
            The results of the linear response calculation.
        """

        response_functions = lr_results.get('response_functions', None)
        keys = list(response_functions.keys())
        is_complex_response = (type(response_functions[keys[0]]) == np.dtype('complex128'))

        if (is_complex_response != self.is_complex):
            error_text = 'Mismatch between LR results and polgrad settings!'
            error_text += 'One is complex, the other is not.'
            raise ValueError(error_text)

    def _freq_sanity_check(self):
        """
       Checks if a zero frequency has been input together with resonance Raman.

       This check is due to convergence/singularity issues in the cphf
       subspace solver for some molecules.
       """

        if self.is_complex:
            try:
                idx0 = self.frequencies.index(0.0)
                warn_msg = 'Zero in frequency list for resonance Raman!\n'
                if len(self.frequencies) == 1:
                    error_msg += 'No other frequencies requested.'
                    erro_msg += 'Will continue with normal Raman.'
                else:
                    self.frequencies.pop(idx0)
                    warn_msg += 'It has been removed from the list.'
                    warn_msg += 'Complex pol. gradient will be calculated for:\n'
                    warn_msg += str(obj.frequencies)
                    self.ostream.print_warning(warn_msg)
            except ValueError:
                pass

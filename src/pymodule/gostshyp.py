#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
from mpi4py import MPI
from pathlib import Path

# from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import mpi_master
# from .veloxchemlib import ThreeCenterOverlapDriver
# from .veloxchemlib import ThreeCenterOverlapGradientDriver
# from .veloxchemlib import ThreeCenterOverlapGeom001Driver
# from .veloxchemlib import ThreeCenterOverlapGeom100Driver
# from .veloxchemlib import ThreeCenterOverlapGradientGeom100Driver
# from .veloxchemlib import ThreeCenterOverlapGradientGeom001Driver
# from .veloxchemlib import ThreeCenterR2Driver
# from .veloxchemlib import ThreeCenterRR2Driver
from .veloxchemlib import compute_tco_s_fock
from .veloxchemlib import compute_screened_tco_s_fock
from .veloxchemlib import compute_tco_s_values
from .veloxchemlib import compute_screened_tco_s_values
from .veloxchemlib import compute_tco_p_fock
from .veloxchemlib import compute_screened_tco_p_fock
from .veloxchemlib import compute_tco_p_values
from .veloxchemlib import compute_screened_tco_p_values
from .veloxchemlib import compute_tco_d_values
from .veloxchemlib import compute_screened_tco_d_values
from .veloxchemlib import compute_tco_f_values
from .veloxchemlib import compute_screened_tco_f_values
from .veloxchemlib import compute_tco_s_gradient
from .veloxchemlib import compute_screened_tco_s_gradient
from .veloxchemlib import compute_tco_p_gradient
from .veloxchemlib import compute_screened_tco_p_gradient
from .outputstream import OutputStream
from .tessellation import TessellationDriver
from .inputparser import (parse_input, print_keywords)
from .errorhandler import assert_msg_critical

class GostshypDriver:
    """
    Implements the GOSTSHYP method for applying hydrostatic pressure to a
    molecular system.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param pressure:
        The applied hydrostatic pressure in GPa.
    :param num_leb_points:
        The number of Lebedev grid points per van der Waals sphere.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - molecule: The molecule.
        - basis: The AO basis set.
        - pressure: The applied hydrostatic pressure.
        - pressure_units: The units of the applied pressure.
        - num_tes_points: The number of points on the tessellated surface.
        - tessellation: The tessellated surface object.
        - comm: The MPI communicator.
        - ostream: The output stream.
    """

    def __init__(self, molecule, basis, pressure, pressure_units, 
        tco_tol, comm=None, ostream=None):
        """
        Initializes the GOSTSHYP method for applying hydrostatic pressure.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.molecule = molecule
        self.basis = basis

        # GOSTSHYP setup
        self.pressure = pressure
        self.pressure_units = pressure_units
        self.num_tes_points = 0
        self.tessellation = None
        self._neg_p_amp = 0
        self.tco_tol = tco_tol
        self._tco_p_stats_printed = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # input keywords
        self._input_keywords = {
            'method_settings': {
                'pressure': ('float', 'applied pressure (default in  MPa)'),
                'pressure_units': ('str', 'the units of the given pressure'),
                'tco_tol': ('float', 'screening threshold for three-center overlap integrals')
            }
        }

    def _print_tco_p_screening_stats(self, screened_counts, total_counts):
        """
        Prints the primitive-pair screening statistics for the screened
        three-center p-type overlap integrals. Since the screening decision
        does not depend on the density matrix, the statistics are the same
        for every SCF iteration on a given geometry/tessellation, so this is
        only printed once per driver instance.

        :param screened_counts:
            The number of screened primitive pairs per angular-momentum block.
        :param total_counts:
            The total number of primitive pairs per angular-momentum block.
        """

        if self._tco_p_stats_printed:
            return

        labels = ['S-S', 'S-P', 'S-D', 'S-F', 'P-P', 'P-D', 'P-F', 'D-D', 'D-F', 'F-F']

        self.ostream.print_info(
            'Screened three-center overlap integrals: primitive-pair screening statistics (screened / total):')
        for label, screened, total in zip(labels, screened_counts, total_counts):
            if total > 0:
                self.ostream.print_info(f'  {label}: {screened} / {total}')
        self.ostream.flush()

        self._tco_p_stats_printed = True

    def screened_gostshyp_contrib(self, den_mat, tessellation_settings=None):
        """
        Computes contributions to the total energy and Fock matrix from
        GOSTSHYP method.

        :param den_mat:
            The density matrix.
        :param tessellation_settings:
            The dictionary of tessellation settings

        :return:
            The GOSTSHYP contribution to energy and Fock matrix.
        """

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j,
        # normalization constants, N_j
        initial_exponents = np.pi * np.log(2.0) / self.tessellation[3]
        initial_centers = (self.tessellation[:3].T).copy()
        initial_norms = (self.tessellation[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        #compute f_tilde vector
        f_tilde, screened_counts, total_counts = compute_screened_tco_p_values(self.molecule,
                                                self.basis,
                                                initial_centers,
                                                initial_exponents,
                                                np.full((self.num_tes_points), 1.0),
                                                initial_norms,
                                                initial_norm_consts,
                                                den_mat,
                                                self.tco_tol)
        self._print_tco_p_screening_stats(screened_counts, total_counts)

        # compute amplitudes
        initial_amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0
        np.savetxt('amps_mask.txt', amplitudes_mask, fmt="%5i") #shouldn't go into main
        
        # compute number of grid points associated with negative amplitudes
        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)

        # update gaussian parameters by removing grid points associated 
        # with negative amplitudes
        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        norm_consts = (initial_norm_consts[amplitudes_mask]).copy()
        
        # compute energy contribution
        p_times_g_tilde = compute_screened_tco_s_values(self.molecule,
                                                        self.basis,
                                                        centers,
                                                        exponents,
                                                        amplitudes,
                                                        norm_consts,
                                                        den_mat,
                                                        self.tco_tol)
        
        e_pr = np.sum(p_times_g_tilde)

        # compute Fock matrix contribution
        V1_pr = compute_screened_tco_s_fock(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            amplitudes,
                                            norm_consts,
                                            self.tco_tol)
        
        pre_fac_V2_pr = (p_times_g_tilde / f_tilde[amplitudes_mask])
        
        V2_pr = compute_screened_tco_p_fock(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            pre_fac_V2_pr,
                                            norms,
                                            norm_consts,
                                            self.tco_tol)
        
        V_pr = V1_pr - V2_pr
        
        return e_pr, V_pr

    def gostshyp_contrib(self, den_mat, tessellation_settings=None):
        """
        Computes contributions to the total energy and Fock matrix from
        GOSTSHYP method.

        :param den_mat:
            The density matrix.
        :param tessellation_settings:
            The dictionary of tessellation settings

        :return:
            The GOSTSHYP contribution to energy and Fock matrix.
        """

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j
        initial_exponents = np.pi * np.log(2.0) / self.tessellation[3]
        initial_centers = (self.tessellation[:3].T).copy()
        initial_norms = (self.tessellation[8:11].T).copy()

        #compute f_tilde vector
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       initial_centers,
                                       initial_exponents, 
                                       np.full((self.num_tes_points), 1.0), 
                                       initial_norms, 
                                       den_mat)
        
        # compute amplitudes
        initial_amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0
        np.savetxt('amps_mask.txt', amplitudes_mask, fmt="%5i") #shouldn't go into main
        
        # compute number of grid points associated with negative amplitudes
        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)

        # update gaussian parameters by removing grid points associated 
        # with negative amplitudes
        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        
        # compute energy contribution
        p_times_g_tilde = compute_tco_s_values(self.molecule, 
                                               self.basis, 
                                               centers, 
                                               exponents, 
                                               amplitudes, 
                                               den_mat)
        
        e_pr = np.sum(p_times_g_tilde)

        # compute Fock matrix contribution
        V1_pr = compute_tco_s_fock(self.molecule, 
                                   self.basis, 
                                   centers, 
                                   exponents, 
                                   amplitudes)
        
        pre_fac_V2_pr = (p_times_g_tilde / f_tilde[amplitudes_mask])
        
        V2_pr = compute_tco_p_fock(self.molecule, 
                                   self.basis, 
                                   centers,
                                   exponents, 
                                   pre_fac_V2_pr,
                                   norms)
        
        V_pr = V1_pr - V2_pr
        
        return e_pr, V_pr
    
    def screened_gostshyp_resp_contrib(self, gs_den_mat, trans_den_mat, tessellation_settings=None):
        """
        Computes linear response contributions as second energy derivative
        wrt to density matrix elements. Can also be considered as first 
        derivative of Fock matrix wrt to density matrix element.

        :param gs_den_mat:
            The ground state density matrix.
        :param trans_den_mat:
            The transition (perturbed) density matrix.
        :param tessellation_settings:
            The dictionary of tessellation settings

        :return:
            The GOSTSHYP response contribution.
        """

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j
        # normalization constants, N_j
        initial_exponents = np.pi * np.log(2.0) / self.tessellation[3]
        initial_centers = (self.tessellation[:3].T).copy()
        initial_norms = (self.tessellation[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        # compute f_tilde vector
        f_tilde, screened_counts, total_counts = compute_screened_tco_p_values(self.molecule,
                                                self.basis,
                                                initial_centers,
                                                initial_exponents,
                                                np.full((self.num_tes_points), 1.0),
                                                initial_norms,
                                                initial_norm_consts,
                                                gs_den_mat,
                                                self.tco_tol)
        self._print_tco_p_screening_stats(screened_counts, total_counts)

        # compute amplitudes and remove grid points associated with
        # negative amplitudes
        initial_amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0
        
        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)
        num_points = np.sum(amplitudes_mask)

        # save grid information and auxiliary f_tilde vector 
        # for kept grid points

        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        norm_consts = (initial_norm_consts[amplitudes_mask]).copy()
        f_tilde = (f_tilde[amplitudes_mask]).copy()

        # compute f_tilde at perturbed density (hence prime)
        f_tilde_prime, screened_counts, total_counts = compute_screened_tco_p_values(self.molecule,
                                                      self.basis,
                                                      centers,
                                                      exponents,
                                                      np.full((num_points), 1.0),
                                                      norms,
                                                      norm_consts,
                                                      trans_den_mat,
                                                      self.tco_tol)
        self._print_tco_p_screening_stats(screened_counts, total_counts)

        # compute g_tilde at gs density
        g_tilde = compute_screened_tco_s_values(self.molecule,
                                                self.basis,
                                                centers,
                                                exponents,
                                                np.full((num_points), 1.0),
                                                norm_consts,
                                                gs_den_mat,
                                                self.tco_tol)

        # compute g_tilde at perturbed density (hence prime)
        g_tilde_prime = compute_screened_tco_s_values(self.molecule,
                                                      self.basis,
                                                      centers,
                                                      exponents,
                                                      np.full((num_points), 1.0),
                                                      norm_consts,
                                                      trans_den_mat,
                                                      self.tco_tol)
        
        # compute s-type TCO contribution with prefactor
        prefac_g_mat = - amplitudes * f_tilde_prime / f_tilde

        g_mat_contrib = compute_screened_tco_s_fock(self.molecule,
                                                    self.basis,
                                                    centers,
                                                    exponents,
                                                    prefac_g_mat,
                                                    norm_consts,
                                                    self.tco_tol)
        
        # compute p-type TCO contribution with prefactor
        prefac_f_mat = (amplitudes / f_tilde * 
                        (2 * g_tilde * f_tilde_prime / f_tilde - g_tilde_prime))
        
        f_mat_contrib = compute_screened_tco_p_fock(self.molecule,
                                                    self.basis,
                                                    centers,
                                                    exponents,
                                                    prefac_f_mat,
                                                    norms,
                                                    norm_consts,
                                                    self.tco_tol)
        
        resp_contrib = g_mat_contrib + f_mat_contrib
        
        return resp_contrib

    def gostshyp_resp_contrib(self, gs_den_mat, trans_den_mat, tessellation_settings=None):
        """
        Computes linear response contributions as second energy derivative
        wrt to density matrix elements. Can also be considered as first 
        derivative of Fock matrix wrt to density matrix element.

        :param gs_den_mat:
            The ground state density matrix.
        :param trans_den_mat:
            The transition (perturbed) density matrix.
        :param tessellation_settings:
            The dictionary of tessellation settings

        :return:
            The GOSTSHYP response contribution.
        """

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j
        initial_exponents = np.pi * np.log(2.0) / self.tessellation[3]
        initial_centers = (self.tessellation[:3].T).copy()
        initial_norms = (self.tessellation[8:11].T).copy()

        # compute f_tilde vector
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       initial_centers, 
                                       initial_exponents, 
                                       np.full((self.num_tes_points), 1.0), 
                                       initial_norms, 
                                       gs_den_mat)
        
        # compute amplitudes and remove grid points associated with
        # negative amplitudes
        initial_amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0
        
        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)
        num_points = np.sum(amplitudes_mask)

        # save grid information and auxiliary f_tilde vector
        # for kept grid points

        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        f_tilde = (f_tilde[amplitudes_mask]).copy()

        # compute f_tilde at perturbed density (hence prime)
        f_tilde_prime = compute_tco_p_values(self.molecule, 
                                             self.basis, 
                                             centers, 
                                             exponents, 
                                             np.full((num_points), 1.0), 
                                             norms, 
                                             trans_den_mat)
        
        # compute g_tilde at gs density
        g_tilde = compute_tco_s_values(self.molecule,
                                       self.basis,
                                       centers,
                                       exponents,
                                       np.full((num_points), 1.0),
                                       gs_den_mat)

        # compute g_tilde at perturbed density (hence prime)
        g_tilde_prime = compute_tco_s_values(self.molecule,
                                             self.basis,
                                             centers,
                                             exponents,
                                             np.full((num_points), 1.0),
                                             trans_den_mat)
        
        # compute s-type TCO contribution with prefactor
        prefac_g_mat = - amplitudes * f_tilde_prime / f_tilde

        g_mat_contrib = compute_tco_s_fock(self.molecule,
                                           self.basis,
                                           centers,
                                           exponents,
                                           prefac_g_mat)
        
        # compute p-type TCO contribution with prefactor
        prefac_f_mat = (amplitudes / f_tilde * 
                        (2 * g_tilde * f_tilde_prime / f_tilde - g_tilde_prime))
        
        f_mat_contrib = compute_tco_p_fock(self.molecule,
                                           self.basis,
                                           centers,
                                           exponents,
                                           prefac_f_mat,
                                           norms)
        
        resp_contrib = g_mat_contrib + f_mat_contrib
        
        return resp_contrib

    def screened_gostshyp_grad_contrib(self, den_mat, tessellation_settings=None):
        """
        Computes GOSTSHYP contribution to the molecular gradient

        Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761. Eq. (30)
        
        :param den_mat:
            The density matrix
        :param tessellation_settings:
            The dictionary of tessellation settings
            
        :return:
            The GOSTSHYP molecular gradient contribution 
        """

        # tessellation driver needed for area gradients
        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        # tessellation data is generated
        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        #tesserae areas are extracted
        initial_areas = self.tessellation[3].copy()
        initial_centers = (self.tessellation[:3].T).copy()
        initial_exponents = (np.pi * np.log(2.0) / initial_areas)
        initial_norms = (self.tessellation[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        # f_tilde is computed for all tesserae
        f_tilde, screened_counts, total_counts = compute_screened_tco_p_values(self.molecule,
                                                self.basis,
                                                initial_centers,
                                                initial_exponents,
                                                np.full((self.num_tes_points), 1.0),
                                                initial_norms,
                                                initial_norm_consts,
                                                den_mat,
                                                self.tco_tol)
        self._print_tco_p_screening_stats(screened_counts, total_counts)

        # amplitudes are computed
        initial_amplitudes = self.pressure * initial_areas / f_tilde
        amps_mask = initial_amplitudes >= 0.0
        num_pos_amps = np.sum(amps_mask)
        self._neg_p_amp = self.num_tes_points - num_pos_amps

        # gaussian information is extracted
        areas = (initial_areas[amps_mask]).copy()
        centers = (initial_centers[amps_mask]).copy()
        exponents = (initial_exponents[amps_mask]).copy()
        amplitudes = (initial_amplitudes[amps_mask]).copy()
        norms = (initial_norms[amps_mask]).copy()
        norm_consts = (initial_norm_consts[amps_mask]).copy()

        # amplitudes divided by areas to be used as prefactor
        amps_areas_ratio = amplitudes / areas

        # terms proportional to the area derivative
        g_tilde_term = compute_screened_tco_s_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     amps_areas_ratio,
                                                     norm_consts,
                                                     den_mat,
                                                     self.tco_tol)

        d_tilde_term = compute_screened_tco_d_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     (amps_areas_ratio * exponents),
                                                     norm_consts,
                                                     den_mat,
                                                     self.tco_tol)

        e_tilde_term = compute_screened_tco_f_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     (g_tilde_term * exponents / f_tilde[amps_mask]),
                                                     norms,
                                                     norm_consts,
                                                     den_mat,
                                                     self.tco_tol)
        
        # terms of constant gaussian exponents
        # sum of bra and ket side gradients
        g_tilde_grad = np.array(compute_screened_tco_s_gradient(self.molecule,
                                                                self.basis,
                                                                centers,
                                                                exponents,
                                                                amplitudes,
                                                                norm_consts,
                                                                den_mat,
                                                                self.tco_tol))

        # sum of bra and ket side gradients
        f_tilde_grad = np.array(compute_screened_tco_p_gradient(self.molecule,
                                                                self.basis,
                                                                centers,
                                                                exponents,
                                                                (g_tilde_term * areas / f_tilde[amps_mask]),
                                                                norms,
                                                                norm_consts,
                                                                den_mat,
                                                                self.tco_tol))
        
        natoms = self.molecule.number_of_atoms()
        a_grads = np.zeros((num_pos_amps, natoms, 3))
        
        for atom in range(natoms):
            
            # area gradients
            a_grads[:, atom] = tessellation_drv.comp_area_grad(self.molecule,
                                                               self.tessellation[:, amps_mask],
                                                               atom).T
            
            # gaussian center gradient is added using translational invariance
            tess_ids = (self.tessellation[11, amps_mask] == atom)
            
            g_tilde_grad[tess_ids, atom] -= np.sum(g_tilde_grad, axis = 1)[tess_ids]
            f_tilde_grad[tess_ids, atom] -= np.sum(f_tilde_grad, axis = 1)[tess_ids]
        
        # terms proportional to area gradient multiplied by the area gradient
        area_grad_term = np.sum((2 * g_tilde_term - d_tilde_term + e_tilde_term)[:, np.newaxis, np.newaxis] * a_grads, axis = 0)
        
        # collecting terms for the gradient
        grad = area_grad_term + np.sum(g_tilde_grad - f_tilde_grad, axis = 0)
        
        return grad
    
    def gostshyp_grad_contrib(self, den_mat, tessellation_settings=None): 
        """
        Computes GOSTSHYP contribution to the molecular gradient

        Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761. Eq. (30)
        
        :param den_mat:
            The density matrix
        :param tessellation_settings:
            The dictionary of tessellation settings
            
        :return:
            The GOSTSHYP molecular gradient contribution 
        """

        # tessellation driver needed for area gradients
        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        # tessellation data is generated
        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        #tesserae areas are extracted
        initial_areas = self.tessellation[3].copy()
        initial_centers = (self.tessellation[:3].T).copy()
        initial_exponents = (np.pi * np.log(2.0) / initial_areas)
        initial_norms = (self.tessellation[8:11].T).copy()

        # f_tilde is computed for all tesserae
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       initial_centers, 
                                       initial_exponents, 
                                       np.full((self.num_tes_points), 1.0), 
                                       initial_norms, 
                                       den_mat)
        
        # amplitudes are computed
        initial_amplitudes = self.pressure * initial_areas / f_tilde
        amps_mask = initial_amplitudes >= 0.0
        num_pos_amps = np.sum(amps_mask)
        self._neg_p_amp = self.num_tes_points - num_pos_amps

        # gaussian information is extracted
        areas = (initial_areas[amps_mask]).copy()
        centers = (initial_centers[amps_mask]).copy()
        exponents = (initial_exponents[amps_mask]).copy()
        amplitudes = (initial_amplitudes[amps_mask]).copy()
        norms = (initial_norms[amps_mask]).copy()

        # amplitudes divided by areas to be used as prefactor
        amps_areas_ratio = amplitudes / areas

        # terms proportional to the area derivative
        g_tilde_term = compute_tco_s_values(self.molecule, 
                                            self.basis, 
                                            centers, 
                                            exponents, 
                                            amps_areas_ratio, 
                                            den_mat)
    
        d_tilde_term = compute_tco_d_values(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            (amps_areas_ratio * exponents),
                                            den_mat)
        
        e_tilde_term = compute_tco_f_values(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            (g_tilde_term * exponents / f_tilde[amps_mask]),
                                            norms,
                                            den_mat)
        
        # terms of constant gaussian exponents
        # sum of bra and ket side gradients
        g_tilde_grad = np.array(compute_tco_s_gradient(self.molecule,
                                                       self.basis,
                                                       centers,
                                                       exponents,
                                                       amplitudes,
                                                       den_mat))
        
        # sum of bra and ket side gradients
        f_tilde_grad = np.array(compute_tco_p_gradient(self.molecule,
                                                       self.basis,
                                                       centers,
                                                       exponents,
                                                       (g_tilde_term * areas / f_tilde[amps_mask]),
                                                       norms,
                                                       den_mat))
        
        natoms = self.molecule.number_of_atoms()
        a_grads = np.zeros((num_pos_amps, natoms, 3))
        
        for atom in range(natoms):
            
            # area gradients
            a_grads[:, atom] = tessellation_drv.comp_area_grad(self.molecule,
                                                               self.tessellation[:, amps_mask],
                                                               atom).T
            
            # gaussian center gradient is added using translational invariance
            tess_ids = (self.tessellation[11, amps_mask] == atom)
            
            g_tilde_grad[tess_ids, atom] -= np.sum(g_tilde_grad, axis = 1)[tess_ids]
            f_tilde_grad[tess_ids, atom] -= np.sum(f_tilde_grad, axis = 1)[tess_ids]
        
        # terms proportional to area gradient multiplied by the area gradient
        area_grad_term = np.sum((2 * g_tilde_term - d_tilde_term + e_tilde_term)[:, np.newaxis, np.newaxis] * a_grads, axis = 0)
        
        # collecting terms for the gradient
        grad = area_grad_term + np.sum(g_tilde_grad - f_tilde_grad, axis = 0)
        
        return grad

    def generate_tessellation(self, tessellation_settings={}):
        """
        Initiates the surface tessellation using a Lebedev grid.

        :param tessellation_settings:
            The dictionary of method settings for the tessellation.
        :return:
            The coordinates, surface area, normal vector coordinates and
            reference atoms of the grid points.
        """

        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        self.tessellation = tessellation_drv.compute(self.molecule)

        # TODO error message if an empty tessellation is returned!
        self.num_tes_points = self.tessellation.shape[1]

        return self.tessellation
    
def parse_pressure_units(pressure, units):
    """
    Checks the input given for the units of the applied hydrostatic pressure
    and converts it to atomic units.

    :param pressure:
        The applied hydrostatic pressure.
    :param units:
        The unit in which the pressure is given.
    :return:
        The applied pressure in atomic units.
    """

    assert_msg_critical(units.lower() in [
        'pa', 'pascal', 'hpa', 'hectopascal', 'kpa', 'kilopascal', 'bar', 'mpa',
        'megapascal', 'gpa', 'gigapascal', 'atm', 'atmosphere', 'atmospheric',
        'torr', 'au', 'atomic', 'atomic units',
        ],
        'GOSTSHYP: Invalid unit for pressure')

    #print(units, pressure)
    # implement those in the C++ layer:
    hartree_per_cubic_bohr_in_pascal = 2.942101569713e13
    pascal_in_hartree_per_cubic_bohr = 1.0 / 2.942101569713e13
    atm_in_pascal = 1.01325e5
    torr_in_pascal = 1.33322368421e2

    if units.lower() in ['pa', 'pascal']:
        pressure *= pascal_in_hartree_per_cubic_bohr
    elif units.lower() in ['hpa', 'hectopascal']:
        pressure *= 1.0e2 * pascal_in_hartree_per_cubic_bohr
    elif units.lower() in ['kpa', 'kilopascal']:
        pressure *= 1.0e3 * pascal_in_hartree_per_cubic_bohr
    elif units.lower() == 'bar':
        pressure *= 1.0e5 * pascal_in_hartree_per_cubic_bohr
    elif units.lower() in ['mpa', 'megapascal']:
        pressure *= 1.0e6 * pascal_in_hartree_per_cubic_bohr
    elif units.lower() in ['gpa', 'gigapascal']:
        pressure *= 1.0e9 * pascal_in_hartree_per_cubic_bohr
    elif units.lower() in ['atm', 'atmosphere', 'atmospheric']:
        pressure *= atm_in_pascal * pascal_in_hartree_per_cubic_bohr
    elif units.lower() == 'torr':
        pressure *= torr_in_pascal * pascal_in_hartree_per_cubic_bohr

    return pressure
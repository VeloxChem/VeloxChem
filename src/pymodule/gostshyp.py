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

from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import mpi_master
from .veloxchemlib import ThreeCenterOverlapDriver
from .veloxchemlib import ThreeCenterOverlapGradientDriver
from .veloxchemlib import ThreeCenterOverlapGeom001Driver
from .veloxchemlib import ThreeCenterOverlapGeom100Driver
from .veloxchemlib import ThreeCenterOverlapGradientGeom100Driver
from .veloxchemlib import ThreeCenterOverlapGradientGeom001Driver
from .veloxchemlib import ThreeCenterR2Driver
from .veloxchemlib import ThreeCenterRR2Driver
from .veloxchemlib import compute_tco_s_fock
from .veloxchemlib import compute_tco_s_values
from .veloxchemlib import compute_tco_p_fock
from .veloxchemlib import compute_tco_p_values
from .veloxchemlib import compute_tco_d_values
from .veloxchemlib import compute_tco_f_values
from .veloxchemlib import compute_tco_s_gradient
from .veloxchemlib import compute_tco_p_gradient
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
    :param ostream: # Is this really needed?
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

    def __init__(self, molecule, basis, pressure, pressure_units, comm=None,
            ostream=None):
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
            }
        }


    def get_gostshyp_contribution(self, den_mat, tessellation_settings=None):
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

        # width parameters w_j
        width_params = np.pi * np.log(2.0) / self.tessellation[3]

        #compute f_tilde vector
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       (self.tessellation[:3].T).copy(), 
                                       width_params, 
                                       np.full((self.num_tes_points), 1.0), 
                                       (self.tessellation[4:7].T).copy(), 
                                       den_mat)
        
        amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = amplitudes >= 0.0

        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)

        centers = (self.tessellation[:3].T[amplitudes_mask]).copy()
        exponents = width_params[amplitudes_mask]
        amps = amplitudes[amplitudes_mask]
        norms = (self.tessellation[4:7].T[amplitudes_mask]).copy()
        
        p_times_g_tilde = compute_tco_s_values(self.molecule, 
                                               self.basis, 
                                               centers, 
                                               exponents, 
                                               amps, 
                                               den_mat)
        
        e_pr = np.sum(p_times_g_tilde)

        V1_pr = compute_tco_s_fock(self.molecule, 
                                   self.basis, 
                                   centers, 
                                   exponents, 
                                   amps)
        
        pre_fac_V2_pr = (p_times_g_tilde / f_tilde[amplitudes_mask])
        
        V2_pr = compute_tco_p_fock(self.molecule, 
                                   self.basis, 
                                   centers, 
                                   exponents, 
                                   pre_fac_V2_pr,
                                   norms)
        
        V_pr = V1_pr - V2_pr
        
        return e_pr, V_pr
    
    def get_gostshyp_contribution_occ(self, den_mat, tessellation_settings=None):
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

        # width parameters w_j
        width_params = np.pi * np.log(2.0) / self.tessellation[3]

        #compute f_tilde vector
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       (self.tessellation[:3].T).copy(), 
                                       width_params, 
                                       np.full((self.num_tes_points), 1.0), 
                                       (self.tessellation[8:11].T).copy(), 
                                       den_mat)
        
        amplitudes = self.pressure * self.tessellation[3] / f_tilde

        amplitudes_mask = amplitudes >= 0.0

        self._neg_p_amp = self.num_tes_points - np.sum(amplitudes_mask)

        centers = (self.tessellation[:3].T[amplitudes_mask]).copy()
        exponents = width_params[amplitudes_mask]
        amps = amplitudes[amplitudes_mask]
        norms = (self.tessellation[8:11].T[amplitudes_mask]).copy()
        
        p_times_g_tilde = compute_tco_s_values(self.molecule, 
                                               self.basis, 
                                               centers, 
                                               exponents, 
                                               amps, 
                                               den_mat)
        
        e_pr = np.sum(p_times_g_tilde)

        V1_pr = compute_tco_s_fock(self.molecule, 
                                   self.basis, 
                                   centers, 
                                   exponents, 
                                   amps)
        
        pre_fac_V2_pr = (p_times_g_tilde / f_tilde[amplitudes_mask])
        
        V2_pr = compute_tco_p_fock(self.molecule, 
                                   self.basis, 
                                   centers, 
                                   exponents, 
                                   pre_fac_V2_pr,
                                   norms)
        
        V_pr = V1_pr - V2_pr
        
        return e_pr, V_pr
    

    def get_gostshyp_grad_new_ints(self, den_mat, tessellation_settings=None):
        """
        Computes gradient of SCF energy with respect to nuclear coordinates using the GOSTSHYP model

        Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761 
        
        :param den_mat:
            The density matrix
        :param tessellation_settings:
            The dictionary of tessellation settings
            
        :return:
            The SCF energy gradient at the input pressure 
        """

        # tessellation driver needed for area gradients
        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        # tessellation data is generated
        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        #tesserae areas are extracted
        areas = self.tessellation[3]

        # f_tilde is computed for all tesserae
        f_tilde = compute_tco_p_values(self.molecule, 
                                       self.basis, 
                                       (self.tessellation[:3].T).copy(), 
                                       (np.pi * np.log(2.0) / areas), 
                                       np.full((self.num_tes_points), 1.0), 
                                       (self.tessellation[4:7].T).copy(), 
                                       den_mat)
        
        # amplitudes are computed
        amplitudes = self.pressure * areas / f_tilde
        amps_mask = amplitudes >= 0.0
        num_pos_amps = np.sum(amps_mask)
        self._neg_p_amp = self.num_tes_points - num_pos_amps

        # gaussian information is extracted
        centers = (self.tessellation[:3].T[amps_mask]).copy()
        exponents = (np.pi * np.log(2.0) / areas)[amps_mask]
        amps = amplitudes[amps_mask]
        norms = (self.tessellation[4:7].T[amps_mask]).copy()

        # implemention of gradient
        # J. Chem. Theory Comput. 2025, 21, 747-761, eq. (30)

        # amplitudes divided by areas
        amps_div_areas = (amplitudes / areas)[amps_mask]

        # terms proportional to the area derivative
        g_tilde_term = compute_tco_s_values(self.molecule, 
                                            self.basis, 
                                            centers, 
                                            exponents, 
                                            amps_div_areas, 
                                            den_mat)
    
        d_tilde_term = compute_tco_d_values(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            (amps_div_areas * exponents),
                                            den_mat)
        
        e_tilde_term = compute_tco_f_values(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            (g_tilde_term * exponents / f_tilde[amps_mask]),
                                            norms,
                                            den_mat)
        
        # terms of constant gaussian exponent
        # sum of bra and ket side gradients
        g_tilde_grad = np.array(compute_tco_s_gradient(self.molecule,
                                                       self.basis,
                                                       centers,
                                                       exponents,
                                                       amps,
                                                       den_mat))
        
        # sum of bra and ket side gradients
        f_tilde_grad = np.array(compute_tco_p_gradient(self.molecule,
                                                       self.basis,
                                                       centers,
                                                       exponents,
                                                       (g_tilde_term * areas[amps_mask] / f_tilde[amps_mask]),
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
            tess_ids = (self.tessellation[7, amps_mask] == atom)
            
            g_tilde_grad[tess_ids, atom] -= np.sum(g_tilde_grad, axis = 1)[tess_ids]
            f_tilde_grad[tess_ids, atom] -= np.sum(f_tilde_grad, axis = 1)[tess_ids]
        
        # terms proportional to area gradient multiplied by the area gradient
        area_grad_term = np.sum(np.reshape(np.repeat((2 * g_tilde_term - d_tilde_term + e_tilde_term), natoms * 3), (num_pos_amps, natoms, 3)) * a_grads, axis = 0)
        
        # collecting terms for the gradient
        grad = area_grad_term + np.sum(g_tilde_grad - f_tilde_grad, axis = 0)
        
        return grad
    
    # def get_gostshyp_contribution(self, den_mat, tessellation_settings=None):
    #     """
    #     Computes contributions to the total energy and Fock matrix from
    #     GOSTSHYP method.

    #     :param den_mat:
    #         The density matrix.
    #     :param tessellation_settings:
    #         The dictionary of tessellation settings

    #     :return:
    #         The GOSTSHYP contribution to energy and Fock matrix.
    #     """

    #     if self.num_tes_points == 0:
    #         self.generate_tessellation(tessellation_settings)

    #     # set up needed components:

    #     # width parameters w_j
    #     width_params = np.pi * np.log(2.0) / self.tessellation[3]

    #     # initialize three center overlap integral drivers
    #     tco_drv = ThreeCenterOverlapDriver()
    #     tcog_drv = ThreeCenterOverlapGradientDriver()

    #     # initialize energy and Hamiltonian contribution:
    #     # (definition of V_pr unnecessarily complicated)
    #     e_pr, V_pr = 0.0, np.zeros(den_mat.shape)

    #     # timing
    #     g_time, f_time = 0.0, 0.0

    #     # first try: loop over tessellation points j
    #     for j in range(self.num_tes_points):

    #         # define Gaussian functions and normal vector components
    #         exp = [width_params[j]]
    #         center = [self.tessellation[:3, j].tolist()]
    #         norm_vec = self.tessellation[4:7, j]
    #         pre_fac = [1]
            
    #         g_t0 = tm.time()
    #         tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
    #         g_time += (tm.time() - g_t0)

    #         f_t0 = tm.time()
    #         tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
    #         f_time += (tm.time() - f_t0)

    #         g_mat = tco_mat.to_numpy()

    #         x_mat = tcog_mats.matrix('X').to_numpy()
    #         y_mat = tcog_mats.matrix('Y').to_numpy()
    #         z_mat = tcog_mats.matrix('Z').to_numpy()

    #         f_mat =     ( norm_vec[0] * x_mat
    #                     + norm_vec[1] * y_mat
    #                     + norm_vec[2] * z_mat ) * (-1)

    #         f_tilde = np.sum(den_mat * f_mat)

    #         # calculate parameter p_j
    #         p = self.pressure * self.tessellation[3, j] / f_tilde

    #         g_tilde = np.sum(den_mat * g_mat)

    #         e_contrib = p * g_tilde

    #         pos_amp = (e_contrib >= 0.0)

    #         if pos_amp:
    #             e_pr += e_contrib
    #         else:
    #             print('NEGATIVE AMPLITUDE DETECTED')
    #             continue

    #         # components from above enable the caluclation of the contribution
    #         # to the Hamiltonian

    #         V_pr += (p * g_mat - p * g_tilde / f_tilde * f_mat)

    #     return e_pr, V_pr, g_time, f_time

    # def get_vectorized_gostshyp_contribution(self, den_mat, tessellation_settings=None):

    #     if self.num_tes_points == 0:
    #         self.generate_tessellation(tessellation_settings)

    #     e_contribs = np.apply_along_axis(func1d = self.get_single_point_gostshyp_energy_contribution, 
    #                                 axis = 0, 
    #                                 arr = self.tessellation,
    #                                 den_mat = den_mat)
        
    #     self._neg_p_amp = np.sum(np.isnan(e_contribs))
    #     e_pr = np.nansum(e_contribs, axis = 0)
        
    #     fock_pr = np.sum(   np.apply_along_axis(func1d = self.get_single_point_gostshyp_fock_contribution, 
    #                                 axis = 0, 
    #                                 arr = self.tessellation,
    #                                 den_mat = den_mat), 
    #                         axis = 2    )

    #     return e_pr, fock_pr

    # def get_single_point_gostshyp_energy_contribution(self, gaussian_information, den_mat):

    #     """
    #     Computes contribution to the energy from a single tessellation point.

    #     Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761 

    #     :param den_mat:
    #         The density matrix.
    #     :param tessellation_settings:
    #         The dictionary of tessellation settings

    #     :return:
    #         The GOSTSHYP contribution to energy.
    #     """
        
    #     a = gaussian_information[3]

    #     f_tilde = self.comp_f_tilde(gaussian_information, den_mat)
    #     g_tilde = self.comp_g_tilde(gaussian_information, den_mat)

    #     # compute Gaussian amplitude
    #     p_amp = self.pressure * a / f_tilde

    #     # compute energy contribution
    #     e_contrib = p_amp * g_tilde

    #     # TODO: raise the detection of negative amplitudes properly
    #     #       (counter and proper output)
    #     if (p_amp >= 0.0):
    #         return e_contrib
    #     else:
    #         #print('NEGATIVE AMPLITUDE DETECTED')
    #         return np.nan
    
    # def get_single_point_gostshyp_fock_contribution(self, gaussian_information, den_mat):

    #     """
    #     Computes contribution to the Fock matrix from a single tessellation point.

    #     Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761 

    #     :param den_mat:
    #         The density matrix.
    #     :param tessellation_settings:
    #         The dictionary of tessellation settings

    #     :return:
    #         The GOSTSHYP contribution to the Fock matrix.
    #     """
    
    #     a = gaussian_information[3]
    #     width_params = np.pi * np.log(2.0) / gaussian_information[3]
    #     exp = [width_params]
    #     center = [gaussian_information[:3].tolist()]
    #     norm_vec = gaussian_information[4:7]
    #     pre_fac = [1]
        
    #     tco_drv = ThreeCenterOverlapDriver()
    #     tcog_drv = ThreeCenterOverlapGradientDriver()

    #     # compute the three center overlap matrix and three center overlap gradient matrices
    #     tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)

    #     x_mat_grad = tcog_mats.matrix_to_numpy('X')
    #     y_mat_grad = tcog_mats.matrix_to_numpy('Y')
    #     z_mat_grad = tcog_mats.matrix_to_numpy('Z')

    #     # compute dot product between normal vector and 
    #     # three center overlap gradient matrices
    #     f_mat = ( norm_vec[0] * x_mat_grad
    #             + norm_vec[1] * y_mat_grad
    #             + norm_vec[2] * z_mat_grad  ) * (-1)
        
    #     f_tilde = np.einsum('pq,pq->', den_mat, f_mat)
        
    #     p_amp = self.pressure * a / f_tilde

    #     if (p_amp >= 0.0):
    #         tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
    #         g_tilde = np.einsum('pq,pq->', den_mat, tco_mat.to_numpy())
            
    #         V_pr = ( p_amp * tco_mat.to_numpy() 
    #                 - g_tilde * p_amp / f_tilde * f_mat )
            
    #         return V_pr
    #     else:
    #         return np.zeros(den_mat.shape)

    # def get_gostshyp_contribution(self, den_mat, tessellation_settings=None):

    #     """
    #     Computes contributions to the total energy and Fock matrix from
    #     GOSTSHYP method.

    #     Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761 

    #     :param den_mat:
    #         The density matrix.
    #     :param tessellation_settings:
    #         The dictionary of tessellation settings

    #     :return:
    #         The GOSTSHYP contribution to energy and Fock matrix.
    #     """

    #     if self.num_tes_points == 0:
    #         self.generate_tessellation(tessellation_settings)

    #     # initialize energy and Hamiltonian contribution:
    #     e_pr, V_pr = 0.0, np.zeros(den_mat.shape)
    #     neg_p_amp = 0

    #     # first try: loop over tessellation points j
    #     for j in range(self.num_tes_points):
            
    #         gaussian_information = self.tessellation[:, j]

    #         a = gaussian_information[3]
    #         width_params = np.pi * np.log(2.0) / gaussian_information[3]
    #         exp = [width_params]
    #         center = [gaussian_information[:3].tolist()]
    #         norm_vec = gaussian_information[4:7]
    #         pre_fac = [1]
            
    #         f_tilde = self.comp_f_tilde(gaussian_information, den_mat)
    #         g_tilde = self.comp_g_tilde(gaussian_information, den_mat)
            
    #         # compute Gaussian amplitude
    #         p_amp = self.pressure * a / f_tilde
            
    #         # compute energy contribution
    #         e_contrib = p_amp * g_tilde

    #         # TODO: raise the detection of negative amplitudes properly
    #         #       (counter and proper output)
    #         if (p_amp >= 0.0):
    #             e_pr += e_contrib
    #         else:
    #             print('NEGATIVE AMPLITUDE DETECTED')
    #             neg_p_amp += 1
    #             continue

    #         # drivers for g_{mu nu, j}, f_{mu nu, j} matrices
    #         tco_drv = ThreeCenterOverlapDriver()
    #         tcog_drv = ThreeCenterOverlapGradientDriver()
        
    #         # compute the three center overlap matrix and three center overlap gradient matrices
    #         tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
    #         tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
        
    #         x_mat_grad = tcog_mats.matrix_to_numpy('X')
    #         y_mat_grad = tcog_mats.matrix_to_numpy('Y')
    #         z_mat_grad = tcog_mats.matrix_to_numpy('Z')
            
    #         # compute dot product between normal vector and 
    #         # three center overlap gradient matrices
    #         f_mat = ( norm_vec[0] * x_mat_grad
    #                 + norm_vec[1] * y_mat_grad
    #                 + norm_vec[2] * z_mat_grad  ) * (-1)
            
    #         V_pr += ( p_amp * tco_mat.to_numpy() 
    #                 - g_tilde * p_amp / f_tilde * f_mat )
            
    #     return e_pr, V_pr #, neg_p_amp

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

        # TODO remove if-statement (added for testing of gradient with fixed cavity)

        if tessellation_settings['homemade'] == True:
            tess_data_file = Path('.', tessellation_settings['tess_file'])
            self.tessellation = np.genfromtxt(tess_data_file)

        else:
            self.tessellation = tessellation_drv.compute(self.molecule)

        # TODO error message if an empty tessellation is returned!
        self.num_tes_points = self.tessellation.shape[1]

        return self.tessellation

    def comp_const_exp_010(self, gaussian_information):
        """
        Computes the center derivative of the g_{mu nu, j} and f_{mu nu, j} integrals 
        under constant exponents using translational invariance
        
        :param gaussian_information:
            The tessellation data associated with the grid point
            
        :return:
            The center derivative of the gaussian matrix with respect to the x, y, z center coordinates
            The center derivative of the gaussian gradient matrices with respect to the x, y, z center coordinates
        """
        
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        #norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        n_basis = self.basis.get_dimensions_of_basis()

        dG_dx = np.zeros((n_basis, n_basis))
        dG_dy = np.zeros((n_basis, n_basis))
        dG_dz = np.zeros((n_basis, n_basis))
        
        # drivers
        tco_drv_100 = ThreeCenterOverlapGeom100Driver()
        tco_drv_001 = ThreeCenterOverlapGeom001Driver()
        
        dGx_dx = np.zeros((n_basis, n_basis))
        dGy_dx = np.zeros((n_basis, n_basis))
        dGz_dx = np.zeros((n_basis, n_basis))
        dGx_dy = np.zeros((n_basis, n_basis))
        dGy_dy = np.zeros((n_basis, n_basis))
        dGz_dy = np.zeros((n_basis, n_basis))
        dGx_dz = np.zeros((n_basis, n_basis))
        dGy_dz = np.zeros((n_basis, n_basis))
        dGz_dz = np.zeros((n_basis, n_basis))
        
        # drivers
        tcog_drv_100 = ThreeCenterOverlapGradientGeom100Driver()
        tcog_drv_001 = ThreeCenterOverlapGradientGeom001Driver()
        
        for atom in range(self.molecule.number_of_atoms()):

            tco_100_mats = tco_drv_100.compute(self.molecule, self.basis, exp, pre_fac, center, atom)
            tco_001_mats = tco_drv_001.compute(self.molecule, self.basis, exp, pre_fac, center, atom)

            x_geom_100 = tco_100_mats.matrix_to_numpy('X') #bra atom x-grad
            x_geom_001 = tco_001_mats.matrix_to_numpy('X') #ket atom x-grad

            y_geom_100 = tco_100_mats.matrix_to_numpy('Y') #bra atom y-grad
            y_geom_001 = tco_001_mats.matrix_to_numpy('Y') #ket atom y-grad

            z_geom_100 = tco_100_mats.matrix_to_numpy('Z') #bra atom z-grad
            z_geom_001 = tco_001_mats.matrix_to_numpy('Z') #ket atom z-grad
            
            # use translational invariance to compute the gaussian center derivatives
            dG_dx -= x_geom_100 + x_geom_001
            dG_dy -= y_geom_100 + y_geom_001
            dG_dz -= z_geom_100 + z_geom_001
            
            tcog_100_mats = tcog_drv_100.compute(self.molecule, self.basis, exp, pre_fac, center, atom)
            tcog_001_mats = tcog_drv_001.compute(self.molecule, self.basis, exp, pre_fac, center, atom)

            x_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('X_X') #gaussian x-grad, bra atom x-grad
            x_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('X_Y') #gaussian y-grad, bra atom x-grad
            x_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('X_Z') #gaussian z-grad, bra atom x-grad

            x_geom_001_x_grad = tcog_001_mats.matrix_to_numpy('X_X') #gaussian x-grad, ket atom x-grad
            x_geom_001_y_grad = tcog_001_mats.matrix_to_numpy('X_Y') #gaussian y-grad, ket atom x-grad
            x_geom_001_z_grad = tcog_001_mats.matrix_to_numpy('X_Z') #gaussian z-grad, ket atom x-grad

            y_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('Y_X') #gaussian x-grad, bra atom y-grad
            y_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('Y_Y') #gaussian y-grad, bra atom y-grad
            y_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('Y_Z') #gaussian z-grad, bra atom y-grad

            y_geom_001_x_grad = tcog_001_mats.matrix_to_numpy('Y_X') #gaussian x-grad, ket atom y-grad
            y_geom_001_y_grad = tcog_001_mats.matrix_to_numpy('Y_Y') #gaussian y-grad, ket atom y-grad
            y_geom_001_z_grad = tcog_001_mats.matrix_to_numpy('Y_Z') #gaussian z-grad, ket atom y-grad

            z_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('Z_X') #gaussian x-grad, bra atom z-grad
            z_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('Z_Y') #gaussian y-grad, bra atom z-grad
            z_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('Z_Z') #gaussian z-grad, bra atom z-grad

            z_geom_001_x_grad = tcog_001_mats.matrix_to_numpy('Z_X') #gaussian x-grad, ket atom z-grad
            z_geom_001_y_grad = tcog_001_mats.matrix_to_numpy('Z_Y') #gaussian y-grad, ket atom z-grad
            z_geom_001_z_grad = tcog_001_mats.matrix_to_numpy('Z_Z') #gaussian z-grad, ket atom z-grad
            
            # use translational invariance to compute the gaussian gradient center derivatives 
            dGx_dx -= x_geom_100_x_grad + x_geom_001_x_grad
            dGy_dx -= x_geom_100_y_grad + x_geom_001_y_grad
            dGz_dx -= x_geom_100_z_grad + x_geom_001_z_grad

            dGx_dy -= y_geom_100_x_grad + y_geom_001_x_grad
            dGy_dy -= y_geom_100_y_grad + y_geom_001_y_grad
            dGz_dy -= y_geom_100_z_grad + y_geom_001_z_grad

            dGx_dz -= z_geom_100_x_grad + z_geom_001_x_grad
            dGy_dz -= z_geom_100_y_grad + z_geom_001_y_grad
            dGz_dz -= z_geom_100_z_grad + z_geom_001_z_grad
        
        return np.array([dG_dx, dG_dy, dG_dz]), np.array([dGx_dx, dGy_dx, dGz_dx, dGx_dy, dGy_dy, dGz_dy, dGx_dz, dGy_dz, dGz_dz])

    def comp_grad_g_tilde(self, gaussian_information, den_mat, atom):
        
        """
        Computes the bra, ket and center gradient wrt atom coordinates of the three center overlap 
        under constant Gaussian exponents.
        
        :param gaussian_information:
            The tessellation data associated with the grid point
        :param den_mat:
            The density matrix
        :param atom:
            The current nucleus index
        
        :return:
            The gradient of the gaussian overlap matrices with respect to nuclear coordinates of atom:
            sum_{mu nu} frac{g_{mu nu, j}}{d R_{I alpha}} |_{w_j}
        """
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        #norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        n_basis = self.basis.get_dimensions_of_basis()
        
        # add contribution from gaussian center if grid point belongs to current atom
        if atom == gaussian_information[7]:
            const_exp_tco_010, const_exp_tcog_010 = self.comp_const_exp_010(gaussian_information)
        else:
            const_exp_tco_010 = np.zeros((3, n_basis, n_basis))
        
        # driver
        tco_drv_100 = ThreeCenterOverlapGeom100Driver()
        
        # compute bra derivative of gaussian matrix
        tco_mats_100 = tco_drv_100.compute(self.molecule, self.basis, exp, pre_fac, center, atom)
        
        x_mat_100 = tco_mats_100.matrix_to_numpy('X')
        y_mat_100 = tco_mats_100.matrix_to_numpy('Y')
        z_mat_100 = tco_mats_100.matrix_to_numpy('Z')
        
        # bra, ket and tcr2 derivatives
        x_contrib = np.einsum('pq,pq->', x_mat_100 + x_mat_100.T + const_exp_tco_010[0], den_mat)
        y_contrib = np.einsum('pq,pq->', y_mat_100 + y_mat_100.T + const_exp_tco_010[1], den_mat)
        z_contrib = np.einsum('pq,pq->', z_mat_100 + z_mat_100.T + const_exp_tco_010[2], den_mat)
        
        return np.array([x_contrib, y_contrib, z_contrib])

    def comp_grad_f_tilde(self, gaussian_information, den_mat, atom):
        
        """
        Computes the bra, ket and center gradient wrt atom coordinates of the three center overlap gradient 
        under constant Gaussian exponents. Contracts with density matrix.
        
        :param gaussian_information:
            The tessellation data associated with the grid point
        :param den_mat:
            The density matrix
        :param atom:
            The current nucleus index
        
        :return:
            The gradient of the gaussian overlap gradient matrices with respect to nuclear coordinates of atom:
            sum_{mu nu} frac{f_{mu nu, j}}{d R_{I alpha}} |_{w_j}
        """
        
        # driver
        tcog_100_drv = ThreeCenterOverlapGradientGeom100Driver()
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        n_basis = self.basis.get_dimensions_of_basis()
        
        # add contribution from gaussian center if grid point belongs to current atom
        if atom == gaussian_information[7]:
            const_exp_tco_010, const_exp_tcog_010 = self.comp_const_exp_010(gaussian_information)
        else:
            const_exp_tcog_010 = np.zeros((9, n_basis, n_basis))
        
        # compute bra side gradients of gaussian overlap gradients
        tcog_100_mats = tcog_100_drv.compute(self.molecule, self.basis, exp, pre_fac, center, atom)
        
        x_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('X_X') #gaussian x-grad, bra atom x-grad
        x_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('X_Y') #gaussian y-grad, bra atom x-grad
        x_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('X_Z') #gaussian z-grad, bra atom x-grad
        
        y_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('Y_X') #gaussian x-grad, bra atom y-grad
        y_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('Y_Y') #gaussian y-grad, bra atom y-grad
        y_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('Y_Z') #gaussian z-grad, bra atom y-grad
        
        z_geom_100_x_grad = tcog_100_mats.matrix_to_numpy('Z_X') #gaussian x-grad, bra atom z-grad
        z_geom_100_y_grad = tcog_100_mats.matrix_to_numpy('Z_Y') #gaussian y-grad, bra atom z-grad
        z_geom_100_z_grad = tcog_100_mats.matrix_to_numpy('Z_Z') #gaussian z-grad, bra atom z-grad
        
        x_contrib = (  norm_vec[0] * (x_geom_100_x_grad + x_geom_100_x_grad.T + const_exp_tcog_010[0]) 
                    +  norm_vec[1] * (x_geom_100_y_grad + x_geom_100_y_grad.T + const_exp_tcog_010[1]) 
                    +  norm_vec[2] * (x_geom_100_z_grad + x_geom_100_z_grad.T + const_exp_tcog_010[2])  ) * (-1)
        
        grad_f_tilde_x_contrib = np.einsum('pq,pq->', den_mat, x_contrib)
        
        y_contrib = (  norm_vec[0] * (y_geom_100_x_grad + y_geom_100_x_grad.T + const_exp_tcog_010[3]) 
                    +  norm_vec[1] * (y_geom_100_y_grad + y_geom_100_y_grad.T + const_exp_tcog_010[4]) 
                    +  norm_vec[2] * (y_geom_100_z_grad + y_geom_100_z_grad.T + const_exp_tcog_010[5])  ) * (-1)
        
        grad_f_tilde_y_contrib = np.einsum('pq,pq->', den_mat, y_contrib)
        
        z_contrib = (  norm_vec[0] * (z_geom_100_x_grad + z_geom_100_x_grad.T + const_exp_tcog_010[6]) 
                    +  norm_vec[1] * (z_geom_100_y_grad + z_geom_100_y_grad.T + const_exp_tcog_010[7]) 
                    +  norm_vec[2] * (z_geom_100_z_grad + z_geom_100_z_grad.T + const_exp_tcog_010[8])  ) * (-1)
        
        grad_f_tilde_z_contrib = np.einsum('pq,pq->', den_mat, z_contrib)
        
        return np.array([grad_f_tilde_x_contrib, grad_f_tilde_y_contrib, grad_f_tilde_z_contrib])

    def comp_g_tilde(self, gaussian_information, den_mat):
        """
        Computes the three center overlap matrix
        Contracts it with the density matrix
        
        :param den_mat:
            The density matrix
        :param gaussian_information:
            The tessellation data associated with the grid point
        
        :return:
            three center overlap matrix contracted with density matrix
            sum_{mu nu} D_{mu nu} ( mu | tilde{G}_j | nu )
        """
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        #norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        # driver
        tco_drv = ThreeCenterOverlapDriver()
        
        # compute g_{\mu \nu, j} matrix
        tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
        
        # compute \tilde{g}_j by contracting g_{\mu \nu, j} with density matrix
        g_tilde = np.einsum('pq,pq->', den_mat, tco_mat.to_numpy()) #scalar
        
        return g_tilde

    def comp_f_tilde(self, gaussian_information, den_mat):
        
        """
        Computes the three center overlap gradient matrices, 
        Takes dot product with normal vector
        Contracts it with the density matrix
        
        :param gaussian_information:
            The tessellation data associated with the grid point
        :param den_mat:
            The density matrix
        
        :return:
            sum_{mu nu} D_{mu nu} f_{mu nu, j}
            f_{mu nu, j} = -2 w_j [n_jx ( mu | (x - x_j) tilde{G}_j | nu ) 
                                 + n_jy ( mu | (y - y_j) tilde{G}_j | nu ) 
                                 + n_jz ( mu | (z - z_j) tilde{G}_j | nu )]
        """
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        # driver
        tcog_drv = ThreeCenterOverlapGradientDriver()
        
        # compute gaussian cartesian gradients
        tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
        
        x_mat_grad = tcog_mats.matrix_to_numpy('X')
        y_mat_grad = tcog_mats.matrix_to_numpy('Y')
        z_mat_grad = tcog_mats.matrix_to_numpy('Z')
        
        # compute f_{\mu \nu, j} matrix
        f_mat = ( norm_vec[0] * x_mat_grad
                + norm_vec[1] * y_mat_grad
                + norm_vec[2] * z_mat_grad  ) * (-1)
        
        # compute f_tilde by contracting f_{\mu \nu, j} with density matrix
        f_tilde = np.einsum('pq,pq->', den_mat, f_mat)
        
        return f_tilde

    def comp_d_tilde(self, gaussian_information, den_mat):
        """
        Computes three center overlap with s-type Cartesian d-function 
        Contracts it with density matrix

        :param gaussian_information:
            The tessellation data associated with the grid point
        :param den_mat:
            The density matrix
        
        :return:
            sum D_{mu nu} d_{mu nu, j}
            d_{mu nu, j} = -( mu | (r - r_j)^2 tilde{G}_j | nu ) 
        """
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        #norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        # driver
        tcr2_drv = ThreeCenterR2Driver()
        
        # compute the d_{\mu \nu, j} matrix
        tcr2_mat = tcr2_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
        
        # compute \tilde{d}_j by contracting d_{\mu \nu, j} with density matrix
        d_tilde = np.einsum('pq,pq->', den_mat, tcr2_mat.to_numpy() * (-1)) #scalar
        
        return d_tilde

    def comp_e_tilde(self, gaussian_information, den_mat):
        """
        Computes linear combination of three center overlaps of p-type Cartesian f-functions
        Contracts it with density matrix

        :param gaussian_information:
            The tessellation data associated with the grid point
        :param den_mat:
            The density matrix
        
        :return:
            sum D_{mu nu} e_{mu nu, j}
            e_{mu nu, j} = 2 w_j [n_jx ( mu | (x - x_j) (r - r_j)^2 tilde{G}_j | nu ) 
                                + n_jy ( mu | (y - y_j) (r - r_j)^2 tilde{G}_j | nu ) 
                                + n_jz ( mu | (z - z_j) (r - r_j)^2 tilde{G}_j | nu )]
        """
        
        # gaussian parameters
        width_params = np.pi * np.log(2.0) / gaussian_information[3]
        exp = [width_params]
        center = [gaussian_information[:3].tolist()]
        norm_vec = gaussian_information[4:7]
        pre_fac = [1]
        
        # driver
        tcrr2_drv = ThreeCenterRR2Driver() 
        
        # compute gaussian cartesian gradients
        tcrr2_mats = tcrr2_drv.compute(self.molecule, self.basis, exp, pre_fac, center)
        
        x_mat = tcrr2_mats.matrix_to_numpy('X')
        y_mat = tcrr2_mats.matrix_to_numpy('Y')
        z_mat = tcrr2_mats.matrix_to_numpy('Z')
        
        # compute f_italic
        e_mat = 2 * width_params * ( norm_vec[0] * x_mat
                                   + norm_vec[1] * y_mat
                                   + norm_vec[2] * z_mat ) #matrix
        
        # compute f_tilde by contracting f_italic with density matrix
        e_tilde = np.einsum('pq,pq->', den_mat, e_mat) #scalar
        
        return e_tilde

    def get_vectorized_gostshyp_grad(self, den_mat, tessellation_settings=None):
    
        grad = np.zeros((self.molecule.number_of_atoms(), 3))

        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)
        
        for atom in range(self.molecule.number_of_atoms()):
            
            a_grads = tessellation_drv.comp_area_grad(self.molecule, self.tessellation, atom)
            
            tessellation_with_a_grad = np.vstack((self.tessellation, a_grads))
            
            grad_contribs = np.apply_along_axis(func1d=self.grad_single_point_contrib, 
                                                    axis=0, 
                                                    arr=tessellation_with_a_grad,  
                                                    den_mat=den_mat, 
                                                    atom=atom)
            
            # TODO negative amplitudes
            self._neg_p_amp = np.sum(np.isnan(grad_contribs)) // 3
            
            grad[atom] = np.nansum(grad_contribs, axis=1)
            
        return grad

    def grad_single_point_contrib(self, gaussian_information, den_mat, atom):
        
        a = gaussian_information[3]
        w = np.pi * np.log(2.0) / gaussian_information[3]

        a_grad = gaussian_information[-3:]

        f_tilde = self.comp_f_tilde(gaussian_information, den_mat)
        g_tilde = self.comp_g_tilde(gaussian_information, den_mat)
        d_tilde = self.comp_d_tilde(gaussian_information, den_mat)
        e_tilde = self.comp_e_tilde(gaussian_information, den_mat)

        p_amp = self.pressure * a / f_tilde

        # mark gradient of negative ampltitudes with nan to exclude the point
        if p_amp < 0:
            return np.full((3), np.nan)

        # compute contribution from bra, ket and center derivatives of the unscaled gaussian j (constant exponents)
        grad_g_tilde = self.comp_grad_g_tilde(gaussian_information, den_mat, atom)

        # compute contribution from bra, ket and center derivatives of f part of amplitude (constant exponents)
        grad_f_tilde = self.comp_grad_f_tilde(gaussian_information, den_mat, atom)

        # gradient contribution
        grad_contrib = (  ( 2 * g_tilde * p_amp / a - p_amp * w * d_tilde / a
                        + g_tilde * p_amp * w * e_tilde / (f_tilde * a) ) * a_grad # terms of non-constant exponents
                        + p_amp * (grad_g_tilde - g_tilde / f_tilde * grad_f_tilde   ) # terms of constant exponents
                        )
        
        return grad_contrib

    def get_gostshyp_grad(self, den_mat, tessellation_settings=None):
        """
        Computes gradient of SCF energy with respect to nuclear coordinates using the GOSTSHYP model

        Pausch, Zeller, Neudecker: J. Chem. Theory Comput. 2025, 21, 747-761 
        
        :param den_mat:
            The density matrix
        :param tessellation_settings:
            The dictionary of tessellation settings
            
        :return:
            The SCF energy gradient at the input pressure 
        """

        grad = np.zeros((self.molecule.number_of_atoms(), 3))

        # tessellation driver needed for area gradients
        tessellation_drv = TessellationDriver(self.comm, self.ostream)
        tessellation_drv.update_settings(tessellation_settings)

        # tessellation data is generated
        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)
        
        for atom in range(self.molecule.number_of_atoms()):
            
            # compute area (cavity) gradients for all tessellation points
            a_grad = tessellation_drv.comp_area_grad(self.molecule, self.tessellation, atom)
        
            for j in range(self.num_tes_points):
                
                gaussian_information = self.tessellation[:, j]
                a = gaussian_information[3]
                w = np.pi * np.log(2.0) / gaussian_information[3]
                
                # compute gaussian amplitudes
                f_tilde = self.comp_f_tilde(gaussian_information, den_mat)
                p_amp = self.pressure * a / f_tilde
                
                if p_amp < 0:
                    print('NEGATIVE AMPLITUDE DETECTED')
                    continue

                # compute remaining contracted quantities
                g_tilde = self.comp_g_tilde(gaussian_information, den_mat)
                d_tilde = self.comp_d_tilde(gaussian_information, den_mat)
                e_tilde = self.comp_e_tilde(gaussian_information, den_mat)
                
                # compute contribution from bra, ket and center derivatives of the unscaled gaussian j (constant exponents)
                grad_g_tilde = self.comp_grad_g_tilde(gaussian_information, den_mat, atom)

                # compute contribution from bra, ket and center derivatives of f part of amplitude (constant exponents)
                grad_f_tilde = self.comp_grad_f_tilde(gaussian_information, den_mat, atom)

                # collect gradient contributions
                grad[atom] += (  ( 2 * g_tilde * p_amp / a - p_amp * w * d_tilde / a
                                + g_tilde * p_amp * w * e_tilde / (f_tilde * a) ) * a_grad[:, j] # terms of non-constant exponents
                                + p_amp * (grad_g_tilde - g_tilde / f_tilde * grad_f_tilde   ) # terms of constant exponents
                            )
                
        return grad

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

        # def get_gostshyp_contribution(self, den_mat, tessellation_settings=None):
    #     """
    #     Computes contributions to the total energy and Fock matrix from
    #     GOSTSHYP method.

    #     :param den_mat:
    #         The density matrix.
    #     :param tessellation_settings:
    #         The dictionary of tessellation settings

    #     :return:
    #         The GOSTSHYP contribution to energy and Fock matrix.
    #     """

    #     if self.num_tes_points == 0:
    #         self.generate_tessellation(tessellation_settings)

    #     # set up needed components:

    #     # width parameters w_j
    #     width_params = np.pi * np.log(2.0) / self.tessellation[3]

    #     # initialize three center overlap integral drivers
    #     tco_drv = ThreeCenterOverlapDriver()
    #     tcog_drv = ThreeCenterOverlapGradientDriver()

    #     # initialize energy and Hamiltonian contribution:
    #     # (definition of V_pr unnecessarily complicated)
    #     e_pr, V_pr = 0.0, np.zeros(den_mat.shape)

    #     # exponents = width_params.tolist()
    #     # centers = self.tessellation[:3, :].T.tolist()
    #     # norm_vecs = self.tessellation[4:7, :].T.tolist()
    #     # pre_factors = np.ones(self.num_tes_points).tolist()

    #     # tcog_mats = tcog_drv.compute(self.molecule, self.basis, exponents, pre_factors, centers)

    #     # x_mat = tcog_mats.matrix('X').to_numpy()
    #     # y_mat = tcog_mats.matrix('Y').to_numpy()
    #     # z_mat = tcog_mats.matrix('Z').to_numpy()

    #     # first try: loop over tessellation points j
    #     for j in range(self.num_tes_points):

    #         # define Gaussian functions and normal vector components
    #         exp = [width_params[j]]
    #         center = [self.tessellation[:3, j].tolist()]
    #         norm_vec = self.tessellation[4:7, j]
    #         pre_fac = [1]

    #         # first naive try to calculate F^tilde (with doing everything in
    #         # numpy and ignoring object functions that may exist (but most
    #         # likely don't)

    #         #tcog_mats = tcog_drv.compute(exp, pre_fac, center, self.basis, self.molecule)

    #         tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)

    #         # x_mat = tcog_mats.get_matrix('x').get_full_matrix().to_numpy()
    #         # y_mat = tcog_mats.get_matrix('y').get_full_matrix().to_numpy()
    #         # z_mat = tcog_mats.get_matrix('z').get_full_matrix().to_numpy()

    #         x_mat = tcog_mats.matrix('X').to_numpy() #matrix
    #         y_mat = tcog_mats.matrix('Y').to_numpy() #matrix
    #         z_mat = tcog_mats.matrix('Z').to_numpy() #matrix

    #         f_italic = (  norm_vec[0] * x_mat
    #                     + norm_vec[1] * y_mat
    #                     + norm_vec[2] * z_mat) #matrix

    #         f_tilde = np.einsum('pq,pq->', den_mat, f_italic) #scalar

    #         # calculate parameter p_j
    #         p = -self.pressure * self.tessellation[3, j] / f_tilde #scalar

    #         # same naive try to calculate the energy contribution
    #         tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)

    #         gauss_mat = p * tco_mat.to_numpy()

    #         amplitude = np.einsum('pq,pq->', den_mat, gauss_mat) #* self.tessellation[-1, j]

    #         pos_amp = (amplitude >= 0.0)

    #         # TODO: raise the detection of negative amplitudes properly
    #         #       (counter and proper output)
    #         if pos_amp:
    #             e_pr += amplitude
    #         else:
    #             print('NEGATIVE AMPLITUDE DETECTED')
    #             continue

    #         # components from above enable the caluclation of the contribution
    #         # to the Hamiltonian

    #         # 'correction' term in two steps:
    #         g_tilde_contr = np.einsum('rs,rs->', den_mat,
    #                                   tco_mat.to_numpy())

    #         correction_term = (f_italic * g_tilde_contr * self.pressure * 
    #                            self.tessellation[3, j] * (1.0 / f_tilde**2))

    #         V_pr += (gauss_mat + correction_term) #* self.tessellation[-1, j]

    #     return e_pr, V_pr
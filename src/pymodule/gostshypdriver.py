#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import sys
from mpi4py import MPI

from .veloxchemlib import mpi_master
from .veloxchemlib import compute_screened_tco_s_fock
from .veloxchemlib import compute_screened_tco_s_values
from .veloxchemlib import compute_screened_tco_p_fock
from .veloxchemlib import compute_screened_tco_p_values
from .veloxchemlib import compute_screened_tco_d_values
from .veloxchemlib import compute_screened_tco_f_values
from .veloxchemlib import compute_screened_contracted_tco_s_gradient
from .veloxchemlib import compute_screened_contracted_tco_p_gradient
from .outputstream import OutputStream
from .tessellation import TessellationDriver
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
        The applied hydrostatic pressure in the units given by pressure_units.
    :param pressure_units:
        The units of the given pressure.
    :param tco_tol:
        The screening threshold for the three-center overlap integrals.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - molecule: The molecule.
        - basis: The AO basis set.
        - pressure_au: The applied hydrostatic pressure in atomic units.
        - num_tes_points: The number of points on the tessellated surface.
        - tessellation: The tessellated surface object.
        - num_neg_amp: The number of tessellation points with negative amplitudes.
        - tco_tol: The screening threshold for the three-center overlap integrals.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
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
        self.pressure_au = self.parse_pressure_units(pressure, pressure_units)
        self.num_tes_points = 0
        self.tessellation = None
        self.num_neg_amp = 0
        self.tco_tol = tco_tol

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

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

        # distribute tesserae over ranks
        ave, rem = divmod(self.num_tes_points, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        tess_local = self.tessellation[:, start:end]
        num_local = tess_local.shape[1]

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j and
        # normalization constants, N_j
        initial_areas = tess_local[3].copy()
        initial_exponents = np.pi * np.log(2.0) / initial_areas
        initial_centers = (tess_local[:3].T).copy()
        initial_norms = (tess_local[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        # set up global normalization constant for screening
        initial_global_norm_const = (np.log(2.0) / self.tessellation[3]) ** 1.5
        initial_global_norm_const_max = np.max(initial_global_norm_const)

        # compute f_tilde vector
        initial_f_tilde = compute_screened_tco_p_values(self.molecule,
                                                        self.basis,
                                                        initial_centers,
                                                        initial_exponents,
                                                        np.full((num_local), 1.0),
                                                        initial_norms,
                                                        initial_norm_consts,
                                                        initial_global_norm_const_max,
                                                        den_mat,
                                                        self.tco_tol)

        # compute amplitudes
        initial_amplitudes = self.pressure_au * initial_areas / initial_f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0

        # compute number of grid points associated with negative amplitudes
        local_neg_p_amp = np.sum(~amplitudes_mask)

        # update gaussian parameters by removing grid points associated
        # with negative amplitudes
        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        norm_consts = (initial_norm_consts[amplitudes_mask]).copy()
        f_tilde = (initial_f_tilde[amplitudes_mask]).copy()

        # compute global maximum of normalization constants within kept grid points
        # for screening
        local_norm_const_max = np.max(norm_consts) if norm_consts.size > 0 else -np.inf
        global_norm_const_max = self.comm.allreduce(local_norm_const_max, op=MPI.MAX)

        # compute energy contribution
        p_times_g_tilde = compute_screened_tco_s_values(self.molecule,
                                                        self.basis,
                                                        centers,
                                                        exponents,
                                                        amplitudes,
                                                        norm_consts,
                                                        global_norm_const_max,
                                                        den_mat,
                                                        self.tco_tol)

        local_e_pr = np.sum(p_times_g_tilde)

        # compute Fock matrix contribution
        V1_pr = compute_screened_tco_s_fock(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            amplitudes,
                                            norm_consts,
                                            global_norm_const_max,
                                            self.tco_tol)

        pre_fac_V2_pr = (p_times_g_tilde / f_tilde)

        V2_pr = compute_screened_tco_p_fock(self.molecule,
                                            self.basis,
                                            centers,
                                            exponents,
                                            pre_fac_V2_pr,
                                            norms,
                                            norm_consts,
                                            global_norm_const_max,
                                            self.tco_tol)

        local_V_pr = V1_pr - V2_pr

        e_pr = self.comm.reduce(local_e_pr, root=mpi_master())
        V_pr = self.comm.reduce(local_V_pr, root=mpi_master())

        self.num_neg_amp = self.comm.reduce(local_neg_p_amp, root=mpi_master())

        return e_pr, V_pr

    def gostshyp_resp_contrib(self, gs_den_mat, trans_den_mat, tessellation_settings=None):
        """
        Computes linear response contributions as energy derivative
        wrt to two density matrix elements.

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

        # distribute tesserae over ranks
        ave, rem = divmod(self.num_tes_points, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        # local column-slice of the tessellation: rows are properties,
        # columns are tesserae, so slice the second axis
        tess_local = self.tessellation[:, start:end]
        num_local = tess_local.shape[1]

        # set up needed components:

        # width parameters w_j, gaussian centers r_j, surface normals n_j
        # normalization constants, N_j
        initial_areas = tess_local[3].copy()
        initial_exponents = np.pi * np.log(2.0) / initial_areas
        initial_centers = (tess_local[:3].T).copy()
        initial_norms = (tess_local[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        # compute global maximum of normalization constants for screening
        initial_global_norm_const = (np.log(2.0) / self.tessellation[3]) ** 1.5
        initial_global_norm_const_max = np.max(initial_global_norm_const)

        # compute f_tilde vector
        initial_f_tilde = compute_screened_tco_p_values(self.molecule,
                                                        self.basis,
                                                        initial_centers,
                                                        initial_exponents,
                                                        np.full((num_local), 1.0),
                                                        initial_norms,
                                                        initial_norm_consts,
                                                        initial_global_norm_const_max,
                                                        gs_den_mat,
                                                        self.tco_tol)

        # compute amplitudes and remove grid points associated with
        # negative amplitudes
        initial_amplitudes = self.pressure_au * initial_areas / initial_f_tilde

        amplitudes_mask = initial_amplitudes >= 0.0

        local_neg_p_amp = np.sum(~amplitudes_mask)
        local_num_points = np.sum(amplitudes_mask)

        # save grid information and auxiliary f_tilde vector
        # for kept grid points

        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        norm_consts = (initial_norm_consts[amplitudes_mask]).copy()
        f_tilde = (initial_f_tilde[amplitudes_mask]).copy()

        # compute global maximum of normalization constants within kept grid points for screening
        local_norm_const_max = np.max(norm_consts) if norm_consts.size > 0 else -np.inf
        global_norm_const_max = self.comm.allreduce(local_norm_const_max, op=MPI.MAX)

        # compute f_tilde at perturbed density (hence prime)
        f_tilde_prime = compute_screened_tco_p_values(self.molecule,
                                                      self.basis,
                                                      centers,
                                                      exponents,
                                                      np.full((local_num_points), 1.0),
                                                      norms,
                                                      norm_consts,
                                                      global_norm_const_max,
                                                      trans_den_mat,
                                                      self.tco_tol)

        # compute g_tilde at gs density
        g_tilde = compute_screened_tco_s_values(self.molecule,
                                                self.basis,
                                                centers,
                                                exponents,
                                                np.full((local_num_points), 1.0),
                                                norm_consts,
                                                global_norm_const_max,
                                                gs_den_mat,
                                                self.tco_tol)

        # compute g_tilde at perturbed density (hence prime)
        g_tilde_prime = compute_screened_tco_s_values(self.molecule,
                                                      self.basis,
                                                      centers,
                                                      exponents,
                                                      np.full((local_num_points), 1.0),
                                                      norm_consts,
                                                      global_norm_const_max,
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
                                                    global_norm_const_max,
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
                                                    global_norm_const_max,
                                                    self.tco_tol)

        local_resp_contrib = g_mat_contrib + f_mat_contrib

        resp_contrib = self.comm.reduce(local_resp_contrib, root=mpi_master())

        self.num_neg_amp = self.comm.reduce(local_neg_p_amp, root=mpi_master())

        return resp_contrib

    def gostshyp_grad_contrib(self, den_mat, tessellation_settings=None):
        """
        Computes GOSTSHYP contribution to the molecular gradient.

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

        ave, rem = divmod(self.num_tes_points, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        tess_local = self.tessellation[:, start:end]
        num_local = tess_local.shape[1]

        # tesserae areas are extracted
        initial_areas = tess_local[3].copy()
        initial_centers = (tess_local[:3].T).copy()
        initial_exponents = (np.pi * np.log(2.0) / initial_areas)
        initial_norms = (tess_local[8:11].T).copy()
        initial_norm_consts = (initial_exponents / np.pi) ** (1.5)

        initial_global_norm_const_max = np.max(initial_norm_consts)

        # set up global normalization constant for screening
        initial_global_norm_const = (np.log(2.0) / self.tessellation[3]) ** 1.5
        initial_global_norm_const_max = np.max(initial_global_norm_const)

        # f_tilde is computed for all tesserae
        initial_f_tilde = compute_screened_tco_p_values(self.molecule,
                                                        self.basis,
                                                        initial_centers,
                                                        initial_exponents,
                                                        np.full((num_local), 1.0),
                                                        initial_norms,
                                                        initial_norm_consts,
                                                        initial_global_norm_const_max,
                                                        den_mat,
                                                        self.tco_tol)

        # compute amplitudes
        initial_amplitudes = self.pressure_au * initial_areas / initial_f_tilde
        amplitudes_mask = initial_amplitudes >= 0.0
        local_neg_p_amp = np.sum(~amplitudes_mask)

        # gaussian information is extracted
        areas = (initial_areas[amplitudes_mask]).copy()
        centers = (initial_centers[amplitudes_mask]).copy()
        exponents = (initial_exponents[amplitudes_mask]).copy()
        amplitudes = (initial_amplitudes[amplitudes_mask]).copy()
        norms = (initial_norms[amplitudes_mask]).copy()
        norm_consts = (initial_norm_consts[amplitudes_mask]).copy()
        f_tilde = (initial_f_tilde[amplitudes_mask]).copy()
        tess_masked = tess_local[:, amplitudes_mask].copy()

        # compute global maximum of normalization constants within kept tesserae for screening
        local_norm_const_max = np.max(norm_consts) if norm_consts.size > 0 else -np.inf
        global_norm_const_max = self.comm.allreduce(local_norm_const_max, op=MPI.MAX)

        # amplitudes divided by areas to be used as prefactor
        amps_areas_ratio = amplitudes / areas

        # terms proportional to the area derivative
        g_tilde_term = compute_screened_tco_s_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     amps_areas_ratio,
                                                     norm_consts,
                                                     global_norm_const_max,
                                                     den_mat,
                                                     self.tco_tol)

        d_tilde_term = compute_screened_tco_d_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     (amps_areas_ratio * exponents),
                                                     norm_consts,
                                                     global_norm_const_max,
                                                     den_mat,
                                                     self.tco_tol)

        e_tilde_term = compute_screened_tco_f_values(self.molecule,
                                                     self.basis,
                                                     centers,
                                                     exponents,
                                                     (g_tilde_term * exponents / f_tilde),
                                                     norms,
                                                     norm_consts,
                                                     global_norm_const_max,
                                                     den_mat,
                                                     self.tco_tol)

        # list of corresponding parent atom indices for each tessera
        parent_atom_ids = np.ascontiguousarray(tess_local[11, amplitudes_mask].astype(np.int32))

        # terms of constant gaussian exponents
        # bra, ket and center gradients contracted over tessera axis
        g_tilde_grad = compute_screened_contracted_tco_s_gradient(self.molecule,
                                                                  self.basis,
                                                                  centers,
                                                                  exponents,
                                                                  amplitudes,
                                                                  norm_consts,
                                                                  global_norm_const_max,
                                                                  parent_atom_ids,
                                                                  den_mat,
                                                                  self.tco_tol)

        f_tilde_grad = compute_screened_contracted_tco_p_gradient(self.molecule,
                                                                  self.basis,
                                                                  centers,
                                                                  exponents,
                                                                  (g_tilde_term * areas / f_tilde),
                                                                  norms,
                                                                  norm_consts,
                                                                  global_norm_const_max,
                                                                  parent_atom_ids,
                                                                  den_mat,
                                                                  self.tco_tol)

        # area gradients contracted over tessera axis
        area_grad_term = tessellation_drv.compute_area_grad(
            self.molecule, tess_masked, 2 * g_tilde_term - d_tilde_term + e_tilde_term)

        local_grad = area_grad_term + g_tilde_grad - f_tilde_grad

        self.num_neg_amp = self.comm.allreduce(local_neg_p_amp, op=MPI.SUM)

        # the ScfGradientDriver handles the reduce operation for the gradient
        return local_grad

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

        assert_msg_critical(self.tessellation.shape[1] != 0,
                            'GOSTSHYP: No tessellation points generated. Check tessellation settings and/or molecular structure.')

        self.num_tes_points = self.tessellation.shape[1]

        return self.tessellation

    @staticmethod
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
            'torr', 'au', 'atomic', 'atomic units'],
            'GOSTSHYP: Invalid unit for pressure')

        # TODO: implement those in the C++ layer:
        # hartree_per_cubic_bohr_in_pascal = 2.942101569713e13
        pascal_in_hartree_per_cubic_bohr = 1.0 / 2.942101569713e13
        atm_in_pascal = 1.01325e5
        torr_in_pascal = 1.33322368421e2

        if units.lower() in ['pa', 'pascal']:
            pressure_au = pressure * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['hpa', 'hectopascal']:
            pressure_au = pressure * 1.0e2 * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['kpa', 'kilopascal']:
            pressure_au = pressure * 1.0e3 * pascal_in_hartree_per_cubic_bohr
        elif units.lower() == 'bar':
            pressure_au = pressure * 1.0e5 * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['mpa', 'megapascal']:
            pressure_au = pressure * 1.0e6 * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['gpa', 'gigapascal']:
            pressure_au = pressure * 1.0e9 * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['atm', 'atmosphere', 'atmospheric']:
            pressure_au = pressure * atm_in_pascal * pascal_in_hartree_per_cubic_bohr
        elif units.lower() == 'torr':
            pressure_au = pressure * torr_in_pascal * pascal_in_hartree_per_cubic_bohr
        elif units.lower() in ['au', 'atomic', 'atomic units']:
            pressure_au = float(pressure)
        else:
            assert_msg_critical(False, 'GOSTSHYP: Invalid unit for pressure')

        return pressure_au

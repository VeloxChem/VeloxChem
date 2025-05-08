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
                'pressure_units': ('str', 'the units of the given pressre'),
            }
        }

    def get_gostshyp_contribution(self, den_mat, tessellation_settings=None):
        """
        Computes contributions to the total energy and Fock matrix from
        GOSTSHYP method.

        :param den_mat:
            The density matrix.

        :return:
            The GOSTSHYP contribution to energy and Fock matrix.
        """

        if self.num_tes_points == 0:
            self.generate_tessellation(tessellation_settings)

        # set up needed components:

        # width parameters w_j
        width_params = np.pi * np.log(2.0) / self.tessellation[3]

        # initialize three center overlap integral drivers
        tco_drv = ThreeCenterOverlapDriver()
        tcog_drv = ThreeCenterOverlapGradientDriver()

        # initialize energy and Hamiltonian contribution:
        # (definition of V_pr unnecessarily complicated)
        e_pr, V_pr = 0.0, np.zeros(den_mat.shape)

        # exponents = width_params.tolist()
        # centers = self.tessellation[:3, :].T.tolist()
        # norm_vecs = self.tessellation[4:7, :].T.tolist()
        # pre_factors = np.ones(self.num_tes_points).tolist()

        # tcog_mats = tcog_drv.compute(self.molecule, self.basis, exponents, pre_factors, centers)

        # x_mat = tcog_mats.matrix('X').to_numpy()
        # y_mat = tcog_mats.matrix('Y').to_numpy()
        # z_mat = tcog_mats.matrix('Z').to_numpy()

        # first try: loop over tessellation points j
        for j in range(self.num_tes_points):

            # define Gaussian functions and normal vector components
            exp = [width_params[j]]
            center = [self.tessellation[:3, j].tolist()]
            norm_vec = self.tessellation[4:7, j]
            pre_fac = [1]

            # first naive try to calculate F^tilde (with doing everything in
            # numpy and ignoring object functions that may exist (but most
            # likely don't)

            #tcog_mats = tcog_drv.compute(exp, pre_fac, center, self.basis, self.molecule)

            tcog_mats = tcog_drv.compute(self.molecule, self.basis, exp, pre_fac, center)

            # x_mat = tcog_mats.get_matrix('x').get_full_matrix().to_numpy()
            # y_mat = tcog_mats.get_matrix('y').get_full_matrix().to_numpy()
            # z_mat = tcog_mats.get_matrix('z').get_full_matrix().to_numpy()

            x_mat = tcog_mats.matrix('X').to_numpy() #matrix
            y_mat = tcog_mats.matrix('Y').to_numpy() #matrix
            z_mat = tcog_mats.matrix('Z').to_numpy() #matrix

            f_italic = (  norm_vec[0] * x_mat
                        + norm_vec[1] * y_mat
                        + norm_vec[2] * z_mat) #matrix

            f_tilde = np.einsum('pq,pq->', den_mat, f_italic) #scalar

            # calculate parameter p_j
            p = -self.pressure * self.tessellation[3, j] / f_tilde #scalar

            # same naive try to calculate the energy contribution
            tco_mat = tco_drv.compute(self.molecule, self.basis, exp, pre_fac, center)

            gauss_mat = p * tco_mat.to_numpy()

            amplitude = np.einsum('pq,pq->', den_mat, gauss_mat) #* self.tessellation[-1, j]

            pos_amp = (amplitude >= 0.0)

            # TODO: raise the detection of negative amplitudes properly
            #       (counter and proper output)
            if pos_amp:
                e_pr += amplitude
            else:
                print('NEGATIVE AMPLITUDE DETECTED')
                continue

            # components from above enable the caluclation of the contribution
            # to the Hamiltonian

            # 'correction' term in two steps:
            g_tilde_contr = np.einsum('rs,rs->', den_mat,
                                      tco_mat.to_numpy())

            correction_term = (f_italic * g_tilde_contr * self.pressure * 
                               self.tessellation[3, j] * (1.0 / f_tilde**2))

            V_pr += (gauss_mat + correction_term) #* self.tessellation[-1, j]

        return e_pr, V_pr

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

        # TODO remove if-statement (added for testing purposes)

        if tessellation_settings['homemade'] == True:

            tess_data_file = Path('.', tessellation_settings['tess_file'])
            self.tessellation = np.genfromtxt(tess_data_file)

        else:

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

    print(units, pressure)
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

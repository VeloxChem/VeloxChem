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
import math
import sys
from mpi4py import MPI

from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom
from .outputstream import OutputStream
from .inputparser import (parse_input, print_keywords,
                          get_random_string_parallel)

class TessellationDriver:
    """
    Implements the discretization of the van der Waals surface
    into a grid with octahedral symmetry.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - num_leb_points: The number of Lebedev points per van der Waals sphere.
        - tssf: The tessellation sphere scaling factor.
        - discretization: The surface discretization method.
        - switching_thresh: The (I)SWIG switching function threshold.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the surface tessellation to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # Lebedev grid setup
        self.num_leb_points = 110
        self.tssf = 1.2
        self.discretization = 'fixed'
        self.switching_thresh = 1.0e-8

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # filename for file writing
        self.filename = get_random_string_parallel(self.comm)

        # mpi information

        # input keywords
        self._input_keywords = {
            'method_settings': {
                'num_leb_points': ('int', 'number of points per sphere'),
                'tssf': ('float', 'tessellation sphere scaling factor'),
                'discretization': ('str', 'surface discretization method'),
                'switching_thresh': ('float', 'switching function threshold'),
                'filename': ('str', 'filename for writing the tessellation'),
            }
        }

    def print_keywords(self):
        """
        Prints input keywords in tessellation driver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def update_settings(self, method_dict=None):
        """
        Updates the settings in tessellation driver.

        :param method_dict:
            The dictionary of method settings.
        """

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        if method_dict is not None:
            self.method_dict = dict(method_dict)

        parse_input(self, method_keywords, method_dict)

    def compute(self, molecule):
        """
        Solves for the van der Waals surface discretized with a Lebedev grid.

        :param molecule:
            The molecule.

        :return:
            The coordinates, surface area, normal vector coordinates and
            reference atoms of the grid points.
        """

        # check wether the requested number of points is valid
        self.update_num_points()

        unit_sphere = self.generate_lebedev_grid()

        norm_vecs = self.get_norm_vecs(unit_sphere)

        vdw_radii = molecule.vdw_radii_to_numpy() * self.tssf

        # create a scaled sphere for every unique atom type
        vdw_spheres = {
            i: self.scale_sphere(unit_sphere, i) for i in np.unique(vdw_radii)
        }

        unscaled_w = unit_sphere[:, 3] * 4.0 * np.pi

        # find neighboring atoms for a more efficient removal of points
        neighbors = self.find_neighbors(molecule, vdw_radii, unscaled_w.max())

        # initialize the van der Waals surface object containing the following
        # information for every grid point: x, y, z, A, n_x, n_y, n_z, atom_num, F
        vdw_surface = np.empty((9,0))

        for i, rad in enumerate(vdw_radii):

            # select the correctly scaled sphere, add the corresponding normal
            # vector components  and mark it with the current atom number
            sphere = vdw_spheres[rad].copy()
            sphere = np.vstack((sphere, norm_vecs.copy(),
                                np.full(sphere.shape[1], i)))

            # translate sphere to the center of the respective atom
            sphere[:3, :] += np.array(molecule.get_atom_coordinates(i))[:, None]

            # remove non-contributing points
            contribution_mask, sw_functions = self.get_contribution_mask(sphere,
                i, neighbors, molecule, vdw_radii, unscaled_w)

            contribution = np.vstack((sphere[:, contribution_mask],
                                      sw_functions[contribution_mask]))

            vdw_surface = np.hstack((vdw_surface, contribution))

        # correct the surface areas
        vdw_surface[3, :] *= vdw_surface[-1, :]

        # write surface to file for testing purposes:
        self.write_grid_to_file(vdw_surface)

        # output for testing
        ######################################################
        tot_surf_str = 'The total area of the van der Waals surface is '
        tot_surf_str += str(np.sum(vdw_surface[3, :]) * bohr_in_angstrom()**2)
        tot_surf_str += ' angstrom^2'
        print(tot_surf_str)
        print('with {} tessellation points'.format(vdw_surface.shape[1]))
        ######################################################

        return vdw_surface

    def generate_lebedev_grid(self):
        """
        Generates Lebedev grid on a sphere centered at the origin.

        :return:
            The list of cartesian coordinates and surface areas of the grid
            points.
        """

        # Get grid generated from the C++ class. The spheres are not scaled.
        # Generating from tabulated angles would allow a wider range of
        # different grid points and an on-the-fly scaling.
        leb_grid = gen_lebedev_grid(self.num_leb_points)
        #leb_grid.generate_grid()

        return leb_grid

    def get_norm_vecs(self, grid):
        """
        Gets the normal vector components of grid points centered at the origin.

        :param grid:
            The grid points.

        :return:
            The normal vectors for each grid point.
        """

        return np.vstack((grid[:, 0], grid[:, 1],
                          grid[:, 2]))

    def get_contribution_mask(self, sphere, idx, neighbors, molecule, vdw_radii,
                              weights):
        """
        Masks grid points of a sphere that don't contribute to the van der Waals
        surface.

        :param sphere:
            The grid points of the full sphere.
        :param idx:
            The index of the current atom.
        :param neighbors:
            The dictionary of neighboring atoms.
        :param molecule:
            The molecule.
        :param vdw_radii:
            The scaled van der Waals radii.
        :param weights:
            The Lebedev integration weights.

        :return:
            The mask for the remaining grid points that contribute to the van
            der Waals surface.
        """

        mask = [True] * sphere.shape[1]
        
        sw_funcs = np.ones(sphere.shape[1])

        if self.discretization.lower() == 'fixed':
            for j in neighbors[idx]:
                diff = sphere[:3, :] - np.array(
                                      molecule.get_atom_coordinates(j))[:, None]
                distances = np.linalg.norm(diff, axis=0)
                neighbor_mask = distances > vdw_radii[j]
                mask = [a and b for a, b in zip(mask, neighbor_mask)]

            sw_funcs[~np.array(mask)] = 0.0

        elif self.discretization.lower() in ['swig', 'iswig']:

            r_i = vdw_radii[idx]

            gamma = np.sqrt(14.0 / self.num_leb_points)
            alpha = 0.5 + 1.0 / gamma - np.sqrt(1.0 / gamma**2 - 1.0 / 28.0)

            # TODO: give error for largest grid "combo X grid points and "iswig" option not compatible"
            if self.discretization.lower() == 'iswig':
                zeta = self.get_zeta(self.num_leb_points)

            for j in neighbors[idx]:

                r_j = vdw_radii[j]

                diff = sphere[:3, :] - np.array(
                                      molecule.get_atom_coordinates(j))[:, None]
                distances = np.linalg.norm(diff, axis=0)

                if self.discretization.lower() == 'swig':
                    sw_radius = r_j * gamma
                    inner_j = r_j - alpha * sw_radius

                    elem_sw_funcs = [self.swig_elem_sw_func(
                                  (d - inner_j) / sw_radius) for d in distances]

                elif self.discretization.lower() == 'iswig':
                    zetas = zeta / (r_i * np.sqrt(weights))

                    elem_sw_funcs = [self.iswig_elem_sw_func(
                                   x, z, r_j) for x, z in zip(distances, zetas)]

                sw_funcs *= elem_sw_funcs

            mask = sw_funcs > self.switching_thresh

        return mask, sw_funcs

    def swig_elem_sw_func(self, x):
        """
        Returns the value for the SWIG elementary switching function.

        :param x:
            The function argument.

        :return:
            The function value.
        """

        if x < 0.0:
            return 0.0
        elif x > 1.0:
            return 1.0
        else:
            return (10.0 * x**3 - 15.0 * x**4 + 6.0 * x**5)

    def iswig_elem_sw_func(self, x, zeta, r_j):
        """
        Returns the value for the ISWIG elementary switching function.

        :param x:
            The function argument.
        :param zeta:
            The zeta value.
        :param r_j:
            The radius of atom j.
        
        :return:
            The function value.
        """

        return (1.0 - 0.5 * (math.erf(zeta * (r_j - x)) +
                             math.erf(zeta * (r_j + x))))

    def scale_sphere(self, grid, radius):
        """
        Scales the unit sphere Lebedev grid to a sphere of a given radius.

        :param grid:
            The grid points on the unit sphere.
        :param radius:
            The radius of the scaled sphere.

        :return:
            The grid points of the scaled sphere.
        """

        scaled_sphere = np.zeros((4, self.num_leb_points))

        scaled_sphere[0, :] = grid[:, 0] * radius
        scaled_sphere[1, :] = grid[:, 1] * radius
        scaled_sphere[2, :] = grid[:, 2] * radius

        # make weights sum to the total surface area in bohr^2
        surface_area = 4.0 * np.pi * radius**2
        scaled_sphere[3, :] = grid[:, 3] * surface_area

        return scaled_sphere

    def find_neighbors(self, molecule, radii, weight):
        """
        Creates a dictionary of atom pairs that are considered neighbors with
        respect to the chosen surface discretization method.

        :param molecule:
            The molecule.
        :param radii:
            The list of van der Waals radii.
        :param weight:
            The greatest possible integration weight of the Lebedev grid.

        :return:
            The dictionary of neighboring atoms.
        """

        # coordinates in bohr
        coords = molecule.get_coordinates_in_bohr()

        n_atoms = coords.shape[0]

        gamma = np.sqrt(14.0 / self.num_leb_points)
        alpha = 0.5 + 1.0 / gamma - np.sqrt(1.0 / gamma**2 - 1.0 / 28.0)

        if self.discretization.lower() == 'iswig':
            zeta = self.get_zeta(self.num_leb_points)
        else:
            zeta = 0.0

        neighbors = {n: [] for n in range(n_atoms)}

        for i in range(n_atoms):
            r_i = radii[i]

            # for SWIG
            inner_i = r_i - alpha * r_i * gamma

            # for ISWIG
            zeta_i = zeta * (1.0 / (r_i * np.sqrt(weight)))

            for j in range(i + 1, n_atoms):
                r_j = radii[j]
                d_ij = self.distance(coords[i], coords[j])

                # for SWIG
                inner_j = r_j - alpha * r_j * gamma

                # for ISWIG
                zeta_j = zeta * (1.0 / (r_j * np.sqrt(weight)))


                # TODO: OR statements?
                if (self.discretization.lower() == 'fixed' and
                    d_ij <= (r_i + r_j)):
                    neighbors[i].append(j)
                    neighbors[j].append(i)

                elif (self.discretization.lower() == 'swig' and
                     ((d_ij - r_i - inner_j) * (1.0 / r_j * gamma) < 1.0
                       and (d_ij - r_j - inner_i) * (1.0 / r_i * gamma) < 1.0)):
                    neighbors[i].append(j)
                    neighbors[j].append(i)

                elif (self.discretization.lower() == 'iswig' and
                      (self.iswig_elem_sw_func((d_ij - r_i), zeta_i, r_j) < 1.0
                       and self.iswig_elem_sw_func((d_ij - r_j), zeta_j, r_i) < 1.0)):
                    neighbors[i].append(j)
                    neighbors[j].append(i)

        return neighbors

    def get_zeta(self, zeta):
        """
        Returns the zeta value from the dictionary for different sizes of the
        Lebedev grid.

        :return:
            The zeta value.
        """

        zeta_dict =  {
            50: 4.893,
            110: 4.901,
            194: 4.903,
            302: 4.905,
            434: 4.906,
            590: 4.905,
            770: 4.899,
            974: 4.907,
        }

        return zeta_dict[zeta]

    def distance(self, coords1, coords2):
        """
        Calculates the distance between two points.

        :param coords1:
            The cartesian coordinates of the first point.
        :param coords2:
            The cartesian coordinates of the second point.

        :return:
            The distance between the points.
        """

        return np.linalg.norm(coords1 - coords2)

    def update_num_points(self):
        """
        Checks if requested number of points is valid and updates it to
        the closest valid one if necessary.
        """

        # the available grids in the C++ class:
        num_points_avail = np.array(
                              [6, 50, 110, 194, 302, 434, 590, 770, 974, 2030])

        if self.num_leb_points not in num_points_avail:
            warn_text = '*** Warning: Requested number of '
            warn_text += str(self.num_leb_points)
            warn_text += ' points for the Lebedev grid is invalid.'

            self.ostream.print_header(warn_text.ljust(97))

            warn_text = '***' + ' ' * 10
            warn_text += 'Valid numbers of grid points are: '
            warn_text += '6, 50, 110, 194, 302, 434, 590, 770, 974 and 2030.'

            self.ostream.print_header(warn_text.ljust(97))

            self.num_leb_points = num_points_avail[np.abs(
                              num_points_avail - self.num_leb_points).argmin()]

            warn_text = '***' + ' ' * 10
            warn_text += 'A number of '
            warn_text += str(self.num_leb_points)
            warn_text += ' points is used instead as the closest valid number.'

            self.ostream.print_header(warn_text.ljust(97))
            self.ostream.print_blank()

            self.ostream.flush()

        return self.num_leb_points

    def write_grid_to_file(self, vdw_surface):
        """
        Writes the cartesian coordinates of the surface tessellation points to
        an .xyz file.

        :param vdw_surface:
            The surface tessellation object.
        """
        with open(f"{self.filename}_tessellation.xyz","w") as f:
            f.write(f"{vdw_surface.shape[1]}  \n\n")
            for i in vdw_surface.T:
                f.write(f"  x   {i[0]}  {i[1]}  {i[2]}    \n")

    ### remove later:
    def print_num_points(self):
        """
        Dummy function that prints the number of points included in the Lebedev
        grid.
        """

        cur_str = 'The Lebedev grid contains: '
        cur_str += str(self.num_leb_points)
        cur_str += ' points!!'

        self.ostream.print_header(cur_str)

        self.ostream.print_blank()
        self.ostream.flush()




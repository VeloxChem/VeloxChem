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
import math
import sys
from mpi4py import MPI

try:
    import scipy
except ImportError:
    pass

from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom
from .outputstream import OutputStream
from .inputparser import (parse_input, print_keywords,
                          get_random_string_parallel)
from .errorhandler import assert_msg_critical


class TessellationDriver:
    """
    Implements the discretization of the van der Waals surface
    into a grid with octahedral symmetry including grid point areas geometric derivatives

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - num_lebedev_points: The number of Lebedev points per van der Waals sphere.
        - tssf: The tessellation sphere scaling factor.
        - discretization: The surface discretization method.
        - switching_thresh: The (I)SWIG switching function threshold.
        - r_ext: The extension radius for the outer cavity correction in angstrom.
        - atom: Area gradient is computed with respect to atom's position.
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
        self.num_lebedev_points = 110
        self.tssf = 1.2
        self.discretization = 'fixed'
        self.switching_thresh = 1.0e-8
        self.r_ext = 0.0

        # area gradient setup
        self.atom = None

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # filename for file writing
        self.filename = get_random_string_parallel(self.comm)

        # input keywords
        self._input_keywords = {
            'method_settings': {
                'num_lebedev_points': ('int', 'number of points per sphere'),
                'tssf': ('float', 'tessellation sphere scaling factor'),
                'discretization': ('str', 'surface discretization method'),
                'switching_thresh': ('float', 'switching function threshold'),
                'filename': ('str', 'filename for writing the tessellation'),
                'r_ext': ('float', 'extension radius for outer cavity correction in angstrom')
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
        # information for every grid point: x, y, z, A, x_oc, y_oc, z_oc, A_oc,
        # n_x, n_y, n_z, atom_num, unscaled w, vdw_rad of atom, F_oc
        # NOTE: unscaled w and vdw_rad are used in area gradient

        vdw_surface = np.empty((15,0))

        for i, rad in enumerate(vdw_radii):

            # select the correctly scaled sphere, add the corresponding normal
            # vector components and mark it with the current atom number,
            # unscaled weight and vdw radii
            sphere = vdw_spheres[rad].copy()
            sphere = np.vstack((sphere, norm_vecs.copy(),
                                np.full(sphere.shape[1], i), unscaled_w,
                                np.full(sphere.shape[1], rad)))

            # translate sphere to the center of the respective atom
            sphere[:3, :] += np.array(molecule.get_atom_coordinates(i))[:, None]
            sphere[4:7, :] += np.array(molecule.get_atom_coordinates(i))[:, None]

            # remove non-contributing points
            contribution_mask, sw_functions = self.get_contribution_mask(
                sphere, i, neighbors, molecule, vdw_radii, unscaled_w)

            contribution = np.vstack((sphere[:, contribution_mask],
                                      sw_functions[contribution_mask]))

            vdw_surface = np.hstack((vdw_surface, contribution))

        # correct the surface areas
        vdw_surface[3, :] *= vdw_surface[-1, :]

        return vdw_surface

    def compute_area_grad(self, molecule, tessellation, coeff):

        """Calculate the cavity (area/switching function) gradient.

        Herbert, Lange, J. Chem. Phys. 133, 244111 (2010), Appendix C

        :param molecule             : the molecule.
        :param tessellation         : the tessellation data.
        :param coeff                : the area gradient coefficients for all grid points.

        :return                     : the area gradient contribution (natoms, 3).
        """

        assert_msg_critical(
            self.discretization.lower() != 'fixed',
            'GOSTSHYP: Molecular gradient not available with the fixed discretization scheme. Use SWIG or ISWIG.')

        radii = molecule.vdw_radii_to_numpy() * self.tssf + self.r_ext / bohr_in_angstrom()
        coords = molecule.get_coordinates_in_bohr()
        natoms = molecule.number_of_atoms()
        num_tes_points = tessellation.shape[1]
        areas = tessellation[3]
        c = coeff * areas
        parent_atom_ids = tessellation[11].astype(int)

        # define parameters used in the swig scheme
        if self.discretization.lower() == 'swig':
            gamma = np.sqrt(14.0 / self.num_lebedev_points)
            alpha = 0.5 + 1.0 / gamma - np.sqrt(1.0 / gamma**2 - 1.0 / 28.0)

        # define parameters used in the iswig scheme
        elif self.discretization.lower() == 'iswig':
            iswig_input = tessellation[12:14, :].copy()
            iswig_input[1] += self.r_ext / bohr_in_angstrom()
            zetas = (self.get_zeta(self.num_lebedev_points)
                     / (iswig_input[1] * np.sqrt(iswig_input[0])))

        # initialize array to store the area gradient
        area_grad = np.zeros((natoms, 3))
        T = np.zeros((3, num_tes_points))

        for J in range(natoms):

            # distance vector between tessellation point and nuclei positions
            diff = tessellation[4:7] - coords[J][:, None]      # (3, npoints)
            distances = np.linalg.norm(diff, axis=0)           # (npoints,)
            norm_diff = diff / distances

            if self.discretization.lower() == 'swig':
                sw_radius = radii[J] * gamma
                inner_J = radii[J] - alpha * sw_radius
                x = (distances - inner_J) / sw_radius
                t_J = (self.swig_elem_sw_func_derivative_array(x)
                       / (self.swig_elem_sw_func_array(x) * sw_radius)) * norm_diff

            elif self.discretization.lower() == 'iswig':
                t_J = (self.iswig_elem_sw_func_derivative(distances, zetas, radii[J])
                       / self.iswig_elem_sw_func_array(distances, zetas, radii[J])) * norm_diff

            T += t_J

            area_grad[J] -= np.sum(c * t_J, axis=1)  # summed over tessera axis

        weighted_T = c * T  # (3, npoints)

        # scatter contributions onto parent atom gradient
        for d in range(3):
            area_grad[:, d] += np.bincount(parent_atom_ids,
                                           weights=weighted_T[d],
                                           minlength=natoms)

        return area_grad

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
        lebedev_grid = gen_lebedev_grid(self.num_lebedev_points)

        return lebedev_grid

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

            r_i = vdw_radii[idx] + self.r_ext / bohr_in_angstrom()

            gamma = np.sqrt(14.0 / self.num_lebedev_points)
            alpha = 0.5 + 1.0 / gamma - np.sqrt(1.0 / gamma**2 - 1.0 / 28.0)

            if self.discretization.lower() == 'iswig':
                zeta = self.get_zeta(self.num_lebedev_points)

            for j in neighbors[idx]:

                r_j = vdw_radii[j] + self.r_ext / bohr_in_angstrom()

                diff = sphere[4:7, :] - np.array(
                    molecule.get_atom_coordinates(j))[:, None]
                distances = np.linalg.norm(diff, axis=0)

                if self.discretization.lower() == 'swig':
                    sw_radius = r_j * gamma
                    inner_j = r_j - alpha * sw_radius

                    elem_sw_funcs = [self.swig_elem_sw_func(
                        (d - inner_j) / sw_radius) for d in distances]

                elif self.discretization.lower() == 'iswig':
                    zetas = zeta / (r_i * np.sqrt(weights))

                    elem_sw_funcs = [self.iswig_elem_sw_func(x, z, r_j)
                                     for x, z in zip(distances, zetas)]

                sw_funcs *= elem_sw_funcs

            mask = sw_funcs > self.switching_thresh

        return mask, sw_funcs

    @staticmethod
    def erf_array(array):
        """
        Returns erf of an array.
        """

        if 'scipy' in sys.modules:
            return scipy.special.erf(array)
        else:
            # slow alternative in case scipy is not available
            return np.vectorize(math.erf)(array)

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

    def swig_elem_sw_func_array(self, x):
        """
        Returns the value for the SWIG elementary switching function of an array.

        :param x:
            The function argument.

        :return:
            The function value.
        """

        x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)

        return x**3 * (10.0 - 15.0 * x + 6.0 * x**2)

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

    def iswig_elem_sw_func_array(self, x, zeta, r_j):
        """
        Returns the value for the ISWIG elementary switching function of an array.

        :param x:
            The function argument.
        :param zeta:
            The zeta value.
        :param r_j:
            The radius of atom j.

        :return:
            The function value.
        """

        return (1.0 - 0.5 * (self.erf_array(zeta * (r_j - x)) +
                             self.erf_array(zeta * (r_j + x))))

    def swig_elem_sw_func_derivative(self, x):
        """
        Returns the value for the SWIG elementary switching function derivative wrt. x.

        :param x:
            The function argument.

        :return:
            The function value.
        """

        if x < 0.0:
            return 0.0
        elif x > 1.0:
            return 0.0
        else:
            return (30 * x**4 - 60 * x**3 + 30 * x**2)

    def swig_elem_sw_func_derivative_array(self, x):
        """
        Returns the value for the SWIG elementary switching function derivative wrt. x of an array.

        :param x:
            The function argument.

        :return:
            The function value.
        """

        x = np.asarray(x, dtype=float)

        inside = (x >= 0.0) & (x <= 1.0)

        return np.where(inside, 30.0 * x**2 * (1.0 - x)**2, 0.0)

    def iswig_elem_sw_func_derivative(self, x, zeta, r_j):

        """
        Returns the value for the ISWIG elementary switching function derivative wrt. x.
        Works for both scalar and array inputs.

        :param x:
            The function argument.
        :param zeta:
            The zeta value.
        :param r_j:
            The radius of atom j.

        :return:
            The function value.
        """

        return zeta / np.sqrt(np.pi) * (
            np.exp(-zeta**2 * (r_j - x)**2) -
            np.exp(-zeta**2 * (r_j + x)**2))

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

        scaled_sphere = np.zeros((8, self.num_lebedev_points))

        scaled_sphere[0, :] = grid[:, 0] * radius
        scaled_sphere[1, :] = grid[:, 1] * radius
        scaled_sphere[2, :] = grid[:, 2] * radius

        # make weights sum to the total surface area in bohr^2
        surface_area = 4.0 * np.pi * radius**2
        scaled_sphere[3, :] = grid[:, 3] * surface_area

        scaled_sphere[4, :] = grid[:, 0] * (radius + self.r_ext / bohr_in_angstrom())
        scaled_sphere[5, :] = grid[:, 1] * (radius + self.r_ext / bohr_in_angstrom())
        scaled_sphere[6, :] = grid[:, 2] * (radius + self.r_ext / bohr_in_angstrom())

        surface_area_occ = 4.0 * np.pi * (radius + self.r_ext / bohr_in_angstrom())**2
        scaled_sphere[7, :] = grid[:, 3] * surface_area_occ

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

        gamma = np.sqrt(14.0 / self.num_lebedev_points)
        alpha = 0.5 + 1.0 / gamma - np.sqrt(1.0 / gamma**2 - 1.0 / 28.0)

        if self.discretization.lower() == 'iswig':
            zeta = self.get_zeta(self.num_lebedev_points)
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

        zeta_dict = {
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
        Checks if the requested number of points is valid and updates it to
        the closest valid one if necessary.

        :return:
            The number of Lebedev points that will be used.
        """

        self.num_lebedev_points = validate_num_lebedev_points(
            self.num_lebedev_points, self.discretization, self.ostream)

        return self.num_lebedev_points

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

    def visualize_grid(self, molecule):
        """
        Visualizes grid for surface discretization.

        :param molecule:
            The molecule.
        """

        grid = self.compute(molecule).T

        try:
            import py3Dmol as p3d
        except ImportError:
            raise ImportError('Unable to import py3Dmol.')

        assert_msg_critical(grid.shape[1] == 15,
                            'TessellationDriver.visualize_grid: Invalid grid size')

        num_points = grid.shape[0]

        grid_in_angstrom = grid[:, :3] * bohr_in_angstrom()
        grid_xyz_string = f'{num_points}\n\n'

        for i in range(num_points):
            x, y, z = grid_in_angstrom[i]
            grid_xyz_string += f'He {x} {y} {z}\n'

        v = p3d.view(width=600, height=600)

        v.addModel(molecule.get_xyz_string(), 'xyz')
        v.setStyle({'stick': {}})

        v.addModel(grid_xyz_string, 'xyz')
        v.setStyle({'elem': 'He'},
                   {'sphere': {
                       'radius': 0.05,
                       'color': 'red',
                       'opacity': 0.7
                   }})

        v.zoomTo()
        v.show()

def validate_num_lebedev_points(num_lebedev_points, discretization='fixed',
                                ostream=None):
    """
    Checks that the requested number of Lebedev points is available and
    returns the closest valid alternative if it is not.

    :param num_lebedev_points:
        The number of Lebedev points per sphere.
    :param discretization:
        The surface discretization method.
    :param ostream:
        The output stream used for warnings; no warning is printed if None.

    :return:
        The valid number of Lebedev points.
    """

    def print_warning(header, valid_points, new_value):
        if ostream is None:
            return
        ostream.print_header(header.ljust(97))
        warn_text = '***' + ' ' * 10 + 'Valid numbers of grid points are: '
        warn_text += ', '.join(str(p) for p in valid_points[:-1])
        warn_text += f' and {valid_points[-1]}.'
        ostream.print_header(warn_text.ljust(97))
        warn_text = '***' + ' ' * 10 + f'A number of {new_value} points is '
        warn_text += 'used instead as the closest valid number.'
        ostream.print_header(warn_text.ljust(97))
        ostream.print_blank()
        ostream.flush()

    # the Lebedev grids available in the C++ class
    num_lebedev_points_available = (6, 50, 110, 194, 302, 434, 590, 770, 974, 2030)

    # ISWIG is not parameterized for the smallest and largest grid
    num_lebedev_points_iswig_available = (50, 110, 194, 302, 434, 590, 770, 974)

    num_lebedev_points = int(num_lebedev_points)
    is_iswig = str(discretization).lower() == 'iswig'

    # select the grid list applicable to the discretization method
    if is_iswig:
        available = np.array(num_lebedev_points_iswig_available)
        valid_points = num_lebedev_points_iswig_available
        warning_tail = 'invalid with ISWIG'
    else:
        available = np.array(num_lebedev_points_available)
        valid_points = num_lebedev_points_available
        warning_tail = 'invalid'

    # snap to the closest available grid
    if num_lebedev_points not in available:
        closest = int(available[np.abs(available - num_lebedev_points).argmin()])
        print_warning(
            f'*** Warning: Requested number of {num_lebedev_points} points '
            f'for the Lebedev grid is {warning_tail}.',
            valid_points, closest)
        num_lebedev_points = closest

    return num_lebedev_points

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

from pathlib import Path
import numpy as np
import math
import h5py

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .visualizationdriver import VisualizationDriver
from .cubicgrid import CubicGrid
from .errorhandler import assert_msg_critical


class DensityViewerAlternative:
    """
    An approximate but very quick module to compute a densities on a
    grid. It does this by computing the AOs of each unique atom type in an
    atomic grid. These AOs are then translated as needed to the nearest grid
    point of the actual atomic center to compute on-the-fly the density
    on the full grid.

    Instance variables
        - grid_density: The density of grid points per a.u.
        - grid_margins: Minimal distance in a.u. between the box edge and all
          atomic nuclei
        - atombox_radius: Radius in a.u. of the atomix grid
        - atom_centers: List of atomic indices around which the box is focused
          (default all)
        - mo_threshold: Do not compute AO contribution if corresponding MO
          coefficient is below this value
    """

    def __init__(self):
        """
        Initializes the density viewer.
        """

        # number of grid points per a.u.
        self.grid_density = 4

        # a.u. added around the molecule to create the box
        self.grid_margins = 3.0

        # Radius in a.u. for each atomic box
        self.atombox_radius = 6.0

        # Do not compute AO if the density matrix element is below this value
        self.dm_element_threshold = 1e-4

        # flag for using interpolation when computing densities
        self.use_visualization_driver = False
        self.interpolate = False

        # flag for the type of density
        self.den_type = denmat.rest

        self.density_color_scheme = 'default'
        self.density_isovalue = 0.002
        self.density_opacity = 0.7

        # To focus the grid around only specific atoms (for large systems)
        self.atom_centers = None

        # molecule
        self._molecule = None
        self._basis = None
        self._atomnr = None
        self._coords = None

        # grid
        self._grid = None
        self.origin = None
        self.stepsize = None
        self.npoints = None
        self._atom_origin = None
        self._atom_npoints = None

        # AOs and density
        self._ao_to_atom = None
        self._atom_to_ao = None
        self._ao_dict = None
        self._density_list = None
        self._density_labels = None
        self._is_uhf = False

    def read_hdf5(self, fname):
        """
        Reads the dictionary of densities from a checkpoint file.

        :param fname:
            The checkpoint file name.

        :return:
            A dictionary containing the density matrices as numpy arrays
            and their labels.
        """

        den_dict = {}

        valid_checkpoint = (fname and isinstance(fname, str) and
                            Path(fname).is_file())

        errmsg = f"DensityViewer: Cannot read file {fname}."
        assert_msg_critical(valid_checkpoint, errmsg)

        h5f = h5py.File(fname, "r")

        # TODO: add other density types as they become available
        for key in h5f.keys():
            if "detach" in key or "attach" in key:
                data = np.array(h5f.get(key))
                den_dict[key] = data
            if "hole" in key or "particle" in key:
                data = np.array(h5f.get(key))
                den_dict[key] = data

        return den_dict

    def help_string_widgets_and_display(self):

        return ('Unable to import ipywidgets or IPython.display.\n' +
                'Please install jupyter notebook via pip or conda.')

    def initialize(self, molecule, basis):
        """
        Initializes molecule and cubic grid; pre-computes atomic orbitals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        self._molecule = molecule
        self._basis = basis

        # Define box size
        self._atomnr = np.array(molecule.get_identifiers()) - 1
        self._coords = molecule.get_coordinates_in_bohr()
        if self.atom_centers is None:
            xmin = self._coords[:, 0].min() - self.grid_margins
            xmax = self._coords[:, 0].max() + self.grid_margins
            ymin = self._coords[:, 1].min() - self.grid_margins
            ymax = self._coords[:, 1].max() + self.grid_margins
            zmin = self._coords[:, 2].min() - self.grid_margins
            zmax = self._coords[:, 2].max() + self.grid_margins
        else:
            xmin = self._coords[self.atom_centers, 0].min() - self.grid_margins
            xmax = self._coords[self.atom_centers, 0].max() + self.grid_margins
            ymin = self._coords[self.atom_centers, 1].min() - self.grid_margins
            ymax = self._coords[self.atom_centers, 1].max() + self.grid_margins
            zmin = self._coords[self.atom_centers, 2].min() - self.grid_margins
            zmax = self._coords[self.atom_centers, 2].max() + self.grid_margins

        nx = math.ceil((xmax - xmin) * self.grid_density)
        ny = math.ceil((ymax - ymin) * self.grid_density)
        nz = math.ceil((zmax - zmin) * self.grid_density)
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        dz = (zmax - zmin) / nz

        self.origin = (xmin, ymin, zmin)
        self.stepsize = (dx, dy, dz)
        self.npoints = (nx, ny, nz)
        self._grid = CubicGrid(self.origin, self.stepsize, self.npoints)

        # Create a small box for each atom
        nx_atom = math.ceil(self.atombox_radius / dx * 2)
        ny_atom = math.ceil(self.atombox_radius / dy * 2)
        nz_atom = math.ceil(self.atombox_radius / dz * 2)
        self._atom_origin = (-nx_atom * dx / 2, -ny_atom * dy / 2,
                             -nz_atom * dz / 2)
        self._atom_npoints = (nx_atom, ny_atom, nz_atom)

        atom_grid = CubicGrid(self._atom_origin, self.stepsize,
                              self._atom_npoints)

        # Create atomic map
        vis_drv = VisualizationDriver()
        ao_info = vis_drv.get_atomic_orbital_info(molecule, basis)

        # Create reverse atomic map
        atom_to_ao = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)
        self._atom_to_ao = atom_to_ao
        self._ao_to_atom = [[] for i in range(len(ao_info))]
        for i_atom, atom_orbs in enumerate(atom_to_ao):
            for i_orb, orb in enumerate(atom_orbs):
                self._ao_to_atom[orb] = [i_atom, i_orb]

        # Compute each unique AO on the grid
        self._ao_dict = {}
        for i_atom, atom in enumerate(self._atomnr):
            if atom not in self._ao_dict:
                atomlist = []
                for orb in atom_to_ao[i_atom]:
                    vis_drv.compute_atomic_orbital_for_grid(
                        atom_grid, basis, ao_info[orb])
                    atomlist.append(atom_grid.values_to_numpy())
                self._ao_dict[atom] = atomlist

        self._atom_grid = atom_grid

    def compute_density(self, density_matrix):
        """
        Compute a specific molecular orbital.

        :param density_matrix:
            A numpy array of the density matrix.

        :return:
            A numpy array with the value of the density on the grid.
        """

        if self.use_visualization_driver:
            return self.compute_density_vis_drv(density_matrix)
        else:
            return self.compute_density_simple_loop_atoms(density_matrix)

    def compute_density_vis_drv(self, density_matrix):
        """
        Compute the density matrix on a grid using the
        visualization driver.

        :param density_matrix:
            A numpy array containing the density matrix.

        :return:
            A numpy array with the value of the density on the grid.
        """

        vis_drv = VisualizationDriver()
        dm_ao = AODensityMatrix([density_matrix], self.den_type)
        vis_drv.compute(self._grid, self._molecule, self._basis, dm_ao, 0,
                        'alpha')

        return self._grid.values_to_numpy()

    def compute_density_simple_loop_atoms(self, density_matrix):
        """
        Compute the density by shifting atom grids.
        Loops over atoms and sums all AOs in one shot.

        :param density_matrix:
            A numpy array containing the density matrix.

        :return:
            A numpy array with the value of the density on the grid.
        """

        nx, ny, nz = self._atom_npoints
        natms = self._molecule.number_of_atoms()

        # Initialize the density on the molecular grid
        np_density = np.zeros(self.npoints)

        identifiers = np.array(self._molecule.get_identifiers()) - 1

        ijk_inds = [(i, j, k) for i in [0, 1] for j in [0, 1] for k in [0, 1]]

        # Loop over atoms
        for i_atom in range(natms):
            ao_indices_i = self._atom_to_ao[i_atom]
            atom_id_i = identifiers[i_atom]
            atom_orbs_i = np.array(self._ao_dict[atom_id_i])

            if self.interpolate:
                ti = (self._coords[i_atom] - self.origin +
                      self._atom_origin) / self.stepsize
                ti_floor = np.floor(ti)
                alpha_i = ti - ti_floor

                new_atom_orbs_i = np.zeros(
                    (atom_orbs_i.shape[0], nx + 1, ny + 1, nz + 1))

                for i, j, k in ijk_inds:
                    # coefficient for trilinear interpolation
                    x_coef = (1.0 - alpha_i[0]) if i == 0 else alpha_i[0]
                    y_coef = (1.0 - alpha_i[1]) if j == 0 else alpha_i[1]
                    z_coef = (1.0 - alpha_i[2]) if k == 0 else alpha_i[2]

                    new_atom_orbs_i[
                        :,
                        i:i + nx,
                        j:j + ny,
                        k:k + nz,
                    ] += x_coef * y_coef * z_coef * atom_orbs_i

            else:
                new_atom_orbs_i = atom_orbs_i

            for j_atom in range(natms):
                ao_indices_j = self._atom_to_ao[j_atom]
                atom_id_j = identifiers[j_atom]
                atom_orbs_j = np.array(self._ao_dict[atom_id_j])

                if self.interpolate:
                    tj = (self._coords[j_atom] - self.origin +
                          self._atom_origin) / self.stepsize
                    tj_floor = np.floor(tj)
                    alpha_j = tj - tj_floor

                    new_atom_orbs_j = np.zeros(
                        (atom_orbs_j.shape[0], nx + 1, ny + 1, nz + 1))

                    for i, j, k in ijk_inds:
                        # coefficient for trilinear interpolation
                        x_coef = (1.0 - alpha_j[0]) if i == 0 else alpha_j[0]
                        y_coef = (1.0 - alpha_j[1]) if j == 0 else alpha_j[1]
                        z_coef = (1.0 - alpha_j[2]) if k == 0 else alpha_j[2]

                        new_atom_orbs_j[
                            :,
                            i:i + nx,
                            j:j + ny,
                            k:k + nz,
                        ] += x_coef * y_coef * z_coef * atom_orbs_j

                else:
                    new_atom_orbs_j = atom_orbs_j

                # Translate the atomic grid from the origin to the position of
                # the current atom in the molecular grid and find the indices
                # of the atomic grid origin in the molecular grid.
                if self.interpolate:
                    ti = ti_floor.astype('int')
                    tj = tj_floor.astype('int')
                else:
                    ti = np.round(
                        (self._coords[i_atom] - self.origin + self._atom_origin)
                        / self.stepsize).astype('int')
                    tj = np.round(
                        (self._coords[j_atom] - self.origin + self._atom_origin)
                        / self.stepsize).astype('int')

                # indices of the origin on the atomic grid i
                p1_i = [0, 0, 0]
                # indices of the origin on the atomic grid j
                p1_j = [0, 0, 0]

                # Compare the i and j grids to find the indices
                # of an overlaping grid. t is the position of the origin
                # and no is the number of grid points in the overlap grid.
                t = np.zeros((3), dtype=int)

                if self.interpolate:
                    no = [nx + 1, ny + 1, nz + 1]
                else:
                    no = [nx, ny, nz]

                discard = False

                for x in range(3):
                    if ti[x] < tj[x]:
                        t[x] = tj[x]
                        no[x] += ti[x] - tj[x]
                        p1_i[x] = tj[x] - ti[x]
                    else:
                        t[x] = ti[x]
                        no[x] += tj[x] - ti[x]
                        p1_j[x] = ti[x] - tj[x]
                    # if the atomic grids do not overlap,
                    # discard them.
                    if no[x] <= 0:
                        discard = True
                        break

                # Once we have the overlaping grid, we have to
                # figure out how it looks in relation with the molecular grid
                # and keep track of the indices of the points in the
                # atomic grids i and j.
                for x in range(3):
                    # if the overlaping atomic grid is
                    # completely outside the
                    # molecular grid, discard it.
                    if t[x] >= self.npoints[x]:
                        discard = True
                        break
                    if t[x] + no[x] < 0:
                        discard = True
                        break
                    # if only some points of the grid are outside, readjust
                    if t[x] < 0:
                        p1_i[x] -= t[x]
                        p1_j[x] -= t[x]
                        no[x] += t[x]
                        t[x] = 0
                    if t[x] + no[x] >= self.npoints[x]:
                        no[x] = self.npoints[x] - t[x]

                if discard:
                    continue

                # Calculate the value of the density on the
                # overlaping grid
                for mu in range(new_atom_orbs_i.shape[0]):
                    for nu in range(new_atom_orbs_j.shape[0]):
                        D_mn = density_matrix[ao_indices_i[mu],
                                              ao_indices_j[nu]]
                        np_density[
                            t[0]:t[0] + no[0],
                            t[1]:t[1] + no[1],
                            t[2]:t[2] + no[2],
                        ] += D_mn * new_atom_orbs_i[
                            mu,
                            p1_i[0]:p1_i[0] + no[0],
                            p1_i[1]:p1_i[1] + no[1],
                            p1_i[2]:p1_i[2] + no[2],
                        ] * new_atom_orbs_j[
                            nu,
                            p1_j[0]:p1_j[0] + no[0],
                            p1_j[1]:p1_j[1] + no[1],
                            p1_j[2]:p1_j[2] + no[2],
                        ]

        return np_density

    def get_cube_str(self, density_data, comment=''):
        """
        Writes density data to a string in cube file format.

        :param cubefile:
            Name of the cube file.
        :param grid:
            The cubic grid.
        """

        cube_str_list = []

        coords = self._molecule.get_coordinates_in_bohr()

        x = coords[:, 0]
        y = coords[:, 1]
        z = coords[:, 2]

        natoms = self._molecule.number_of_atoms()
        elem_ids = self._molecule.get_identifiers()

        x0, y0, z0 = self.origin
        dx, dy, dz = self.stepsize
        nx, ny, nz = self.npoints

        cube_str_list.append('Cube file generated by VeloxChem\n')

        # cube for density
        cube_str_list.append('Electron {:s} {:s}\n'.format('density', comment))
        line = '{:5d}{:12.6f}{:12.6f}{:12.6f}{:5d}\n'.format(
            natoms, x0, y0, z0, 1)
        cube_str_list.append(line)

        cube_str_list.append('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(
            nx, dx, 0, 0))
        cube_str_list.append('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(
            ny, 0, dy, 0))
        cube_str_list.append('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(
            nz, 0, 0, dz))

        for a in range(natoms):
            line = '{:5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(
                elem_ids[a], float(elem_ids[a]), x[a], y[a], z[a])
            cube_str_list.append(line)

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    cube_str_list.append(' {:12.5E}'.format(density_data[ix, iy,
                                                                         iz]))
                    if iz % 6 == 5:
                        cube_str_list.append('\n')
                cube_str_list.append('\n')

        return ''.join(cube_str_list)

    def plot(self, molecule, basis, den_inp):
        """
        Plots the densities, with a widget to choose which.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param den_inp:
            The density matrix input (filename of h5 file storing the
            densities, or a dictionary of density matrices (numpy arrays)
            and their labels).
        """

        try:
            import ipywidgets as widgets
        except ImportError:
            raise ImportError(self.help_string_widgets_and_display())

        if isinstance(den_inp, str):
            den_dict = self.read_hdf5(den_inp)
        elif isinstance(den_inp, dict):
            den_dict = den_inp
        else:
            assert_msg_critical(False,
                                'DensityViewer.plot: Invalid density input')

        self._density_dict = den_dict

        self.initialize(molecule, basis)

        den_key_list = []
        for key in self._density_dict:
            den_key_list.append(key)

        widgets.interact(self.draw_density, density_label=den_key_list)

    def draw_density(self, density_label):

        try:
            import py3Dmol
        except ImportError:
            raise ImportError('Unable to import py3Dmol')

        density_np = self._density_dict[density_label]
        density_cube_data = self.compute_density(density_np)
        density_cube_str = self.get_cube_str(density_cube_data)

        viewer = py3Dmol.view(width=600, height=450)

        viewer.addModel(density_cube_str, "cube")
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.1}})

        isovalue = self.density_isovalue
        opacity = self.density_opacity
        if self.density_color_scheme == 'default':
            color = 0x0000ff
        elif self.density_color_scheme == 'alternative':
            color = 0x62a0ea

        viewer.addVolumetricData(density_cube_str, "cube", {
            "isoval": isovalue,
            "color": color,
            "opacity": opacity,
        })

        isovalue = -isovalue
        if self.density_color_scheme == 'default':
            color = 0xff0000
        elif self.density_color_scheme == 'alternative':
            color = 0xe5a50a

        viewer.addVolumetricData(density_cube_str, "cube", {
            "isoval": isovalue,
            "color": color,
            "opacity": opacity,
        })

        viewer.zoomTo()
        viewer.show()

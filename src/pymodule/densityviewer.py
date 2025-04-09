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


class DensityViewer:
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
        self._i_den = None
        self._density_list = None
        self._density_labels = None
        self._is_uhf = False

        # plot
        self._this_plot = None
        self._plt_iso_one = None
        self._plt_iso_two = None

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

    def help_string_k3d(self):

        # TODO: double check if the nbextension part is still needed

        return ('Unable to import k3d. Please install k3d via pip or conda,\n' +
                '  and then run\n' +
                '  jupyter nbextension install --py --sys-prefix k3d\n' +
                '  jupyter nbextension enable --py --sys-prefix k3d\n')

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
            import k3d
        except ImportError:
            raise ImportError(self.help_string_k3d())

        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            raise ImportError(self.help_string_widgets_and_display())

        if isinstance(den_inp, str):
            den_dict = self.read_hdf5(den_inp)
        elif isinstance(den_inp, dict):
            den_dict = den_inp
        else:
            assert_msg_critical(False,
                                'DensityViewer.plot: Invalid density input')

        self.initialize(molecule, basis)

        self._density_labels = []
        self._density_list = []
        for label in den_dict.keys():
            self._density_labels.append(label)
            self._density_list.append(den_dict[label])

        self._i_den = 0

        self._this_plot = k3d.plot(grid_visible=False)
        plt_atoms, plt_bonds = self.draw_molecule(molecule)
        self._this_plot += plt_atoms
        for bonds in plt_bonds:
            self._this_plot += bonds

        density = self.compute_density(self._density_list[self._i_den])
        self._plt_iso_one, self._plt_iso_two = self.draw_density(density)
        self._this_plot += self._plt_iso_one
        self._this_plot += self._plt_iso_two
        self._this_plot.display()

        den_list = []
        for i, label in enumerate(self._density_labels):
            den_list.append((label, i))

        # Add widget
        self.density_selector = widgets.Dropdown(options=den_list,
                                                 value=self._i_den,
                                                 description='Density:')
        display(self.density_selector)
        self.density_selector.observe(self.on_density_index_change,
                                      names='value')

    def on_density_index_change(self, change):
        """
        Registers a widget event to plot a different density.

        :param change:
            A dictionary created by the widget observe.
        """

        i_den = change['new']

        # Do not do anything if 'blank index' is chosen
        if i_den < 0:
            return

        self._i_den = i_den

        # To avoid disturbing the current view
        self._this_plot.camera_auto_fit = False
        self._this_plot -= self._plt_iso_one
        self._this_plot -= self._plt_iso_two

        density = self.compute_density(self._density_list[self._i_den])

        self._plt_iso_one, self._plt_iso_two = self.draw_density(density)
        self._this_plot += self._plt_iso_one
        self._this_plot += self._plt_iso_two
        self._this_plot.render()

    def draw_molecule(self, molecule):
        """
        Draws a ball-and-stick model of the molecule.

        :param molecule:
            The molecule.

        :return:
            A tuple with:
                - a k3d points objects for all atoms
                - a list of k3d line objects, 2 per bonds
        """

        try:
            import k3d
        except ImportError:
            raise ImportError(self.help_string_k3d())

        atomcolor = np.array(
            [
                0xFFFFFF, 0xD9FFFF, 0xCC80FF, 0xC2FF00, 0xFFB5B5, 0x909090,
                0x3050F8, 0xFF0D0D, 0x90E050, 0xB3E3F5, 0xAB5CF2, 0x8AFF00,
                0xBFA6A6, 0xF0C8A0, 0xFF8000, 0xFFFF30, 0x1FF01F, 0x80D1E3,
                0x8F40D4, 0x3DFF00, 0xE6E6E6, 0xBFC2C7, 0xA6A6AB, 0x8A99C7,
                0x9C7AC7, 0xE06633, 0xF090A0, 0x50D050, 0xC88033, 0x7D80B0,
                0xC28F8F, 0x668F8F, 0xBD80E3, 0xFFA100, 0xA62929, 0x5CB8D1,
                0x702EB0, 0x00FF00, 0x94FFFF, 0x94E0E0, 0x73C2C9, 0x54B5B5,
                0x3B9E9E, 0x248F8F, 0x0A7D8C, 0x006985, 0xC0C0C0, 0xFFD98F,
                0xA67573, 0x668080, 0x9E63B5, 0xD47A00, 0x940094, 0x429EB0,
                0x57178F, 0x00C900, 0x70D4FF, 0xFFFFC7, 0xD9FFC7, 0xC7FFC7,
                0xA3FFC7, 0x8FFFC7, 0x61FFC7, 0x45FFC7, 0x30FFC7, 0x1FFFC7,
                0x00FF9C, 0x00E675, 0x00D452, 0x00BF38, 0x00AB24, 0x4DC2FF,
                0x4DA6FF, 0x2194D6, 0x267DAB, 0x266696, 0x175487, 0xD0D0E0,
                0xFFD123, 0xB8B8D0, 0xA6544D, 0x575961, 0x9E4FB5, 0xAB5C00,
                0x754F45, 0x428296, 0x420066, 0x007D00, 0x70ABFA, 0x00BAFF,
                0x00A1FF, 0x008FFF, 0x0080FF, 0x006BFF, 0x545CF2, 0x785CE3,
                0x8A4FE3, 0xA136D4, 0xB31FD4, 0xB31FBA, 0xB30DA6, 0xBD0D87,
                0xC70066, 0xCC0059, 0xD1004F, 0xD90045, 0xE00038, 0xE6002E,
                0xEB0026
            ],
            dtype='uint32',
        )

        connectivity_matrix = molecule.get_connectivity_matrix()

        natoms = molecule.number_of_atoms()
        coords = molecule.get_coordinates_in_bohr().astype('float32')

        # Create a list of colors
        colors = []
        for nr in self._atomnr:
            if nr < len(atomcolor):
                colors.append(atomcolor[nr])
            else:
                colors.append(atomcolor[-1])

        # Balls
        plt_atoms = k3d.points(positions=coords,
                               point_size=0.7,
                               colors=colors,
                               shader='mesh',
                               name='Atoms',
                               group='Molecule')

        # Sticks
        # Create a lines object for each atom type
        bonddict = {}
        labels = molecule.get_labels()
        names = {}
        for i in range(natoms):
            atom_id = self._atomnr[i]
            if atom_id not in bonddict:
                bonddict[atom_id] = []
            if atom_id not in names:
                names[atom_id] = labels[i]

        newcoords = []
        ncoords = natoms
        for i in range(natoms):
            for j in range(i + 1, natoms):
                # Check if there is a bond
                if connectivity_matrix[i, j] != 1:
                    continue
                # If single atom type, just record it
                if self._atomnr[i] == self._atomnr[j]:
                    bonddict[self._atomnr[i]].append([i, j])
                # Else do 2 segments (which means adding a new middle-vertex)
                else:
                    newcoords.append(0.5 * (coords[i, :] + coords[j, :]))
                    bonddict[self._atomnr[i]].append([i, ncoords])
                    bonddict[self._atomnr[j]].append([ncoords, j])
                    ncoords += 1

        finalcoords = np.append(coords, newcoords)
        plt_bonds = []
        for atom, bondlist in bonddict.items():
            plt_bonds.append(
                k3d.lines(finalcoords,
                          bondlist,
                          width=0.35,
                          color=int(atomcolor[atom]),
                          shader='mesh',
                          radial_segments=16,
                          indices_type='segment',
                          name=names[atom] + ' bonds',
                          group='Molecule'))

        return plt_atoms, plt_bonds

    def draw_density(self, density):
        """
        Draws an isosurface of the density.

        :param density:
            The density on the grid.

        :return:
            A tuple with the k3d positive and negative isosurfaces.
        """

        try:
            import k3d
        except ImportError:
            raise ImportError(self.help_string_k3d())

        # There is a mismatch between the x,y,z order of vlx and k3d
        density_k3d = density.swapaxes(0, 2).astype('float32')

        # Plot the density
        xmin, ymin, zmin = self._grid.get_origin()
        dx, dy, dz = self._grid.get_step_size()
        nx, ny, nz = self._grid.get_num_points()
        xmax = xmin + (nx - 1) * dx
        ymax = ymin + (ny - 1) * dy
        zmax = zmin + (nz - 1) * dz

        bounds = [xmin, xmax, ymin, ymax, zmin, zmax]

        # Default settings
        isovalue = self.density_isovalue
        opacity = self.density_opacity
        wireframe = False

        if self.density_color_scheme == 'default':
            color = 0x0000ff
        elif self.density_color_scheme == 'alternative':
            color = 0x62a0ea

        # Find if the user changed the defaults
        if self._plt_iso_one:
            isovalue = self._plt_iso_one.level
            opacity = self._plt_iso_one.opacity
            wireframe = self._plt_iso_one.wireframe
            color = self._plt_iso_one.color

        plt_iso_one = k3d.marching_cubes(density_k3d,
                                         compression_level=9,
                                         bounds=bounds,
                                         level=isovalue,
                                         flat_shading=False,
                                         opacity=opacity,
                                         wireframe=wireframe,
                                         color=color,
                                         name='Positive isosurface',
                                         group='Densities')

        # Find if the user changed the defaults
        isovalue = -isovalue

        if self.density_color_scheme == 'default':
            color = 0xff0000
        elif self.density_color_scheme == 'alternative':
            color = 0xe5a50a

        if self._plt_iso_two:
            isovalue = self._plt_iso_two.level
            opacity = self._plt_iso_two.opacity
            wireframe = self._plt_iso_two.wireframe
            color = self._plt_iso_two.color

        plt_iso_two = k3d.marching_cubes(density_k3d,
                                         compression_level=9,
                                         bounds=bounds,
                                         level=isovalue,
                                         flat_shading=False,
                                         opacity=opacity,
                                         wireframe=wireframe,
                                         color=color,
                                         name='Negative isosurface',
                                         group='Densities')

        return plt_iso_one, plt_iso_two

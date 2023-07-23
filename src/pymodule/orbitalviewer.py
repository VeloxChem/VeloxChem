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
import math

from .veloxchemlib import VisualizationDriver, CubicGrid
from .veloxchemlib import molorb
from .veloxchemlib import bohr_in_angstrom


class OrbitalViewer:
    """
    An approximate but very quick module to compute a molecular orbital on a
    grid. It does this by computing the AOs of each unique atom type in an
    atomic grid. These AOs are then translated as needed to the nearest grid
    point of the actual atomic center to compute on-the-fly the molecular
    orbital on the full grid.

    Instance variables
        - grid_density: The density of grid points per a.u.
        - grid_margins: Minimal distance in a.u. between the box edge and all
          atomic nuclei
        - atombox_radius: Radius in a.u. of the atomix grid
        - mo_threshold: Do not compute AO contribution if corresponding MO
          coefficient is below this value
    """

    def __init__(self):
        """
        Initializes the orbital viewer.
        """

        # number of grid points per a.u.
        self.grid_density = 4

        # a.u. added around the molecule to create the box
        self.grid_margins = 3.0

        # Radius in a.u. for each atomic box
        self.atombox_radius = 6.0

        # Do not compute AO if MO coefficient is below this value
        self.mo_threshold = 0.01

        # molecule
        self.molecule = None
        self.basis = None
        self.atomnr = None
        self.coords = None

        # grid
        self.grid = None
        self.origin = None
        self.stepsize = None
        self.npoints = None
        self.atom_origin = None
        self.atom_npoints = None

        # orbital
        self.ao_to_atom = None
        self.ao_dict = None
        self.i_orb = None
        self.mo_coefs = None
        self.is_uhf = False

        # plot
        self.this_plot = None
        self.plt_iso_one = None
        self.plt_iso_two = None

    def help_string_k3d(self):

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

        self.molecule = molecule
        self.basis = basis

        # Define box size
        self.atomnr = molecule.elem_ids_to_numpy()
        self.coords = molecule.get_coordinates_in_bohr()
        xmin = self.coords[:, 0].min() - self.grid_margins
        xmax = self.coords[:, 0].max() + self.grid_margins
        ymin = self.coords[:, 1].min() - self.grid_margins
        ymax = self.coords[:, 1].max() + self.grid_margins
        zmin = self.coords[:, 2].min() - self.grid_margins
        zmax = self.coords[:, 2].max() + self.grid_margins
        nx = math.ceil((xmax - xmin) * self.grid_density)
        ny = math.ceil((ymax - ymin) * self.grid_density)
        nz = math.ceil((zmax - zmin) * self.grid_density)
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        dz = (zmax - zmin) / nz

        self.origin = (xmin, ymin, zmin)
        self.stepsize = (dx, dy, dz)
        self.npoints = (nx, ny, nz)
        self.grid = CubicGrid(self.origin, self.stepsize, self.npoints)

        # Create a small box for each atom
        nx_atom = math.ceil(self.atombox_radius / dx * 2)
        ny_atom = math.ceil(self.atombox_radius / dy * 2)
        nz_atom = math.ceil(self.atombox_radius / dz * 2)
        self.atom_origin = (-nx_atom * dx / 2, -ny_atom * dy / 2,
                            -nz_atom * dz / 2)
        self.atom_npoints = (nx_atom, ny_atom, nz_atom)

        atom_grid = CubicGrid(self.atom_origin, self.stepsize,
                              self.atom_npoints)

        # Create atomic map
        vis_drv = VisualizationDriver()
        ao_info = vis_drv.get_atomic_orbital_info(molecule, basis)

        # Create reverse atomic map
        atom_to_ao = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)
        self.ao_to_atom = [[] for i in range(len(ao_info))]
        for i_atom, atom_orbs in enumerate(atom_to_ao):
            for i_orb, orb in enumerate(atom_orbs):
                self.ao_to_atom[orb] = [i_atom, i_orb]

        # Compute each unique AO on the grid
        self.ao_dict = {}
        for i_atom, atom in enumerate(self.atomnr):
            if atom not in self.ao_dict:
                atomlist = []
                for orb in atom_to_ao[i_atom]:
                    vis_drv.compute_atomic_orbital_for_grid(
                        atom_grid, basis, ao_info[orb])
                    atomlist.append(atom_grid.values_to_numpy())
                self.ao_dict[atom] = atomlist

    def compute_orbital(self, orbital, index):
        """
        Compute a specific molecular orbital.

        :param orbital:
            A numpy array containing MO coefficients.
        :param index:
            The index of the MO of interest.

        :return:
            A numpy array with the value of the MO on the grid.
        """

        this_orb = orbital[:, index]

        # Create full grid
        nx, ny, nz = self.atom_npoints
        np_orb = np.zeros(self.npoints)

        # Loop over AOs
        for i_coef, coef in enumerate(this_orb):
            if abs(coef) > self.mo_threshold:
                i_atom, i_orb = self.ao_to_atom[i_coef]
                this_atom = self.atomnr[i_atom]
                atom_orb = self.ao_dict[this_atom][i_orb]
                t = np.round(
                    (self.coords[i_atom] - self.origin + self.atom_origin) /
                    self.stepsize).astype('int')
                ncopy = [nx, ny, nz]
                p1 = [0, 0, 0]
                for i in range(3):
                    if t[i] < 0:
                        p1[i] = -t[i]
                        t[i] = 0
                        ncopy[i] -= p1[i]
                    if t[i] + ncopy[i] >= self.npoints[i]:
                        ncopy[i] = self.npoints[i] - t[i]
                np_orb[t[0]:t[0] + ncopy[0], t[1]:t[1] + ncopy[1], t[2]:t[2] +
                       ncopy[2]] += coef * atom_orb[p1[0]:p1[0] + ncopy[0],
                                                    p1[1]:p1[1] + ncopy[1],
                                                    p1[2]:p1[2] + ncopy[2]]

        return np_orb

    def plot(self, molecule, basis, mo_object):
        """
        Plots the orbitals, with a widget to choose which.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mo_object:
            The MolecularOrbitals object.
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

        self.initialize(molecule, basis)

        self.is_uhf = (mo_object.get_orbitals_type() == molorb.unrest)

        # i_orb is an instance variable accessed by MultiPsi
        self.i_orb = self.molecule.number_of_alpha_electrons() - 1
        self.mo_coefs = mo_object.alpha_to_numpy()

        # In some cases (for example NTOs) the number of orbitals is less than
        # the number of electrons. In this case, print the middle orbital
        if self.mo_coefs.shape[1] < self.molecule.number_of_alpha_electrons():
            self.i_orb = self.mo_coefs.shape[1] // 2
        if self.is_uhf:
            self.mo_coefs_beta = mo_object.beta_to_numpy()
        else:
            self.mo_coefs_beta = self.mo_coefs

        self.this_plot = k3d.plot(grid_visible=False)
        plt_atoms, plt_bonds = self.draw_molecule(self.molecule)
        self.this_plot += plt_atoms
        for bonds in plt_bonds:
            self.this_plot += bonds

        orbital = self.compute_orbital(self.mo_coefs, self.i_orb)
        self.plt_iso_one, self.plt_iso_two = self.draw_orbital(orbital)
        self.this_plot += self.plt_iso_one
        self.this_plot += self.plt_iso_two
        self.this_plot.display()

        # Create orbital list:
        orb_ene = mo_object.ea_to_numpy()
        orb_occ = mo_object.occa_to_numpy()
        orb_occ_beta = mo_object.occb_to_numpy()
        if not self.is_uhf:
            # Check if MO is NTO. In case of NTO, only print alpha occupation
            # numbers which are actually lambda's
            is_nto = (np.min(orb_ene) == 0.0 and np.max(orb_ene) == 0.0 and
                      np.min(orb_occ) < 0.0)
            if not is_nto:
                orb_occ += orb_occ_beta
        orblist = []
        for i in range(len(orb_ene)):
            orb_label = f'{i+1:3d} occ={orb_occ[i]:.3f} '
            orb_label += f'ene={orb_ene[i]:.3f}'
            orblist.append((orb_label, i))

        # Also do for beta if UHF
        if self.is_uhf:
            orb_ene_beta = mo_object.eb_to_numpy()
            orblist_beta = [('', -1)]
            for i in range(len(orb_ene_beta)):
                orb_label = f'{i+1:3d} occ={orb_occ_beta[i]:.3f} '
                orb_label += f'ene={orb_ene_beta[i]:.3f}'
                orblist_beta.append((orb_label, i))

            # Add empty space
            orblist.insert(0, ('', -1))
            # Add widget
            self.orbital_selector = widgets.Dropdown(
                options=orblist, value=self.i_orb, description='Alpha orbital:')
            self.orbital_selector_beta = widgets.Dropdown(
                options=orblist_beta, value=-1, description='Beta orbital:')
            hbox = widgets.HBox(
                [self.orbital_selector, self.orbital_selector_beta])
            display(hbox)
            self.orbital_selector_beta.observe(
                self.on_orbital_index_change_beta, names='value')
        else:
            # Add widget
            self.orbital_selector = widgets.Dropdown(options=orblist,
                                                     value=self.i_orb,
                                                     description='Orbital:')
            display(self.orbital_selector)
        self.orbital_selector.observe(self.on_orbital_index_change,
                                      names='value')

    def on_orbital_index_change(self, change, spin='alpha'):
        """
        Registers a widget event to plot a different orbital.

        :param change:
            A dictionary created by the widget observe.
        :param spin:
            The spin of the orbital.
        """

        i_orb = change['new']

        # Do not do anything if 'blank index' is chosen
        if i_orb < 0:
            return

        self.i_orb = i_orb

        # To avoid disturbing the current view
        self.this_plot.camera_auto_fit = False
        self.this_plot -= self.plt_iso_one
        self.this_plot -= self.plt_iso_two

        if spin == 'alpha':
            if self.is_uhf:
                # Reset the beta index to blank if exists
                self.orbital_selector_beta.value = -1
            orbital = self.compute_orbital(self.mo_coefs, self.i_orb)
        else:
            # Reset the alpha index to blank
            self.orbital_selector.value = -1
            orbital = self.compute_orbital(self.mo_coefs_beta, self.i_orb)

        self.plt_iso_one, self.plt_iso_two = self.draw_orbital(orbital)
        self.this_plot += self.plt_iso_one
        self.this_plot += self.plt_iso_two
        self.this_plot.render()

    def on_orbital_index_change_beta(self, change):
        """
        Registers a widget event to plot a different orbital of beta-spin.

        :param change:
            A dictionary created by the widget observe.
        """

        self.on_orbital_index_change(change, 'beta')

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

        atomradius = [
            0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56, 0.48, 0.42, 0.38, 1.90,
            1.45, 1.18, 1.11, 0.98, 0.88, 0.79, 0.71, 2.43, 1.94, 1.84, 1.76,
            1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25, 1.14,
            1.03, 0.94, 0.88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.90, 1.83, 1.78,
            1.73, 1.69, 1.65, 1.61, 1.56, 1.45, 1.33, 1.23, 1.15, 1.08, 2.98,
            2.53, 1.95, 1.85, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33, 2.25, 2.28,
            2.26, 2.26, 2.22, 2.22, 2.17, 2.08, 2.00, 1.93, 1.88, 1.85, 1.80,
            1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, 1.27, 1.20
        ]

        natoms = molecule.number_of_atoms()
        atomnr = molecule.elem_ids_to_numpy() - 1
        coords = molecule.get_coordinates_in_bohr().astype('float32')

        # Create a list of colors and radii
        colors = []
        radii = []
        for nr in atomnr:
            if nr < len(atomcolor):
                colors.append(atomcolor[nr])
            else:
                colors.append(atomcolor[-1])
            if nr < len(atomradius):
                radii.append(atomradius[nr])
            else:
                radii.append(2.0)

        # Balls
        plt_atoms = k3d.points(positions=coords,
                               point_size=0.7,
                               colors=colors,
                               shader='mesh',
                               name='Atoms',
                               group='Molecule')

        # Sticks
        plt_bonds = []
        for i in range(natoms):
            color_i = colors[i]
            for j in range(i + 1, natoms):
                bond = (radii[i] + radii[j]) / bohr_in_angstrom()
                if np.linalg.norm(coords[i, :] - coords[j, :]) > 1.25 * bond:
                    continue
                color_j = colors[j]
                plt_bonds.append(
                    k3d.line(
                        [coords[i, :], 0.5 * (coords[i, :] + coords[j, :])],
                        width=0.35,
                        colors=[color_i, color_i],
                        shader='mesh',
                        radial_segments=16,
                        group='Molecule'))
                plt_bonds.append(
                    k3d.line(
                        [0.5 * (coords[i, :] + coords[j, :]), coords[j, :]],
                        width=0.35,
                        colors=[color_j, color_j],
                        shader='mesh',
                        radial_segments=16,
                        group='Molecule'))

        return plt_atoms, plt_bonds

    def draw_orbital(self, orbital):
        """
        Draws an isosurface of the orbitals.

        :param orbital:
            The molecular orbital on the grid.
        :param isovalue:
            The isovalue for isosurfaces.

        :return:
            A tuple with the k3d positive and negative isosurfaces.
        """

        try:
            import k3d
        except ImportError:
            raise ImportError(self.help_string_k3d())

        # There is a mismatch between the x,y,z order of vlx and k3d
        orbital_k3d = orbital.swapaxes(0, 2).astype('float32')

        # Plot orbitals
        xmin, ymin, zmin = self.grid.get_origin()
        dx, dy, dz = self.grid.get_step_size()
        nx, ny, nz = self.grid.get_num_points()
        xmax = xmin + (nx - 1) * dx
        ymax = ymin + (ny - 1) * dy
        zmax = zmin + (nz - 1) * dz

        bounds = [xmin, xmax, ymin, ymax, zmin, zmax]

        # Default settings
        isovalue = 0.05
        opacity = 1.0
        wireframe = False
        color = 0x0000ff

        # Find if the user changed the defaults
        if self.plt_iso_one:
            isovalue = self.plt_iso_one.level
            opacity = self.plt_iso_one.opacity
            wireframe = self.plt_iso_one.wireframe
            color = self.plt_iso_one.color

        plt_iso_one = k3d.marching_cubes(orbital_k3d,
                                         compression_level=9,
                                         bounds=bounds,
                                         level=isovalue,
                                         flat_shading=False,
                                         opacity=opacity,
                                         wireframe=wireframe,
                                         color=color,
                                         name='Positive isosurface',
                                         group='Orbitals')

        # Find if the user changed the defaults
        isovalue = -isovalue
        color = 0xff0000
        if self.plt_iso_two:
            isovalue = self.plt_iso_two.level
            opacity = self.plt_iso_two.opacity
            wireframe = self.plt_iso_two.wireframe
            color = self.plt_iso_two.color

        plt_iso_two = k3d.marching_cubes(orbital_k3d,
                                         compression_level=9,
                                         bounds=bounds,
                                         level=isovalue,
                                         flat_shading=False,
                                         opacity=opacity,
                                         wireframe=wireframe,
                                         color=color,
                                         name='Negative isosurface',
                                         group='Orbitals')

        return plt_iso_one, plt_iso_two

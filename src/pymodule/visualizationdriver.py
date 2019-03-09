from .veloxchemlib import VisualizationDriver
from .veloxchemlib import CubicGrid

from .molecularorbitals import MolecularOrbitals
from .aodensitymatrix import AODensityMatrix

import multiprocessing as mp
import numpy as np
import os


@staticmethod
def _VisualizationDriver_gen_grid(molecule, n_x=80, n_y=80, n_z=80):

    x = molecule.x_to_numpy()
    y = molecule.y_to_numpy()
    z = molecule.z_to_numpy()

    atom_radii = molecule.vdw_radii_to_numpy()

    x_max = np.max(x + atom_radii * 2.0)
    x_min = np.min(x - atom_radii * 2.0)

    y_max = np.max(y + atom_radii * 2.0)
    y_min = np.min(y - atom_radii * 2.0)

    z_max = np.max(z + atom_radii * 2.0)
    z_min = np.min(z - atom_radii * 2.0)

    x_step = (x_max - x_min) / float(n_x)
    y_step = (y_max - y_min) / float(n_y)
    z_step = (z_max - z_min) / float(n_z)

    return CubicGrid([x_min, y_min, z_min], [x_step, y_step, z_step],
                     [n_x, n_y, n_z])


@staticmethod
def _VisualizationDriver_write_data(cubefile, molecule, flag, index, spin, grid,
                                    data):

    f_cube = open(cubefile, 'w')

    x = molecule.x_to_numpy()
    y = molecule.y_to_numpy()
    z = molecule.z_to_numpy()

    natoms = molecule.number_of_atoms()
    elem_ids = molecule.elem_ids_to_numpy()

    x0, y0, z0 = grid.x_origin(), grid.y_origin(), grid.z_origin()
    dx, dy, dz = grid.x_step_size(), grid.y_step_size(), grid.z_step_size()
    nx, ny, nz = grid.x_num_points(), grid.y_num_points(), grid.z_num_points()

    f_cube.write('VeloxChem Cube File\n')

    if flag == 'mo':
        f_cube.write('MO {:d}, {:s} spin\n'.format(index + 1, spin))
        f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}{:5d}\n'.format(
            -natoms, x0, y0, z0, 1))

    elif flag == 'density':
        f_cube.write('Electron density, {:s} spin\n'.format(spin))
        f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}{:5d}\n'.format(
            natoms, x0, y0, z0, 1))

    f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(nx, dx, 0, 0))
    f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(ny, 0, dy, 0))
    f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(nz, 0, 0, dz))

    for a in range(natoms):
        f_cube.write('{:5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(
            elem_ids[a], float(elem_ids[a]), x[a], y[a], z[a]))

    if flag == 'mo':
        f_cube.write('{:5d}{:5d}\n'.format(1, index + 1))

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                f_cube.write(' {:12.5E}'.format(data[ix, iy, iz]))
                if iz % 6 == 5:
                    f_cube.write("\n")
            f_cube.write("\n")

    f_cube.close()


def _VisualizationDriver_write_cube(self, cubefile, molecule, basis,
                                    mo_or_density, index, spin, grid):

    if isinstance(mo_or_density, MolecularOrbitals):
        flag = 'mo'

    elif isinstance(mo_or_density, AODensityMatrix):
        flag = 'density'

    else:
        errmsg = 'VisualizationDriver.write_cube: invalide argument'
        assert_msg_critical(False, errmsg)

    data = self.compute(molecule, basis, mo_or_density, index, spin, grid)

    self.write_data(cubefile, molecule, flag, index, spin, grid, data)


VisualizationDriver.gen_grid = _VisualizationDriver_gen_grid
VisualizationDriver.write_data = _VisualizationDriver_write_data
VisualizationDriver.write_cube = _VisualizationDriver_write_cube

import numpy as np
import h5py
import re

from .veloxchemlib import VisualizationDriver
from .veloxchemlib import CubicGrid
from .veloxchemlib import DenseMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import ao_matrix_to_veloxchem
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .errorhandler import assert_msg_critical


@staticmethod
def _VisualizationDriver_gen_cubic_grid(molecule, n_x=80, n_y=80, n_z=80):
    """
    Creates cubic grid for a molecule.

    :param molecule:
        The molecule.
    :param n_x:
        Number of grid points in X direction.
    :param n_y:
        Number of grid points in Y direction.
    :param n_z:
        Number of grid points in Z direction.

    :return:
        The cubic grid.
    """

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
def _VisualizationDriver_write_data(cubefile, grid, molecule, flag, index,
                                    spin):
    """
    Writes cubic grid data to cube file.

    :param cubefile:
        Name of the cube file.
    :param grid:
        The cubic grid.
    :param molecule:
        The molecule.
    :param flag:
        The flag for the type of cube data (mo, nto, density).
    :param index:
        Index (0-based) of the molecular orbital or density matrix.
    :param spin:
        Spin of the molecular orbital or density.
    """

    f_cube = open(cubefile, 'w')

    x = molecule.x_to_numpy()
    y = molecule.y_to_numpy()
    z = molecule.z_to_numpy()

    natoms = molecule.number_of_atoms()
    elem_ids = molecule.elem_ids_to_numpy()

    x0, y0, z0 = grid.x_origin(), grid.y_origin(), grid.z_origin()
    dx, dy, dz = grid.x_step_size(), grid.y_step_size(), grid.z_step_size()
    nx, ny, nz = grid.x_num_points(), grid.y_num_points(), grid.z_num_points()

    print('VeloxChem Cube File', file=f_cube)

    if flag in ['mo', 'nto']:
        print('{:s} {:d} ({:s})'.format(flag.upper(), index + 1, spin),
              file=f_cube)
        print('{:5d}{:12.6f}{:12.6f}{:12.6f}{:5d}'.format(
            -natoms, x0, y0, z0, 1),
              file=f_cube)

    elif flag == 'density':
        print('Electron density ({:s})'.format(spin), file=f_cube)
        print('{:5d}{:12.6f}{:12.6f}{:12.6f}{:5d}'.format(
            natoms, x0, y0, z0, 1),
              file=f_cube)

    print('{:5d}{:12.6f}{:12.6f}{:12.6f}'.format(nx, dx, 0, 0), file=f_cube)
    print('{:5d}{:12.6f}{:12.6f}{:12.6f}'.format(ny, 0, dy, 0), file=f_cube)
    print('{:5d}{:12.6f}{:12.6f}{:12.6f}'.format(nz, 0, 0, dz), file=f_cube)

    for a in range(natoms):
        print('{:5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}'.format(
            elem_ids[a], float(elem_ids[a]), x[a], y[a], z[a]),
              file=f_cube)

    if flag in ['mo', 'nto']:
        print('{:5d}{:5d}'.format(1, index + 1), file=f_cube)

    data = grid.values_to_numpy()

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                print(' {:12.5E}'.format(data[ix, iy, iz]), end='', file=f_cube)
                if iz % 6 == 5:
                    print('', file=f_cube)
            print('', file=f_cube)

    f_cube.close()


def _VisualizationDriver_gen_cubes(self, cube_dict, molecule, basis, mol_orbs,
                                   density):
    """
    Computes and writes cube file.

    :param cube_dict:
        The input dictionary of cube options.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param mol_orbs:
        The molecular orbitals.
    :param density:
        The density matrix.
    """

    if 'grid' in cube_dict:
        grid = [int(x) for x in cube_dict['grid'].split(',')]
    else:
        grid = [80, 80, 80]
    cubic_grid = self.gen_cubic_grid(molecule, *grid[:3])

    cubes = [x.strip() for x in cube_dict['cubes'].split(',')]
    if 'files' in cube_dict:
        files = [x.strip() for x in cube_dict['files'].split(',')]
    else:
        files = ['cube_{:d}.cube'.format(i + 1) for i in range(len(cubes))]

    assert_msg_critical(
        len(cubes) == len(files),
        'VisualizationDriver: inconsistent number of cubes')

    for cube, fname in zip(cubes, files):

        m = re.search(r'^(.*)\((.*)\)$', cube)

        assert_msg_critical(m is not None,
                            'VisualizationDriver: failed to read cube inputs')

        cube_type = m.group(1).strip().lower()

        if cube_type in ['mo', 'amo', 'bmo']:

            if cube_type in ['mo', 'amo']:
                spin = 'alpha'
                nelec = molecule.number_of_alpha_electrons()
            else:
                spin = 'beta'
                nelec = molecule.number_of_beta_electrons()

            homo = nelec - 1
            lumo = nelec

            cube_value = m.group(2).strip().lower()
            cube_value = cube_value.replace('homo', str(homo))
            cube_value = cube_value.replace('lumo', str(lumo))
            orb_id = eval(cube_value)

            self.compute(cubic_grid, molecule, basis, mol_orbs, orb_id, spin)

            if self.get_rank() == mpi_master():
                self.write_data(fname, cubic_grid, molecule, 'mo', orb_id, spin)

        elif cube_type == 'density':

            cube_value = m.group(2).strip().lower()
            spin = cube_value

            self.compute(cubic_grid, molecule, basis, density, 0, spin)

            if self.get_rank() == mpi_master():
                self.write_data(fname, cubic_grid, molecule, 'density', 0, spin)

        elif cube_type == 'read_dalton':

            cube_value = m.group(2).strip()

            hf = h5py.File(cube_value, 'r')
            dal_dens = DenseMatrix(np.array(hf.get('DALTON_AO_MATRIX')))
            vlx_dens = ao_matrix_to_veloxchem(dal_dens, basis,
                                              molecule).to_numpy()
            read_density = AODensityMatrix([vlx_dens], denmat.rest)
            hf.close()

            spin = 'alpha'
            self.compute(cubic_grid, molecule, basis, read_density, 0, spin)

            if self.get_rank() == mpi_master():
                self.write_data(fname, cubic_grid, molecule, 'density', 0, spin)


VisualizationDriver.gen_cubic_grid = _VisualizationDriver_gen_cubic_grid
VisualizationDriver.write_data = _VisualizationDriver_write_data
VisualizationDriver.gen_cubes = _VisualizationDriver_gen_cubes

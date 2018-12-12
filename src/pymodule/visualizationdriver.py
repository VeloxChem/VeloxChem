from .veloxchemlib import VisualizationDriver
import numpy as np


@staticmethod
def _gen_grid(mol, n_x=80, n_y=80, n_z=80):

    x = mol.x_to_numpy()
    y = mol.y_to_numpy()
    z = mol.z_to_numpy()

    atom_radii = mol.vdw_radii_to_numpy()

    x_max = np.max(x + atom_radii * 2.0)
    x_min = np.min(x - atom_radii * 2.0)

    y_max = np.max(y + atom_radii * 2.0)
    y_min = np.min(y - atom_radii * 2.0)

    z_max = np.max(z + atom_radii * 2.0)
    z_min = np.min(z - atom_radii * 2.0)

    x_step = (x_max - x_min) / float(n_x)
    y_step = (y_max - y_min) / float(n_y)
    z_step = (z_max - z_min) / float(n_z)

    v = [x_min, y_min, z_min, x_step, y_step, z_step, n_x, n_y, n_z]

    return v


def _write_cube(self, mol, basis, mol_orbs, mo_idx, mo_spin, v):

    # TODO: take cube filename from input parser

    outfile = open("outfilename.cube", 'w')

    x = mol.x_to_numpy()
    y = mol.y_to_numpy()
    z = mol.z_to_numpy()

    natoms = mol.number_of_atoms()

    elem_ids = mol.elem_ids_to_numpy()

    x_min = v[0]
    y_min = v[1]
    z_min = v[2]

    x_step = v[3]
    y_step = v[4]
    z_step = v[5]

    n_x = v[6]
    n_y = v[7]
    n_z = v[8]

    outfile.write('VeloxChem Cube File\n')
    outfile.write('MO %d, %s spin\n' % (mo_idx + 1, mo_spin))
    outfile.write('%5d%12.6f%12.6f%12.6f%5d\n' %
                  (-natoms, x_min, y_min, z_min, 1))

    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_x, x_step, 0, 0))
    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_y, 0, y_step, 0))
    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_z, 0, 0, z_step))

    for a in range(natoms):
        outfile.write('%5d%12.6f%12.6f%12.6f%12.6f\n' %
                      (elem_ids[a], float(elem_ids[a]), x[a], y[a], z[a]))

    outfile.write('%5d%5d\n' % (1, mo_idx + 1))

    for ix in range(n_x):
        rx = x_min + x_step * ix

        for iy in range(n_y):
            ry = y_min + y_step * iy

            for iz in range(n_z):
                rz = z_min + z_step * iz

                psi = self.compute(mol, basis, mol_orbs, mo_idx, mo_spin,
                                   rx, ry, rz)

                outfile.write(' %12.5E' % psi)
                if iz % 6 == 5:
                    outfile.write("\n")

            outfile.write("\n")


def _write_cube_dens(self, mol, basis, density, dens_idx, dens_spin, v):

    # TODO: take cube filename from input parser

    outfile = open("outfilename_dens.cube", 'w')

    x = mol.x_to_numpy()
    y = mol.y_to_numpy()
    z = mol.z_to_numpy()

    natoms = mol.number_of_atoms()

    elem_ids = mol.elem_ids_to_numpy()

    x_min = v[0]
    y_min = v[1]
    z_min = v[2]

    x_step = v[3]
    y_step = v[4]
    z_step = v[5]

    n_x = v[6]
    n_y = v[7]
    n_z = v[8]

    outfile.write('VeloxChem Cube File\n')
    outfile.write('Electron density, %s spin\n' % dens_spin)
    outfile.write('%5d%12.6f%12.6f%12.6f%5d\n' %
                  (natoms, x_min, y_min, z_min, 1))

    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_x, x_step, 0, 0))
    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_y, 0, y_step, 0))
    outfile.write('%5d%12.6f%12.6f%12.6f\n' % (n_z, 0, 0, z_step))

    for a in range(natoms):
        outfile.write('%5d%12.6f%12.6f%12.6f%12.6f\n' %
                      (elem_ids[a], float(elem_ids[a]), x[a], y[a], z[a]))

    for ix in range(n_x):
        rx = x_min + x_step * ix

        for iy in range(n_y):
            ry = y_min + y_step * iy

            for iz in range(n_z):
                rz = z_min + z_step * iz

                dens = self.compute(mol, basis, density, dens_idx, dens_spin,
                                    rx, ry, rz)

                outfile.write(' %12.5E' % (dens))
                if iz % 6 == 5:
                    outfile.write("\n")

            outfile.write("\n")


VisualizationDriver.gen_grid = _gen_grid
VisualizationDriver.write_cube = _write_cube
VisualizationDriver.write_cube_dens = _write_cube_dens

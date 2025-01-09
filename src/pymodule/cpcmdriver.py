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

from mpi4py import MPI
import numpy as np
import math
import sys
from scipy.special import erf

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import gen_lebedev_grid
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .oneeints import compute_nuclear_potential_integrals
from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import NuclearPotentialErfDriver
from .veloxchemlib import NuclearPotentialGeom010Driver
from .veloxchemlib import NuclearPotentialGeom100Driver


class CpcmDriver:
    # add for being able to run with alternative grid
    """
    Implements CPCM driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - grid_per_sphere: Number of Lebedev grid points per sphere.
    """

    def __init__(self, gps, qgrid, comm=None, ostream=None):
        """
        Initializes CPCM driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        # grid information
        self.grid_per_sphere = gps
        self.qgrid = qgrid

        # input keywords
        self.input_keywords = {
            'cosmo': {
                'grid_per_sphere': (int, 'number of Lebedev grid points'),
            },
        }

    def compute_solv_energy(self, Bzvec, Cvec, q, molecule=False):
        """
        Computes (electrostatic component of) C-PCM energy.
        TODO: add other components of the energy
        :param molecule:
            The molecule.
        :param Bzvec:
            The nuclear potential on the grid.
        :param Cvec:
            The electronic potential on the grid.

        :return:
            The C-PCM energy.
        """

        return 0.5 * np.vdot(q, Bzvec + Cvec)

    def generate_cosmo_grid(self, molecule):
        """
        Generates Lebedev grid for surface discretization.

        :param molecule:
            The molecule.

        :return:
            TODO
        """

        valid_grid_numbers = [50, 110, 194, 302, 434, 590, 770, 974, 2030]

        assert_msg_critical(
            self.grid_per_sphere in valid_grid_numbers,
            'CosmoDriver.generate_grid: Invalid grid_per_sphere')

        unit_grid = gen_lebedev_grid(self.grid_per_sphere)
        unit_grid_coords = unit_grid[:, :3]
        unit_grid_weights = unit_grid[:, 3:]

        zeta = self.get_zeta_dict()[self.grid_per_sphere]

        # increase radii by 20%
        atom_radii = molecule.vdw_radii_to_numpy() * 1.2
        atom_coords = molecule.get_coordinates_in_bohr()

        cosmo_grid_raw = np.zeros((0, 6))

        for i in range(molecule.number_of_atoms()):
            # scale and shift unit grid
            atom_grid_coords = unit_grid_coords * atom_radii[i] + atom_coords[i]
            grid_zeta = zeta / (atom_radii[i] * np.sqrt(unit_grid_weights))
            atom_idx = np.full_like(grid_zeta, i)
            atom_grid = np.hstack(
                (atom_grid_coords, unit_grid_weights, grid_zeta, atom_idx))
            cosmo_grid_raw = np.vstack((cosmo_grid_raw, atom_grid))

        sw_func_raw = self.get_switching_function(atom_coords, atom_radii, cosmo_grid_raw)
        sw_mask = (sw_func_raw > 1.0e-8)

        cosmo_grid = cosmo_grid_raw[sw_mask, :]
        sw_func = sw_func_raw[sw_mask]

        return cosmo_grid, sw_func

    def get_switching_function(self, atom_coords, atom_radii, grid):

        assert_msg_critical(
            atom_coords.shape[0] == atom_radii.shape[0],
            'CosmoDriver.get_switching_function: Inconsistent atom_coords ' +
            'and atom_radii')

        assert_msg_critical(
            grid.shape[0] == self.grid_per_sphere * atom_coords.shape[0],
            'CosmoDriver.get_switching_function: Should only be used on ' +
            'raw COSMO grid')

        sw_func = np.zeros(grid.shape[0])

        for g in range(grid.shape[0]):
            gx, gy, gz = grid[g, :3]
            zeta_g = grid[g, 4]

            sw_func[g] = 1.0
            atom_idx = g // self.grid_per_sphere

            for i in range(atom_coords.shape[0]):
                if i == atom_idx:
                    continue

                ax, ay, az = atom_coords[i]
                a_radius = atom_radii[i]
                r_ag = math.sqrt((ax - gx)**2 + (ay - gy)**2 + (az - gz)**2)
                f_ag = 1.0 - 0.5 * (math.erf(zeta_g * (a_radius - r_ag)) +
                                    math.erf(zeta_g * (a_radius + r_ag)))

                sw_func[g] *= f_ag

        return sw_func

    def alt_switching_function(self, molecule, grid):
        grid = self.qgrid
        atom_radii = molecule.vdw_radii_to_numpy() * 1.2
        atom_coords = molecule.get_coordinates_in_bohr()
        #assert_msg_critical(
        #    atom_coords.shape[0] == atom_radii.shape[0],
        #    'CosmoDriver.get_switching_function: Inconsistent atom_coords ' +
        #    'and atom_radii')

        #assert_msg_critical(
        #    grid.shape[0] == self.grid_per_sphere * atom_coords.shape[0],
        #    'CosmoDriver.get_switching_function: Should only be used on ' +
        #    'raw COSMO grid')

        sw_func = np.zeros(grid.shape[0])

        for g in range(grid.shape[0]):
            gx, gy, gz = grid[g, :3]
            zeta_g = grid[g, 4]

            sw_func[g] = 1.0
            atom_idx = g // self.grid_per_sphere

            for i in range(atom_coords.shape[0]):
                if i == atom_idx:
                    continue

                ax, ay, az = atom_coords[i]
                a_radius = atom_radii[i]
                r_ag = math.sqrt((ax - gx)**2 + (ay - gy)**2 + (az - gz)**2)
                f_ag = 1.0 - 0.5 * (math.erf(zeta_g * (a_radius - r_ag)) +
                                    math.erf(zeta_g * (a_radius + r_ag)))

                sw_func[g] *= f_ag

        return sw_func

    def get_zeta_dict(self):
        #ref: B. A. Gregersen and D. M. York, J. Chem. Phys. 122, 194110 (2005)
        
        return {
            50: 4.893,
            110: 4.901,
            194: 4.903,
            302: 4.905,
            434: 4.906,
            590: 4.905,
            770: 4.899,
            974: 4.907,
            2030: 4.907
        }

    def form_matrix_A(self, grid, sw_func):
        grid = self.qgrid

        Amat = np.zeros((grid.shape[0], grid.shape[0]))

        sqrt_2_invpi = math.sqrt(2.0 / math.pi)

        for i in range(grid.shape[0]):
            xi, yi, zi, wi, zeta_i, atom_idx = grid[i]
            Amat[i, i] = zeta_i * sqrt_2_invpi / sw_func[i]
            zeta_i2 = zeta_i**2

            for j in range(i + 1, grid.shape[0]):
                xj, yj, zj, wj, zeta_j, atom_idx = grid[j]

                zeta_j2 = zeta_j**2
                zeta_ij = zeta_i * zeta_j / math.sqrt(zeta_i2 + zeta_j2)

                r_ij = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)

                Aij = math.erf(zeta_ij * r_ij) / r_ij
                Amat[i, j] = Aij
                Amat[j, i] = Aij

        return Amat

    def form_matrix_B(self, grid, molecule):
        grid = self.qgrid

        Bmat = np.zeros((grid.shape[0], molecule.number_of_atoms()))
        natoms = molecule.number_of_atoms()
        atom_coords = molecule.get_coordinates_in_bohr()

        for i in range(grid.shape[0]):
            xi, yi, zi, wi, zeta_i, atom_idx = grid[i]
            for a in range(natoms):
                xa, ya, za = atom_coords[a]
                r_ia = math.sqrt((xi - xa)**2 + (yi - ya)**2 + (zi - za)**2)
                Bmat[i, a] = math.erf(zeta_i * r_ia) / r_ia

        return Bmat

    def form_vector_C(self, grid, molecule, basis, D, erf):
        grid = self.qgrid

        esp = np.zeros(grid.shape[0])
        # electrostatic potential integrals
        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        local_esp = np.zeros(end - start)
        npot_drv = NuclearPotentialDriver()
        nerf_drv = NuclearPotentialErfDriver()
        
        for i in range(start, end):
            if erf:
                epi_matrix = 1.0 * nerf_drv.compute(molecule, basis, [1.0], [grid[i,:3]], grid[i,4]).full_matrix().to_numpy()
            else:
                epi_matrix = 1.0 * npot_drv.compute(molecule, basis, [1.0], [grid[i,:3]]).full_matrix().to_numpy()

            if local_comm.Get_rank() == mpi_master():
                local_esp[i - start] -= np.sum(epi_matrix * D)

        if local_comm.Get_rank() == mpi_master():
            local_esp = cross_comm.gather(local_esp, root=mpi_master())

        if self.rank == mpi_master():
            for i in range(self.nodes):
                start = sum(counts[:i])
                end = sum(counts[:i + 1])
                esp[start:end] += local_esp[i]
            return esp
        else:
            return None

    def get_contribution_to_Fock(self, molecule, basis, grid, q, erf):
        grid = self.qgrid

        grid_coords = grid[:, :3].copy()
        zeta        = grid[:, 4].copy()

        if self.rank == mpi_master():
            if erf:
                nerf_drv = NuclearPotentialErfDriver()
                V_es = -1.0 * nerf_drv.compute(molecule, basis, q, grid_coords, zeta).full_matrix().to_numpy()
            else:
                npot_drv = NuclearPotentialDriver()
                V_es = -1.0 * npot_drv.compute(molecule, basis, q, grid_coords).full_matrix().to_numpy()
        else:
            V_es = None

        return V_es

    def visualize_cosmo_grid(self, molecule, grid):
        """
        Visualizes grid for surface discretization.

        :param molecule:
            The molecule.
        :param grid:
            The grid.
        """

        try:
            import py3Dmol as p3d
        except ImportError:
            raise ImportError('Unable to import py3Dmol.')

        assert_msg_critical(grid.shape[1] == 6,
                            'CosmoDriver.visualize_grid: Invalid grid size')

        grid_in_angstrom = grid[:, :3] * bohr_in_angstrom()

        grid_xyz_string = f'{grid_in_angstrom.shape[0]}\n\n'

        for i in range(grid_in_angstrom.shape[0]):
            x, y, z = grid_in_angstrom[i]
            grid_xyz_string += f'He {x} {y} {z}\n'

        v = p3d.view(width=400, height=400)

        v.addModel(molecule.get_xyz_string(), 'xyz')
        v.setStyle({'stick': {}})

        v.addModel(grid_xyz_string, 'xyz')
        v.setStyle({'elem': 'He'},
                   {'sphere': {
                       'radius': 0.05,
                       'color': 'red',
                       'opacity': 0.5
                   }})

        v.zoomTo()
        v.show()

    def cosmo_grad_contribution(self, molecule, basis, grid, sw_f, q, D, eps, x):
        # get the C-PCM gradient contribution
        transl_inv = True
        # Helper functions
        def dr_rij(ri, rj):
            r_ij = np.array([(ri[0] - rj[0]), (ri[1] - rj[1]), (ri[2] - rj[2])])
            return r_ij / np.linalg.norm(r_ij)

        def delta_ij(M, index_i, index_j):
            delta_iM, delta_jM = 0, 0
            if M == index_i:
                delta_iM = 1
            if M == index_j:
                delta_jM = 1

            return delta_iM - delta_jM
        
        def grad_B(molecule, grid, q):
            # Nuclear-cavity contribution
            atom_coords = molecule.get_coordinates_in_bohr()
            natoms = molecule.number_of_atoms()
            natoms_ = natoms
            two_sqrt_invpi = 2 / math.sqrt(math.pi)
            dB_mat = np.zeros((grid.shape[0], natoms, natoms, 3))

            if transl_inv:
                natoms -= 1
            for I in range(natoms):
                for m in range(grid.shape[0]):
                    xi, yi, zi, wi, zeta_i, atom_idx = grid[m]
                    for a in range(natoms_):
                        xa, ya, za = atom_coords[a]
                        ra = np.array([xa, ya, za])
                        #factor = delta_ij(I, atom_idx, -1)
                        factor = delta_ij(I, atom_idx, a)

                        ri = np.array([xi, yi, zi])
                        r_ia_2 = (xi - xa)**2 + (yi - ya)**2 + (zi - za)**2
                        r_ia = math.sqrt(r_ia_2)

                        dr_ia = factor * dr_rij(ri, ra)

                        dB_dr = -1.0 * (math.erf(zeta_i * r_ia) - two_sqrt_invpi * zeta_i * r_ia * 
                                            math.exp(-1.0 * zeta_i**2 * r_ia_2)) / r_ia_2
                            
                        dB_mat[m, a, I] = dB_dr * dr_ia
            if transl_inv:
                dB_mat[:,:,-1] = -np.sum(dB_mat[:,:,:-1], axis=2)

            gradBvec = np.einsum('m,mzax,z -> ax', q, dB_mat, molecule.get_element_ids())
            return gradBvec

        def grad_C(molecule, basis, grid, q, DM):
            natoms = molecule.number_of_atoms()
            grad_C_nuc = np.zeros((natoms, 3))
            grad_C_cav = np.zeros((natoms, 3))
            atom_indices = grid[:, -1]
            labels = ['X', 'Y', 'Z']
            geom100_drv = NuclearPotentialGeom100Driver()
            geom010_drv = NuclearPotentialGeom010Driver()

            natoms -= 1 # translational invariance
            for a in range(natoms):
                indices = (atom_indices == a)
                geom100_mats, geom010_mats, geom001_mats = [], [], []
                gra = grid[indices, :3]
                
                for i, charge in enumerate(q[indices]):
                    grad_010 = geom010_drv.compute(molecule, basis, [charge], [gra[i].tolist()])
                    temp_ = []
                    for label in labels:
                        temp_.append(-1.0 *grad_010.matrix(label).full_matrix().to_numpy())
                    geom010_mats.append(np.array(temp_))
                    #geom010_mats.append(np.array([-1.0 * grad_010.matrix_to_numpy(label) for label in labels]))

                grad_100 = geom100_drv.compute(basis, molecule, a, q, grid[:, :3])

                for label in labels:
                    geom100_mats.append(-1.0 *grad_100.matrix(label).full_matrix().to_numpy())
                    geom001_mats.append(-1.0 *grad_100.matrix(label).full_matrix().to_numpy().T)
                
                geom010_mats = np.array(geom010_mats)
                geom100_mats = np.array(geom100_mats)
                geom001_mats = np.array(geom001_mats)
                geom100_mats += geom001_mats

                grad_C_nuc[a] = np.einsum('kxij,ij -> x', geom010_mats, DM)
                grad_C_cav[a] = np.einsum("mn,xmn -> x", DM, geom100_mats)

            grad_C_cav[-1] = -np.sum(grad_C_cav[:-1], axis=0)
            grad_C_nuc[-1] = -np.sum(grad_C_nuc[:-1], axis=0)
            return grad_C_nuc + grad_C_cav

        def grad_Aij(molecule, grid, q, eps, x):
            two_sqrt_invpi = 2 / math.sqrt(math.pi)
            natoms         = molecule.number_of_atoms()
            scale_f        = -(eps - 1) / (eps + x)
            grad           = np.zeros((grid.shape[0], grid.shape[0], natoms, 3)) # (grid_pts, grid_pts, natoms, 3)
            grid_coords    = grid[:, :3]
            weights        = grid[:, 3]
            zeta           = grid[:, 4]
            atom_indices   = grid[:, 5]
            zeta_2         = zeta**2

            delta_r = grid_coords[:, None] - grid_coords[None, :]
            r_ij_2  = np.sum(delta_r**2, axis=-1) 
            r_ij    = np.sqrt(r_ij_2)
            np.fill_diagonal(r_ij, 1e-12), np.fill_diagonal(r_ij_2, 1e-12)
            dr_rij  = delta_r / r_ij[..., None]

            zeta_ij  = (zeta[:, None] * zeta[None, :]) / np.sqrt((zeta_2[:, None] + zeta_2[None, :]))
            delta_ij = np.array([(a == atom_indices[:, None]).astype(int) - (a == atom_indices[None, :]).astype(int) for a in range(natoms - 1)])
            dA_dr    = -1.0 * (erf(zeta_ij * r_ij) - two_sqrt_invpi * zeta_ij * r_ij * np.exp(-1.0 * zeta_ij**2 * r_ij_2)) / r_ij_2
            
            for a in range(natoms - 1):
                _delta      = delta_ij[a]
                dr_ij       = _delta[..., None] * dr_rij
                grad[:,:,a] = dA_dr[..., None] * dr_ij
            
            grad[:,:,-1] = -np.sum(grad[:,:,:-1], axis=2)
            gradAij = (-0.5 / scale_f) * np.einsum('i,ijax,j -> ax', q, grad, q, optimize=True)
            return gradAij

        def grad_Aii(molecule, grid, sw_f, q, eps, x):
            natoms = molecule.number_of_atoms()
            natoms_ = natoms
            grad = np.zeros((grid.shape[0], grid.shape[0], natoms, 3))
            sqrt_2_inv_pi = math.sqrt(2 / math.pi)
            atom_radii = molecule.vdw_radii_to_numpy() * 1.2
            atom_coords = molecule.get_coordinates_in_bohr()
            scale_f = -(eps - 1) / (eps + x)
            
            if transl_inv:
                natoms -= 1
            for a in range(natoms):
                for i in range(grid.shape[0]):
                    xi, yi, zi, wi, zeta_i, atom_idx = grid[i]
                    r_i = np.array([xi, yi, zi])
                    F_i = sw_f[i]
                    
                    summed_fi = np.zeros((1,3))
                    for J in range(natoms_):
                        r_J = atom_coords[J]
                        RJ = atom_radii[J]
                        vec_r_iJ = r_i - r_J
                        r_iJ = np.linalg.norm(vec_r_iJ)
                        factor = delta_ij(a, atom_idx, J)
                        
                        dr_iJ = factor * dr_rij(r_i, r_J) 
                        fiJ = 1 - 0.5 * (math.erf(zeta_i * (RJ - r_iJ)) + math.erf(zeta_i * (RJ + r_iJ)))
                        dfiJ_driJ = zeta_i/math.sqrt(math.pi) * (-1.0 * math.exp(-1.0 * zeta_i**2 * (RJ - r_iJ)**2) + math.exp(-1.0 * zeta_i**2 * (RJ + r_iJ)**2))
                        summed_fi += dfiJ_driJ / fiJ * dr_iJ
                    
                    grad_F_i = -1.0 * F_i * summed_fi
                    
                    grad[i,i,a] = -1.0 * zeta_i * sqrt_2_inv_pi / F_i**2 * grad_F_i
            if transl_inv:
                grad[:,:,-1] = -np.sum(grad[:,:,:-1], axis=2)

            grad_Aii = (-0.5 / scale_f) * np.einsum('i, ijax, j -> ax', q, grad, q)
            return grad_Aii

        gradA = grad_Aij(molecule, grid, q, eps, x) + grad_Aii(molecule, grid, sw_f, q, eps, x)
        gradB = grad_B(molecule, grid, q)
        gradC = grad_C(molecule, basis, grid, q, D)
        return gradA + gradB + gradC
    

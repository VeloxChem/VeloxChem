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
from .veloxchemlib import NuclearPotentialErfDriver
from .veloxchemlib import NuclearPotentialErfGeom010Driver
from .veloxchemlib import NuclearPotentialErfGeom100Driver
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .inputparser import parse_input, print_keywords


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
        - epsilon: The dielectric constant of the solvent.
        - x: Parameter used in the (denominator of) scaling function.
    """

    def __init__(self, comm=None, ostream=None):
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

        # model settings
        # standard value for dielectric const. is for that of water
        self.epsilon         = 78.39 
        self.grid_per_sphere = 194
        self.x               = 0
        
        # input keywords
        self.input_keywords = {
            'cpcm': {
                'grid_per_sphere': ('int', 'number of Lebedev grid points'),
                'epsilon': ('float', 'dielectric constant of solvent'),
                'x': ('float', 'parameter for scaling function'),
            },
        }

    def print_keywords(self):
        """
        Prints input keywords in cpcm driver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, cpcm_dict):
        """
        Updates settings in CPCM driver.

        :param cpcm_dict:
            The dictionary of CPCM input.
        """

        cpcm_keywords = {
            key: val[0] for key, val in self.input_keywords['cpcm'].items()
        }

        parse_input(self, cpcm_keywords, cpcm_dict)

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

    def generate_cpcm_grid(self, molecule):
        """
        Generates Lebedev grid for surface discretization.

        :param molecule:
            The molecule.

        :return:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        """

        valid_grid_numbers = [50, 110, 194, 302, 434, 590, 770, 974, 2030]

        assert_msg_critical(
            self.grid_per_sphere in valid_grid_numbers,
            'CpcmDriver.generate_grid: Invalid grid_per_sphere')

        unit_grid = gen_lebedev_grid(self.grid_per_sphere)
        unit_grid_coords = unit_grid[:, :3]
        unit_grid_weights = unit_grid[:, 3:]
        # standard normalization of lebedev weights -- unit sphere surface; 1 -> 4pi
        unit_grid_weights *= 4 * np.pi

        zeta = self.get_zeta_dict()[self.grid_per_sphere]

        # increase radii by 20%
        atom_radii = molecule.vdw_radii_to_numpy() * 1.2
        atom_coords = molecule.get_coordinates_in_bohr()

        cpcm_grid_raw = np.zeros((0, 6))

        for i in range(molecule.number_of_atoms()):
            # scale and shift unit grid
            atom_grid_coords = unit_grid_coords * atom_radii[i] + atom_coords[i]
            grid_zeta = zeta / (atom_radii[i] * np.sqrt(unit_grid_weights))
            atom_idx = np.full_like(grid_zeta, i)
            atom_grid = np.hstack(
                (atom_grid_coords, unit_grid_weights, grid_zeta, atom_idx))
            cpcm_grid_raw = np.vstack((cpcm_grid_raw, atom_grid))

        sw_func_raw = self.get_switching_function(atom_coords, atom_radii, cpcm_grid_raw)
        sw_mask = (sw_func_raw > 1.0e-8)

        cpcm_grid = cpcm_grid_raw[sw_mask, :]
        sw_func = sw_func_raw[sw_mask]

        return cpcm_grid, sw_func

    def get_switching_function(self, atom_coords, atom_radii, grid):
        """
        Construct the switching function.

        :param atom_coords:
            The atomic coordinates (a.u.)
        :param atom_radii:
            The Van der Waals atomic radii (a.u.)
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        
        :return:
            Switching function array of each grid point.
        """
        assert_msg_critical(
            atom_coords.shape[0] == atom_radii.shape[0],
            'CpcmDriver.get_switching_function: Inconsistent atom_coords ' +
            'and atom_radii')

        assert_msg_critical(
            grid.shape[0] == self.grid_per_sphere * atom_coords.shape[0],
            'CpcmDriver.get_switching_function: Should only be used on ' +
            'raw CPCM grid')

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
                r_ag = np.sqrt((ax - gx)**2 + (ay - gy)**2 + (az - gz)**2)
                f_ag = 1.0 - 0.5 * (math.erf(zeta_g * (a_radius - r_ag)) +
                                    math.erf(zeta_g * (a_radius + r_ag)))

                sw_func[g] *= f_ag

        return sw_func
    
    def get_zeta_dict(self):
        """
        Return the dictionary of Gaussian exponents for different grid-levels.
        
        Ref.: B. A. Gregersen and D. M. York, J. Chem. Phys. 122, 194110 (2005)
        """
        return {
            50: 4.89250673295,
            110: 4.90101060987,
            194: 4.90337644248,
            302: 4.90498088169,
            434: 4.90567349080,
            590: 4.90624071359,
            770: 4.90656435779,
            974: 4.90685167998,
            2030: 4.90744499142,
        }

    def form_matrix_A(self, grid, sw_func):
        """
        Forms the cavity-cavity interaction matrix.

        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param sw_func:
            The switching function.
        
        :return:
            The (cavity) electrostatic potential at each grid point due to
            the grid (unweighted by the charges).
        """

        Amat = np.zeros((grid.shape[0], grid.shape[0]))

        sqrt_2_invpi = np.sqrt(2.0 / np.pi)

        for i in range(grid.shape[0]):
            xi, yi, zi, wi, zeta_i, atom_idx = grid[i]
            Amat[i, i] = zeta_i * sqrt_2_invpi / sw_func[i]
            zeta_i2 = zeta_i**2

            for j in range(i + 1, grid.shape[0]):
                xj, yj, zj, wj, zeta_j, atom_idx = grid[j]

                zeta_j2 = zeta_j**2
                zeta_ij = zeta_i * zeta_j / np.sqrt(zeta_i2 + zeta_j2)

                r_ij = np.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)

                Aij = math.erf(zeta_ij * r_ij) / r_ij
                Amat[i, j] = Aij
                Amat[j, i] = Aij

        return Amat

    def form_matrix_B(self, grid, molecule):
        """
        Forms the nuclear-cavity interaction matrix.

        :param molecule:
            The molecule.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        
        :return:
            The (nuclear) electrostatic potential at each grid point due to
            each nucleus (not weighted by nuclear charge).
        """
        Bmat = np.zeros((grid.shape[0], molecule.number_of_atoms()))
        natoms = molecule.number_of_atoms()
        atom_coords = molecule.get_coordinates_in_bohr()

        for i in range(grid.shape[0]):
            xi, yi, zi, wi, zeta_i, atom_idx = grid[i]
            for a in range(natoms):
                xa, ya, za = atom_coords[a]
                r_ia = np.sqrt((xi - xa)**2 + (yi - ya)**2 + (zi - za)**2)
                Bmat[i, a] = math.erf(zeta_i * r_ia) / r_ia

        return Bmat

    def form_vector_C(self, molecule, basis, grid, D):
        """
        Forms the electron-cavity interaction vector.

        :param molecule:
            The molecule.
        :param basis:
            The atomic basis.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param D:
            The density matrix.

        :return:
            The total (electronic) electrostatic potential at each grid point.
        """
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
        nerf_drv = NuclearPotentialErfDriver()
        
        for i in range(start, end):
            epi_matrix = 1.0 * nerf_drv.compute(molecule, basis, [1.0], [grid[i,:3]], grid[i,4]).full_matrix().to_numpy()

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

    def get_contribution_to_Fock(self, molecule, basis, grid, q):

        grid_coords = grid[:, :3].copy()
        zeta        = grid[:, 4].copy()

        if self.rank == mpi_master():
            nerf_drv = NuclearPotentialErfDriver()
            V_es = -1.0 * nerf_drv.compute(molecule, basis, q, grid_coords, zeta).full_matrix().to_numpy()

        else:
            V_es = None

        return V_es

    def visualize_cpcm_grid(self, molecule, grid):
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
                            'CpcmDriver.visualize_grid: Invalid grid size')

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

    def grad_Aij(self, molecule, grid, q, eps, x):
        """
        Calculates the (off-diagonal) cavity-cavity gradient contribution.

        :param molecule: 
            The molecule object.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param q:
            The grid point charges.
        :param eps:
            Dielectric constant.
        :param x: 
            Alternative scale term in the denominator of
            the scaling function f.

        :return:
        The gradient array of each cartesian component -- of shape (nAtoms, 3).
        """
        # Defnine constants
        two_sqrt_invpi = 2.0 / np.sqrt(np.pi)
        natoms         = molecule.number_of_atoms()
        scale_f        = -(eps - 1) / (eps + x)
        avoid_div_0    = 1e-12
        grid_coords    = grid[:, :3]
        zeta           = grid[:, 4]
        atom_indices   = grid[:, 5]
        
        # M = Nr. of grid pts.
        M = grid_coords.shape[0]
        grad = np.zeros((M, M, natoms, 3))

        delta_r = grid_coords[:, None, :] - grid_coords[None, :, :]
        r_ij_2  = np.sum(delta_r**2, axis=-1)
        r_ij_2_avoid = np.where(r_ij_2 == 0.0, avoid_div_0, r_ij_2)
        r_ij        = np.sqrt(r_ij_2_avoid)
        dr_rij      = delta_r / r_ij[..., None]

        # Construct the explicit terms appearing in the gradient
        zeta_2   = zeta**2
        zeta_ij  = (zeta[:, None] * zeta[None, :]) / np.sqrt(zeta_2[:, None] + zeta_2[None, :])
        erf_term = erf(zeta_ij * r_ij)
        exp_term = np.exp(-zeta_ij**2 * r_ij_2_avoid)

        dA_dr = -1.0 * (erf_term - two_sqrt_invpi * zeta_ij * r_ij * exp_term) / r_ij_2_avoid
        
        # Definitions to keep track of which atom each piint belongs to
        I_vals      = np.arange(natoms - 1)[:, None, None]
        atom_idx_m  = atom_indices[None, :, None]
        atom_idx_n  = atom_indices[None, None, :]
        delta_ij    = ((I_vals == atom_idx_m).astype(int)
                    - (I_vals == atom_idx_n).astype(int))

        # Prepare for broadcasting
        _dA_dr = dA_dr[..., None, None]
        _delta_ij = delta_ij.transpose(1, 2, 0)[..., None]
        _dr_rij = dr_rij[..., None, :]

        # Construct the (n-1) gradient terms
        grad[:, :, :-1, :] = _dA_dr * _delta_ij * _dr_rij
        # Translational invariance
        grad[:, :, -1, :] = -np.sum(grad[:, :, :-1, :], axis=2)

        # Perform contractions
        #gradAij = (-0.5 / scale_f) * np.einsum('i,ijax,j->ax', q, grad, q, optimize=True)
        partial = np.tensordot(q, grad, axes=(0, 0))  
        gradAij = np.tensordot(partial, q, axes=(0, 0)) 
        gradAij *= (-0.5 / scale_f)

        return gradAij

    def grad_Aii(self, molecule, grid, sw_f, q, eps, x):
        """
        Calculates the (diagonal) cavity-cavity gradient contribution.

        :param molecule: 
            The molecule object.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param sw_f:
            The switching function.
        :param q:
            The grid point charges.
        :param eps:
            Dielectric constant.
        :param x: 
            Alternative scale factor in the denominator of
            the scaling function f.

        :return:
        The gradient array of each cartesian component -- of shape (nAtoms, 3).
        """
        # Basic constants
        sqrt_2_inv_pi = np.sqrt(2.0 / np.pi)
        scale_f       = -(eps - 1) / (eps + x)
        avoid_div_0   = 1e-12
        natoms        = molecule.number_of_atoms()
        atom_coords   = molecule.get_coordinates_in_bohr()
        atom_radii    = molecule.vdw_radii_to_numpy() * 1.2

        # Grid info
        M = grid.shape[0]
        grid_coords = grid[:, :3]
        zeta_i      = grid[:, 4]
        atom_idx    = grid[:, 5]
        F_i         = sw_f

        grad = np.zeros((M, M, natoms, 3))

        # Basic geometrical vectors
        r_iJ = grid_coords[:, None, :] - atom_coords[None, :, :]
        r_iJ_norm = np.linalg.norm(r_iJ, axis=2)

        # Avoiding division by zero
        r_iJ_safe = np.where(r_iJ_norm == 0.0, avoid_div_0, r_iJ_norm)
        dr_iJ     = r_iJ / r_iJ_safe[..., None]

        # Definitions to keep track of which atom each point belongs to
        a_vals = np.arange(natoms - 1)[:, None, None]
        _i_idx = atom_idx[None, :, None]
        J_vals = np.arange(natoms)[None, None, :]
        delta  = ((a_vals == _i_idx).astype(int)
                - (a_vals == J_vals).astype(int))
        
        # Prep for broadcasting
        RJ = atom_radii[None, :]
        _r_iJ = r_iJ_norm
        _zeta = zeta_i[:, None]

        # Compute fiJ
        term_m = _zeta * (RJ - _r_iJ)
        term_p = _zeta * (RJ + _r_iJ)
        fiJ = 1.0 - 0.5*(erf(term_m) + erf(term_p))

        # Derivative of fiJ: dfiJ_driJ
        z2 = _zeta**2
        dfiJ_driJ = (_zeta / np.sqrt(np.pi)) * (
            -np.exp(-z2 * (RJ - _r_iJ)**2) + np.exp(-z2 * (RJ + _r_iJ)**2))

        # Avoid division by zero
        fiJ_avoid = np.where(fiJ == 0.0, avoid_div_0, fiJ)
        ratio_fiJ = dfiJ_driJ / fiJ_avoid

        # Prep for broadcasting
        _delta = delta[..., None]
        _ratio = ratio_fiJ[None, :, :, None]
        _dr_iJ  = dr_iJ[None, ...]

        partial = _delta * _ratio * _dr_iJ
        summed_fi = np.sum(partial, axis=2)

        _F_i = F_i[None, :, None]
        __zeta = zeta_i[None, :]

        grad_Fi = -_F_i * summed_fi

        factor_i = -__zeta * sqrt_2_inv_pi / (_F_i[..., 0]**2)
        _factor_i = factor_i[..., None]

        final_contribution = grad_Fi * _factor_i

        idx = np.arange(M)
        for a in range(natoms - 1):
            contrib_a = final_contribution[a, :, :]
            # SEt the diagonal
            grad[idx, idx, a, :] = contrib_a

        # Translational invariance
        grad[:, :, -1, :] = -np.sum(grad[:, :, :-1, :], axis=2)

        # Contract
        #grad_Aii = (-0.5 / scale_f) * np.einsum('i, ijax, j -> ax', q, grad, q)
        partial = np.tensordot(q, grad, axes=(0, 0))
        grad_Aii = np.tensordot(partial, q, axes=(0, 0))
        grad_Aii *= (-0.5 / scale_f)
        return grad_Aii

    def grad_B(self, molecule, grid, q):
        """
        Calculates the nuclear-cavity gradient contribution.

        :param molecule:
            The molecule.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param q:
            The grid point charges.

        :return:
            The gradient array of each cartesian component -- of shape (nAtoms, 3).
        """
        atom_coords = molecule.get_coordinates_in_bohr()
        natoms      = molecule.number_of_atoms()
        grid_coords = grid[:, :3]
        zeta_i      = grid[:, 4]
        atom_idx    = grid[:, 5]

        # Define constants
        # M = Nr. of grid points
        M              = grid_coords.shape[0]
        two_sqrt_invpi = 2.0 / np.sqrt(np.pi)
        avoid_div_0    = 1e-12
        dB_mat         = np.zeros((M, natoms, natoms, 3))

        # distances and direction vects.
        r_iA   = grid_coords[:, None, :] - atom_coords[None, :, :]
        r_iA_2 = np.sum(r_iA**2, axis=2)
        d_iA   = np.sqrt(r_iA_2)

        d_iA_avoid = np.where(d_iA== 0.0, avoid_div_0, d_iA)
        dr_iA      = r_iA / d_iA_avoid[..., None]

        zeta_r      = zeta_i[:, None] *d_iA
        erf_term    = erf(zeta_r)
        exp_term    = np.exp(-1.0 * (zeta_i[:, None]**2) * r_iA_2)
        denominator = np.where(r_iA_2 == 0.0, avoid_div_0, r_iA_2)

        dB_dr = -1.0 * (erf_term - two_sqrt_invpi * zeta_r * exp_term) / denominator

        # Set up for broadcasting
        I_vals    = np.arange(natoms - 1)[:, None, None]
        A_vals    = np.arange(natoms)[None, None, :]
        _atom_idx = atom_idx[None, :, None]
        factor    = ((I_vals == _atom_idx).astype(int)
                    - (I_vals == A_vals).astype(int))
        
        _dB_dr  = dB_dr[None, ..., None]
        _dr_iA  = dr_iA[None, ...]
        _factor = factor[..., None]
        _dB_mat = _dB_dr * _factor * _dr_iA

        # Reordeer to (M, nAtoms, nAtoms-1, 3)
        _dB_mat = np.transpose(_dB_mat, (1, 2, 0, 3))
        
        # Translational invariance
        dB_mat[:, :, :-1, :] = _dB_mat
        dB_mat[:, :, -1, :] = -np.sum(dB_mat[:, :, :-1, :], axis=2)
        
        # Contract
        #gradBvec = np.einsum('m,mzax,z -> ax', q, dB_mat, molecule.get_element_ids())
        partial_res = np.tensordot(q, dB_mat, axes=(0, 0))  # -> (z, ax)
        gradBvec = np.tensordot(partial_res, molecule.get_element_ids(), axes=(0, 0))  # -> (a, x)

        return gradBvec
    
    def grad_C(self, molecule, basis, grid, q, DM):
        """
        Calculates the electron-cavity and electron-nuclear gradient contribution.

        :param molecule: 
            The molecule object.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        :param q:
            The grid point charges.
        :param DM:
            The converged density matrix.

        :return:
        The gradient array of each cartesian component -- of shape (nAtoms, 3).
        """
        geom100_drv = NuclearPotentialErfGeom100Driver()
        geom010_drv = NuclearPotentialErfGeom010Driver()
        
        # DEfine constants
        natoms = molecule.number_of_atoms()
        grad_C_nuc = np.zeros((natoms, 3))
        grad_C_cav = np.zeros((natoms, 3))
        grid_coords  = grid[:, :3]
        zeta         = grid[:, 4]
        atom_indices = grid[:, 5].astype(int)
        labels = ['X', 'Y', 'Z']

        # Compute both the nuclear and cavity contributions
        for a in range(natoms - 1):
            # Indices where the grid belongs to atom a
            indices_a = (atom_indices == a)
            q_subset  = q[indices_a]
            grid_a    = grid_coords[indices_a]
            zeta_a    = zeta[indices_a]
            
            geom100_mats, geom010_mats, geom001_mats = [], [], []

            for i, charge in enumerate(q_subset):
                grad_010 = geom010_drv.compute(molecule, basis, [charge], [grid_a[i]], [zeta_a[i]])
                geom010_mats.append(np.array([-1.0 * grad_010.matrix(label).full_matrix().to_numpy() for label in labels]))

            grad_100 = geom100_drv.compute(molecule, basis, a, grid_coords, q, zeta)
            
            for label in labels:
                mat_100 = -1.0 * grad_100.matrix(label).full_matrix().to_numpy()
                geom100_mats.append(mat_100)
                geom001_mats.append(mat_100.T)

            geom100_mats = np.array(geom100_mats)
            geom010_mats = np.array(geom010_mats)
            geom001_mats = np.array(geom001_mats)
            geom100_mats += geom001_mats
            
            partial_nuc = np.tensordot(geom010_mats, DM, axes=([2, 3], [0, 1]))
            grad_C_nuc[a] = np.sum(partial_nuc, axis=0)
            #grad_C_nuc[a] = np.einsum('kxij,ij -> x', geom010_mats, DM)
            grad_C_cav[a] = np.tensordot(DM, geom100_mats, axes=([0, 1], [1, 2]))
            #grad_C_cav[a] = np.einsum('xij,ij-> x', geom100_mats, DM)

        # Translational invariance
        grad_C_cav[-1] = -np.sum(grad_C_cav[:-1], axis=0)
        grad_C_nuc[-1] = -np.sum(grad_C_nuc[:-1], axis=0)
        return grad_C_nuc + grad_C_cav

    def cpcm_grad_contribution(self, molecule, basis, grid, sw_f, q, D):
        """
        Collects the CPCM gradient contribution.
        """
        gradA = self.grad_Aij(molecule, grid, q, self.epsilon, self.x) + self.grad_Aii(molecule, grid, sw_f, q, self.epsilon, self.x)
        gradB = self.grad_B(molecule, grid, q)
        gradC = self.grad_C(molecule, basis, grid, q, D)
        return gradA + gradB + gradC

    def non_vec_cpcm_grad_contribution(self, molecule, basis, grid, sw_f, q, D, eps, x):
        # get the C-PCM gradient contribution (non-vectorized)

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
            two_sqrt_invpi = 2 / np.sqrt(np.pi)
            dB_mat = np.zeros((grid.shape[0], natoms, natoms, 3))

            for I in range(natoms - 1):
                for m in range(grid.shape[0]):
                    xi, yi, zi, wi, zeta_i, atom_idx = grid[m]
                    for a in range(natoms_):
                        xa, ya, za = atom_coords[a]
                        ra = np.array([xa, ya, za])
                        #factor = delta_ij(I, atom_idx, -1)
                        factor = delta_ij(I, atom_idx, a)

                        ri = np.array([xi, yi, zi])
                        r_ia_2 = (xi - xa)**2 + (yi - ya)**2 + (zi - za)**2
                        r_ia = np.sqrt(r_ia_2)

                        dr_ia = factor * dr_rij(ri, ra)

                        dB_dr = -1.0 * (math.erf(zeta_i * r_ia) - two_sqrt_invpi * zeta_i * r_ia * 
                                            math.exp(-1.0 * zeta_i**2 * r_ia_2)) / r_ia_2
                            
                        dB_mat[m, a, I] = dB_dr * dr_ia
            # transl inv
            dB_mat[:,:,-1] = -np.sum(dB_mat[:,:,:-1], axis=2)

            partial_res = np.tensordot(q, dB_mat, axes=(0, 0))  # -> (z, ax)
            gradBvec = np.tensordot(partial_res, molecule.get_element_ids(), axes=(0, 0))
            return gradBvec

        def grad_C(molecule, basis, grid, q, DM):
            natoms = molecule.number_of_atoms()
            grad_C_nuc = np.zeros((natoms, 3))
            grad_C_cav = np.zeros((natoms, 3))
            grid_coords    = grid[:, :3]
            weights        = grid[:, 3]
            zeta           = grid[:, 4]
            atom_indices   = grid[:, 5]
            labels = ['X', 'Y', 'Z']
            geom100_drv = NuclearPotentialErfGeom100Driver()
            geom010_drv = NuclearPotentialErfGeom010Driver()

            # -1: translational invariance
            for a in range(natoms - 1):
                indices = (atom_indices == a)
                geom100_mats, geom010_mats, geom001_mats = [], [], []
                grid_a = grid_coords[indices, :3]
                zeta_a = zeta[indices]
                
                for i, charge in enumerate(q[indices]):
                    grad_010 = geom010_drv.compute(molecule, basis, [charge], [grid_a[i].tolist()], [zeta_a[i]])
                    geom010_mats.append(np.array([-1.0 * grad_010.matrix(label).full_matrix().to_numpy() for label in labels]))

                grad_100 = geom100_drv.compute(molecule, basis, a, grid_coords, q, zeta)

                for label in labels:
                    geom100_mats.append(-1.0 *grad_100.matrix(label).full_matrix().to_numpy())
                    geom001_mats.append(-1.0 *grad_100.matrix(label).full_matrix().to_numpy().T)
                
                geom010_mats = np.array(geom010_mats)
                geom100_mats = np.array(geom100_mats)
                geom001_mats = np.array(geom001_mats)
                geom100_mats += geom001_mats

                partial_nuc = np.tensordot(geom010_mats, DM, axes=([2, 3], [0, 1]))
                grad_C_nuc[a] = np.sum(partial_nuc, axis=0)
                #grad_C_nuc[a] = np.einsum('kxij,ij -> x', geom010_mats, DM)
                grad_C_cav[a] = np.tensordot(DM, geom100_mats, axes=([0, 1], [1, 2]))
                #grad_C_cav[a] = np.einsum('xij,ij-> x', geom100_mats, DM)

            grad_C_cav[-1] = -np.sum(grad_C_cav[:-1], axis=0)
            grad_C_nuc[-1] = -np.sum(grad_C_nuc[:-1], axis=0)
            return grad_C_nuc + grad_C_cav

        def grad_Aij(molecule, grid, q, eps, x):
            two_sqrt_invpi = 2 / np.sqrt(np.pi)
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
            #gradAij = (-0.5 / scale_f) * np.einsum('i,ijax,j -> ax', q, grad, q, optimize=True)
            partial = np.tensordot(q, grad, axes=(0, 0))  
            gradAij = np.tensordot(partial, q, axes=(0, 0)) 
            gradAij *= (-0.5 / scale_f)
            return gradAij

        def grad_Aii(molecule, grid, sw_f, q, eps, x):
            natoms = molecule.number_of_atoms()
            natoms_ = natoms
            grad = np.zeros((grid.shape[0], grid.shape[0], natoms, 3))
            sqrt_2_inv_pi = np.sqrt(2 / np.pi)
            atom_radii = molecule.vdw_radii_to_numpy() * 1.2
            atom_coords = molecule.get_coordinates_in_bohr()
            scale_f = -(eps - 1) / (eps + x)
            
            for a in range(natoms - 1):
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
                        dfiJ_driJ = zeta_i/np.sqrt(np.pi) * (-1.0 * math.exp(-1.0 * zeta_i**2 * (RJ - r_iJ)**2) + math.exp(-1.0 * zeta_i**2 * (RJ + r_iJ)**2))
                        summed_fi += dfiJ_driJ / fiJ * dr_iJ
                    
                    grad_F_i = -1.0 * F_i * summed_fi
                    
                    grad[i,i,a] = -1.0 * zeta_i * sqrt_2_inv_pi / F_i**2 * grad_F_i
            # transl inv
            grad[:,:,-1] = - np.sum(grad[:,:,:-1], axis=2)

            #grad_Aii = (-0.5 / scale_f) * np.einsum('i, ijax, j -> ax', q, grad, q)
            partial = np.tensordot(q, grad, axes=(0, 0))
            grad_Aii = np.tensordot(partial, q, axes=(0, 0))
            grad_Aii *= (-0.5 / scale_f)
            return grad_Aii

        gradA = grad_Aij(molecule, grid, q, eps, x) + grad_Aii(molecule, grid, sw_f, q, eps, x)
        gradB = grad_B(molecule, grid, q)
        gradC = grad_C(molecule, basis, grid, q, D)
        return gradA + gradB + gradC
    

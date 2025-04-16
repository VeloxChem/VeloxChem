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

from mpi4py import MPI
import numpy as np
import math
import time
import sys

from .veloxchemlib import cpcm_form_matrix_A, cpcm_comp_grad_Aij
from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import compute_nuclear_potential_erf_values
from .veloxchemlib import compute_nuclear_potential_erf_gradient
from .veloxchemlib import NuclearPotentialErfDriver
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .inputparser import parse_input, print_keywords

try:
    import scipy
except ImportError:
    pass


class CpcmDriver:
    """
    Implements the CPCM driver.

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

        self.custom_vdw_radii = None
        
        # input keywords
        self.input_keywords = {
            'cpcm': {
                'grid_per_sphere': ('int', 'number of Lebedev grid points'),
                'epsilon': ('float', 'dielectric constant of solvent'),
                'x': ('float', 'parameter for scaling function'),
            },
        }

    @staticmethod
    def erf_array(array):
        """
        Computes erf of an array.
        """

        if 'scipy' in sys.modules:
            return scipy.special.erf(array)
        else:
            # slow alternative in case scipy is not available
            return np.vectorize(math.erf)(array)

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
    
    def get_cpcm_vdw_radii(self, molecule):
        """
        Get C-PCM VDW radii.

        :param molecule:
            The molecule.

        :return:
            The VDW radii of the atoms.
        """
 
        atom_radii = molecule.vdw_radii_to_numpy()
        if self.custom_vdw_radii is not None:
            assert_msg_critical(
                len(self.custom_vdw_radii) % 2 == 0,
                'C-PCM: expecting even number of entries for user-defined C-PCM radii')
 
            keys = self.custom_vdw_radii[0::2]
            vals = self.custom_vdw_radii[1::2]
 
            for key, val in zip(keys, vals):
                val_au = float(val) / bohr_in_angstrom()
                try:
                    idx = int(key) - 1
                    assert_msg_critical(
                        0 <= idx and idx < molecule.number_of_atoms(),
                        'C-PCM: invalid atom index for user-defined C-PCM radii')
                    atom_radii[idx] = val_au
                    self.ostream.print_info(
                        f'Applying user-defined C-PCM radius {val} for atom {key}')

                except ValueError:
                    elem_found = False
                    for idx, label in enumerate(molecule.get_labels()):
                        if label.upper() == key.upper():
                            atom_radii[idx] = val_au
                            elem_found = True
 
                    if elem_found:
                        self.ostream.print_info(
                            f'Applying user-defined C-PCM radius {val} for atom {key}')

            self.ostream.print_blank()

        return atom_radii

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
            'CpcmDriver.generate_cpcm_grid: Invalid grid_per_sphere')

        unit_grid = gen_lebedev_grid(self.grid_per_sphere)
        unit_grid_coords = unit_grid[:, :3]
        unit_grid_weights = unit_grid[:, 3:]

        # standard normalization of lebedev weights -- unit sphere surface; 1 -> 4pi
        # Ref.: B. A. Gregersen and D. M. York, J. Chem. Phys. 122, 194110 (2005)
        unit_grid_weights *= 4 * np.pi

        zeta = self.get_zeta_dict()[self.grid_per_sphere]

        # increase radii by 20%
        atom_radii = self.get_cpcm_vdw_radii(molecule) * 1.2
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

        return cpcm_form_matrix_A(grid, sw_func)

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
            The AO density matrix.

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

        local_esp = compute_nuclear_potential_erf_values(
            molecule, basis, np.copy(grid[start:end, :3]), D, np.copy(grid[start:end, 4]))

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

        natoms = molecule.number_of_atoms()
        scale_f = -(eps - 1.0) / (eps + x)

        grid_coords = np.copy(grid[:, :3])
        zeta = np.copy(grid[:, 4])
        atom_indices = np.copy(grid[:, 5].astype(int))

        grad_Aij = cpcm_comp_grad_Aij(grid_coords, zeta, atom_indices, q, natoms)
        grad_Aij *= (-0.5 / scale_f)

        return grad_Aij

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

        t0 = time.time()

        # Basic constants
        sqrt_2_inv_pi = np.sqrt(2.0 / np.pi)
        scale_f       = -(eps - 1) / (eps + x)
        natoms        = molecule.number_of_atoms()
        atom_coords   = molecule.get_coordinates_in_bohr()
        atom_radii    = self.get_cpcm_vdw_radii(molecule) * 1.2
        
        # Grid info
        M           = grid.shape[0]
        grid_coords = grid[:, :3]
        zeta_i      = grid[:, 4]
        atom_idx    = grid[:, 5]
        n_dim       = 3 # x,y,z
        F_i         = sw_f

        # Initialise grad. object
        grad = np.zeros((M, M, natoms, 3))

        t1 = time.time()

        # Basic geometrical vectors
        r_iJ      = grid_coords[:, np.newaxis, :] - atom_coords[np.newaxis, :, :]
        r_iJ_norm = np.linalg.norm(r_iJ, axis=2)
        dr_iJ     = r_iJ / r_iJ_norm[:, :, np.newaxis]
        RJ        = atom_radii[np.newaxis, :]

        t2 = time.time()

        # Definitions to keep track of which atom each grid point belongs to
        a_vals = np.arange(natoms)[:, np.newaxis, np.newaxis]
        _i_idx = atom_idx[np.newaxis, :, np.newaxis]
        J_vals = np.arange(natoms)[np.newaxis, np.newaxis, :]
        delta  = ((a_vals == _i_idx).astype(int)
                - (a_vals == J_vals).astype(int))

        t3 = time.time()

        # Compute fiJ
        term_m = zeta_i[:, np.newaxis] * (RJ - r_iJ_norm)
        term_p = zeta_i[:, np.newaxis] * (RJ + r_iJ_norm)
        fiJ    = 1.0 - 0.5*(self.erf_array(term_m) + self.erf_array(term_p))

        t4 = time.time()

        # Derivative of fiJ: dfiJ_driJ
        z2 = zeta_i[:, np.newaxis]**2
        dfiJ_driJ = (zeta_i[:, np.newaxis] / np.sqrt(np.pi)) * (
            -np.exp(-z2 * (RJ - r_iJ_norm)**2) + np.exp(-z2 * (RJ + r_iJ_norm)**2))
        ratio_fiJ = dfiJ_driJ / fiJ

        t5 = time.time()

        # Primitive switching func. derivative contribution
        f_i = delta[:, :, :, np.newaxis] * ratio_fiJ[np.newaxis, :, :, np.newaxis] * dr_iJ[np.newaxis, :, :, :]
        summed_fi = np.sum(f_i, axis=2)

        t6 = time.time()

        # Switching func. contribution
        grad_Fi = -F_i[np.newaxis, :, np.newaxis] * summed_fi

        factor_i = -zeta_i[np.newaxis, :] * sqrt_2_inv_pi / (F_i[np.newaxis, :]**2)

        final_contribution = grad_Fi * factor_i[:, :, np.newaxis]

        t7 = time.time()

        idx = np.arange(M)
        for a in range(natoms):
            contrib_a = final_contribution[a, :, :]
            # Set the diagonal
            grad[idx, idx, a, :] = contrib_a

        t8 = time.time()

        # Perform contractions
        partial_contract = np.matmul(q, grad.reshape(M, M * natoms * n_dim)).reshape(M, natoms * n_dim)
        grad_Aii         = np.matmul(partial_contract.T, q).reshape(natoms, n_dim)
        grad_Aii        *= (-0.5 / scale_f)

        t9 = time.time()

        self.ostream.print_info(f'Grad Aii step 1: {t1 - t0:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 2: {t2 - t1:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 3: {t3 - t2:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 4: {t4 - t3:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 5: {t5 - t4:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 6: {t6 - t5:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 7: {t7 - t6:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 8: {t8 - t7:.2f} sec')
        self.ostream.print_info(f'Grad Aii step 9: {t9 - t8:.2f} sec')
        self.ostream.print_blank()
        self.ostream.flush()

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
        n_dim       = 3 # x,y,z

        # Define constants
        # M = Nr. of grid points
        M              = grid_coords.shape[0]
        two_sqrt_invpi = 2.0 / np.sqrt(np.pi)
        dB_mat         = np.zeros((M, natoms, natoms, n_dim))

        # distances and direction vects.
        r_iA   = grid_coords[:, np.newaxis, :] - atom_coords[np.newaxis, :, :]
        r_iA_2 = np.sum(r_iA**2, axis=2)
        d_iA   = np.sqrt(r_iA_2)
        dr_iA  = r_iA / d_iA[:, :, np.newaxis]

        zeta_r   = zeta_i[:, np.newaxis] * d_iA
        erf_term = self.erf_array(zeta_r)
        exp_term = np.exp(-1.0 * (zeta_i[:, np.newaxis]**2) * r_iA_2)
        dB_dr    = -1.0 * (erf_term - two_sqrt_invpi * zeta_r * exp_term) / r_iA_2

        # Set up for broadcasting
        I_vals       = np.arange(natoms)[:, np.newaxis, np.newaxis]
        A_vals       = np.arange(natoms)[np.newaxis, np.newaxis, :]
        atom_idx_exp = atom_idx[np.newaxis, :, np.newaxis]
        factor       = ((I_vals == atom_idx_exp).astype(int)
                    - (I_vals == A_vals).astype(int))
        
        # Calculate for n atoms
        dB_mat = dB_dr[np.newaxis, :, :, np.newaxis] * factor[:, :, :, np.newaxis] * dr_iA[np.newaxis, :, :, :]
        # Reordeer to (M, nAtoms, nAtoms, 3)
        dB_mat = np.transpose(dB_mat, (1, 2, 0, 3))
        
        # Contract
        partial_contract = np.matmul(q, dB_mat.reshape(M, natoms * natoms * n_dim)).reshape(natoms, natoms * n_dim)
        gradB_vec = np.matmul(partial_contract.T, molecule.get_element_ids()).reshape(natoms, n_dim)

        return gradB_vec
    
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
        
        grid_coords  = np.copy(grid[:, :3])
        zeta         = np.copy(grid[:, 4])
        atom_indices = np.copy(grid[:, 5].astype(int))

        # TODO: parallelize over MPI ranks

        return compute_nuclear_potential_erf_gradient(
            molecule, basis, grid_coords, q, DM, zeta, atom_indices)

    def cpcm_grad_contribution(self, molecule, basis, grid, sw_f, q, D):
        """
        Collects the CPCM gradient contribution.
        """

        t0 = time.time()

        gradA = self.grad_Aij(molecule, grid, q, self.epsilon, self.x) + self.grad_Aii(molecule, grid, sw_f, q, self.epsilon, self.x)

        t1 = time.time()

        gradB = self.grad_B(molecule, grid, q)

        t2 = time.time()

        gradC = self.grad_C(molecule, basis, grid, q, D)

        t3 = time.time()

        self.ostream.print_info(f'Time spent in CPCM gradient: {t3 - t0:.2f} sec')
        self.ostream.print_info(f'    Grad A took {t1 - t0:.2f} sec')
        self.ostream.print_info(f'    Grad B took {t2 - t1:.2f} sec')
        self.ostream.print_info(f'    Grad C took {t3 - t2:.2f} sec')
        self.ostream.print_blank()
        self.ostream.flush()

        return gradA + gradB + gradC

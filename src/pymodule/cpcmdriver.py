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

from .veloxchemlib import cpcm_local_matrix_A_diagonals
from .veloxchemlib import cpcm_local_matrix_A_dot_vector
from .veloxchemlib import cpcm_comp_grad_Aii, cpcm_comp_grad_Aij
from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import compute_nuclear_potential_erf_values
from .veloxchemlib import compute_nuclear_potential_erf_gradient
from .veloxchemlib import NuclearPotentialErfDriver
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
        self.grid_per_sphere = (194, 110)
        self.x               = 0

        self.custom_vdw_radii = None
        
        # input keywords
        self.input_keywords = {
            'cpcm': {
                'grid_per_sphere':
                    ('seq_fixed_int', 'number of Lebedev grid points'),
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
            (len(self.grid_per_sphere) == 2) and
                (self.grid_per_sphere[0] in valid_grid_numbers) and
                (self.grid_per_sphere[1] in valid_grid_numbers),
            'CpcmDriver.generate_cpcm_grid: Invalid grid_per_sphere')

        unit_grid = gen_lebedev_grid(self.grid_per_sphere[0])
        unit_grid_coords = unit_grid[:, :3]
        unit_grid_weights = unit_grid[:, 3:]

        unit_hydrogen_grid = gen_lebedev_grid(self.grid_per_sphere[1])
        unit_hydrogen_grid_coords = unit_hydrogen_grid[:, :3]
        unit_hydrogen_grid_weights = unit_hydrogen_grid[:, 3:]

        # standard normalization of lebedev weights -- unit sphere surface; 1 -> 4pi
        # Ref.: B. A. Gregersen and D. M. York, J. Chem. Phys. 122, 194110 (2005)
        unit_grid_weights *= 4 * np.pi
        unit_hydrogen_grid_weights *= 4 * np.pi

        zeta = self.get_zeta_dict()[self.grid_per_sphere[0]]
        hydrogen_zeta = self.get_zeta_dict()[self.grid_per_sphere[1]]

        # increase radii by 20%
        atom_radii = self.get_cpcm_vdw_radii(molecule) * 1.2
        atom_coords = molecule.get_coordinates_in_bohr()
        identifiers = molecule.get_identifiers()

        cpcm_grid_raw = np.zeros((0, 6))

        natoms = molecule.number_of_atoms()
        ave, rem = divmod(natoms, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        atom_start = sum(counts[:self.rank])
        atom_end = sum(counts[:self.rank + 1])

        for i in range(atom_start, atom_end):
            if identifiers[i] == 1:
                # scale and shift unit grid of hydrogen atom
                atom_grid_coords = unit_hydrogen_grid_coords * atom_radii[i] + atom_coords[i]
                grid_zeta = hydrogen_zeta / (atom_radii[i] * np.sqrt(unit_hydrogen_grid_weights))
                atom_idx = np.full_like(grid_zeta, i)
                atom_grid = np.hstack(
                    (atom_grid_coords, unit_hydrogen_grid_weights, grid_zeta, atom_idx))
            else:
                # scale and shift unit grid of non-hydrogen atom
                atom_grid_coords = unit_grid_coords * atom_radii[i] + atom_coords[i]
                grid_zeta = zeta / (atom_radii[i] * np.sqrt(unit_grid_weights))
                atom_idx = np.full_like(grid_zeta, i)
                atom_grid = np.hstack(
                    (atom_grid_coords, unit_grid_weights, grid_zeta, atom_idx))
            cpcm_grid_raw = np.vstack((cpcm_grid_raw, atom_grid))

        gathered_cpcm_grid_raw = self.comm.allgather(cpcm_grid_raw)
        cpcm_grid_raw = np.vstack(gathered_cpcm_grid_raw)

        sw_func_raw = self.get_switching_function(atom_coords, atom_radii, cpcm_grid_raw)

        gathered_sw_func_raw = self.comm.allgather(sw_func_raw)
        sw_func_raw = np.hstack(gathered_sw_func_raw)

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

        npoints = grid.shape[0]
        ave, rem = divmod(npoints, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        sw_func = np.zeros(end - start)

        for g in range(start, end):
            gx, gy, gz = grid[g, :3]
            zeta_g = grid[g, 4]

            sw_func[g - start] = 1.0

            # TODO: save grid atom_idx in another array
            atom_idx = int(grid[g, 5])

            # TODO: consider moving to C++
            for i in range(atom_coords.shape[0]):
                if i == atom_idx:
                    continue

                ax, ay, az = atom_coords[i]
                a_radius = atom_radii[i]
                r_ag = np.sqrt((ax - gx)**2 + (ay - gy)**2 + (az - gz)**2)
                f_ag = 1.0 - 0.5 * (math.erf(zeta_g * (a_radius - r_ag)) +
                                    math.erf(zeta_g * (a_radius + r_ag)))

                sw_func[g - start] *= f_ag

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

    def form_local_precond(self, grid, sw_func):
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

        ave, rem = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        local_Adiag = cpcm_local_matrix_A_diagonals(grid, sw_func, start, end)
        local_precond = 1.0 / local_Adiag

        return local_precond

    def form_vector_Bz(self, grid, molecule):
        """
        Forms the nuclear-cavity interaction vector.

        :param molecule:
            The molecule.
        :param grid:
            The grid object containing the grid positions, weights, 
            the Gaussian exponents, and indices for which atom they belong to.
        
        :return:
            The (nuclear) electrostatic potential at each grid point due to
            each nucleus (not weighted by nuclear charge).
        """

        Bzvec = np.zeros(grid.shape[0])

        atom_coords = molecule.get_coordinates_in_bohr()
        elem_ids = molecule.get_element_ids()

        grid_coords = np.copy(grid[:, :3])
        grid_zeta = np.copy(grid[:, 4])

        for a in range(molecule.number_of_atoms()):

            r_ia = np.sqrt(np.sum((grid_coords - atom_coords[a])**2, axis=1))

            Bzvec += elem_ids[a] * self.erf_array(grid_zeta * r_ia) / r_ia

        return Bzvec

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

        ave, res = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        local_esp = compute_nuclear_potential_erf_values(
            molecule, basis, np.copy(grid[start:end, :3]), D,
            np.copy(grid[start:end, 4]))

        gathered_esp = self.comm.gather(local_esp, root=mpi_master())
        if self.rank == mpi_master():
            esp = np.hstack(gathered_esp)
        else:
            esp = None

        return esp

    def get_contribution_to_Fock(self, molecule, basis, grid, q):

        ave, res = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        grid_coords = grid[start:end, :3].copy()
        zeta        = grid[start:end, 4].copy()

        nerf_drv = NuclearPotentialErfDriver()

        local_V_es = -1.0 * nerf_drv.compute(molecule, basis, q[start:end], grid_coords, zeta).to_numpy()

        V_es = self.comm.reduce(local_V_es, root=mpi_master())

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

        npoints = grid_coords.shape[0]
        ave, rem = divmod(npoints, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        grad_Aij = cpcm_comp_grad_Aij(grid_coords, zeta, atom_indices, q, start, end, natoms)
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

        scale_f       = -(eps - 1) / (eps + x)

        grid_coords = np.copy(grid[:, :3])
        zeta_i      = np.copy(grid[:, 4])
        atom_idx    = np.copy(grid[:, 5])

        atom_coords   = molecule.get_coordinates_in_bohr()
        atom_radii    = self.get_cpcm_vdw_radii(molecule) * 1.2

        npoints = grid_coords.shape[0]
        ave, rem = divmod(npoints, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        # a, b: atoms
        # i: grid points
        # c: Cartesian components
        # np.einsum('aib,ib,ibc,i->ac', delta, ratio_fiJ, dr_iJ, factor_i * q**2)
        grad_Aii = cpcm_comp_grad_Aii(grid_coords, zeta_i, sw_f, atom_idx, q,
                                      start, end, atom_coords, atom_radii)
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

        grid_coords = np.copy(grid[:, :3])
        zeta_i      = np.copy(grid[:, 4])
        atom_idx    = np.copy(grid[:, 5]).astype(int)

        two_sqrt_invpi = 2.0 / np.sqrt(np.pi)

        elem_ids = molecule.get_element_ids()
        gradB_vec = np.zeros((natoms, 3))

        # np.einsum('ia,bia,iac,i,a->bc', dB_dr, factor, dr_iA, q, Z)

        for a in list(range(natoms))[self.rank::self.nodes]:

            r_iA   = grid_coords - atom_coords[a]
            r_iA_2 = np.sum(r_iA**2, axis=1)
            d_iA   = np.sqrt(r_iA_2)

            dr_iA  = r_iA / d_iA.reshape(-1, 1)

            zeta_r   = zeta_i * d_iA
            erf_term = self.erf_array(zeta_r)
            exp_term = np.exp(-1.0 * (zeta_i**2) * r_iA_2)
            dB_dr    = -1.0 * (erf_term - two_sqrt_invpi * zeta_r * exp_term) / r_iA_2

            for b in range(natoms):

                factor = (atom_idx == b).astype(int) - int(a == b)

                gradB_vec[b] += np.dot(dB_dr * factor * q, dr_iA) * elem_ids[a]

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
        
        npoints = grid.shape[0]

        ave, rem = divmod(npoints, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        grid_coords  = np.copy(grid[start:end, :3])
        zeta         = np.copy(grid[start:end, 4])
        atom_indices = np.copy(grid[start:end, 5].astype(int))

        return compute_nuclear_potential_erf_gradient(
            molecule, basis, grid_coords, q[start:end], DM, zeta, atom_indices)

    def cpcm_grad_contribution(self, molecule, basis, grid, sw_f, q, D):
        """
        Collects the CPCM gradient contribution.
        """

        gradA = self.grad_Aij(molecule, grid, q, self.epsilon, self.x)

        gradA += self.grad_Aii(molecule, grid, sw_f, q, self.epsilon, self.x)

        gradB = self.grad_B(molecule, grid, q)

        gradC = self.grad_C(molecule, basis, grid, q, D)

        return gradA + gradB + gradC

    def cg_solve_parallel_direct(self,
                                 grid,
                                 sw_func,
                                 precond,
                                 rhs,
                                 x0=None,
                                 cg_thresh=1.0e-8):
        """
        Solves the C-PCM equations using conjugate gradient.
        """

        try:
            from scipy.sparse import linalg
        except ImportError:
            raise ImportError('Unable to import scipy. Please install scipy ' +
                              'via pip or conda.')

        def matvec(v):
            """
            Matrix-vector product
            """

            ave, res = divmod(grid.shape[0], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
            start = sum(counts[:self.rank])
            end = sum(counts[:self.rank + 1])

            local_v = cpcm_local_matrix_A_dot_vector(grid, sw_func, v, start, end)

            ret = self.comm.allgather(local_v)
            return np.hstack(ret)

        def precond_matvec(v):
            """
            Matrix-vector product for preconditioner using the
            inverse of the diagonal
            """

            return precond * v

        n = grid.shape[0]

        LinOp = linalg.LinearOperator((n, n), matvec=matvec)
        PrecondOp = linalg.LinearOperator((n, n), matvec=precond_matvec)

        b = rhs
        if x0 is None:
            x0 = np.zeros(rhs.shape)

        try:
            cg_solution, cg_conv = linalg.cg(A=LinOp,
                                             b=b,
                                             x0=x0,
                                             M=PrecondOp,
                                             rtol=cg_thresh,
                                             atol=0)
        except TypeError:
            # workaround for scipy < 1.11
            cg_solution, cg_conv = linalg.cg(A=LinOp,
                                             b=b,
                                             x0=x0,
                                             M=PrecondOp,
                                             tol=(cg_thresh * np.linalg.norm(b)),
                                             atol=0)


        assert_msg_critical(cg_conv == 0,
                            'C-PCM: conjugate gradient solver did not converge')

        cg_solution = self.comm.bcast(cg_solution, root=mpi_master())

        return cg_solution

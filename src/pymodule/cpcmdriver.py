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
import sys

from .veloxchemlib import cpcm_local_matrix_A_diagonals
from .veloxchemlib import cpcm_local_matrix_A_dot_vector
from .veloxchemlib import cpcm_comp_grad_Aii, cpcm_comp_grad_Aij
from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import gen_lebedev_grid
from .veloxchemlib import compute_nuclear_potential_erf_values
from .veloxchemlib import compute_nuclear_potential_erf_gradient
from .veloxchemlib import compute_nuclear_potential_erf_hessian_200
from .veloxchemlib import compute_nuclear_potential_erf_hessian_101
from .veloxchemlib import compute_nuclear_potential_erf_hessian_110
from .veloxchemlib import compute_nuclear_potential_erf_hessian_020
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
        self.epsilon = 78.39
        self.optical_epsilon = 1.777849
        self.grid_per_sphere = (194, 110)
        self.x = 0

        self.custom_vdw_radii = None
        self.custom_vdw_radii_verbose = True
        self.radii_scaling = 1.2

        # input keywords
        self.input_keywords = {
            'cpcm': {
                'grid_per_sphere':
                    ('seq_fixed_int', 'number of Lebedev grid points'),
                'epsilon': ('float', 'dielectric constant of solvent'),
                'x': ('float', 'parameter for scaling function'),
            },
        }

    def print_cpcm_info(self):
        """
        Print information and reference.
        """

        cpcm_info = 'Using C-PCM with the ISWIG discretization method.'
        self.ostream.print_info(cpcm_info)
        self.ostream.print_blank()
        iswig_ref = 'A. W. Lange, J. M. Herbert,'
        iswig_ref += ' J. Chem. Phys. 2010, 133, 244111.'
        self.ostream.print_reference(iswig_ref)
        self.ostream.print_blank()
        self.ostream.flush()

    def init(self, molecule, basis, do_nuclear=True):
        """
        Initialize the driver for energy calculations.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param do_nuclear:
            Flag to compute or not the nuclear contribution.
        """

        (self._cpcm_grid,
         self._cpcm_sw_func) = self.generate_cpcm_grid(molecule)

        cpcm_local_precond = self.form_local_precond(self._cpcm_grid,
                                                     self._cpcm_sw_func)

        self._cpcm_precond = self.comm.allgather(cpcm_local_precond)
        self._cpcm_precond = np.hstack(self._cpcm_precond)

        if do_nuclear:
            self._cpcm_Bzvec = self.form_vector_Bz(self._cpcm_grid, molecule,
                                                   basis)
        else:
            self._cpcm_Bzvec = None

        self._cpcm_q = None

    def compute_gs_fock(self, molecule, basis, density, cpcm_cg_thresh):
        """
        Compute the energy and Fock matrix contribution.

        :param molecule:
            The molecule.
        :param basis:
            The atomic basis.
        :param density:
            The density matrix.
        :param cpcm_cg_thresh:
            threshold for solving charges.

        :return:
            The energy and Fock matrix contributions
        """

        Cvec = self.form_vector_C(molecule, basis, self._cpcm_grid, density)

        if self.rank == mpi_master():
            scale_f = -(self.epsilon - 1) / (self.epsilon + self.x)
            rhs = scale_f * (self._cpcm_Bzvec + Cvec)
        else:
            rhs = None
        rhs = self.comm.bcast(rhs, root=mpi_master())

        # in case number of C-PCM grid points do not match between
        # cpcm_q and rhs, such as between previous and current
        # geometries during an optimization, reset cpcm_q

        if self._cpcm_q is not None and self._cpcm_q.size != rhs.size:
            self._cpcm_q = None

        self._cpcm_q = self.cg_solve_parallel_direct(self._cpcm_grid,
                                                     self._cpcm_sw_func,
                                                     self._cpcm_precond, rhs,
                                                     self._cpcm_q,
                                                     cpcm_cg_thresh)

        if self.rank == mpi_master():
            self.cpcm_epol = self.compute_solv_energy(self._cpcm_Bzvec, Cvec,
                                                      self._cpcm_q)
        else:
            self.cpcm_epol = None

        Fock_sol = self.get_contribution_to_Fock(molecule, basis,
                                                 self._cpcm_grid, self._cpcm_q)

        return self.cpcm_epol, Fock_sol

    def compute_response_fock(self, molecule, basis, density, cpcm_cg_thresh,
                              non_equilibrium):
        """
        Compute the Fock matrix contribution to response calculations.

        :param molecule:
            The molecule.
        :param basis:
            The atomic basis.
        :param density:
            The density matrix.
        :param cpcm_cg_thresh:
            threshold for solving charges.
        :param non_equilibrium:
            flag to activate the use of non-equilibrium epsilon.

        :return:
            The Fock matrix contribution.
        """

        Cvec = self.form_vector_C(molecule, basis, self._cpcm_grid, density)

        if self.rank == mpi_master():
            if non_equilibrium:
                scale_f = -(self.optical_epsilon - 1) / (self.optical_epsilon +
                                                         self.x)
            else:
                scale_f = -(self.epsilon - 1) / (self.epsilon + self.x)
            rhs = scale_f * (Cvec)
        else:
            rhs = None

        rhs = self.comm.bcast(rhs, root=mpi_master())

        cpcm_rsp_q = self.cg_solve_parallel_direct(self._cpcm_grid,
                                                   self._cpcm_sw_func,
                                                   self._cpcm_precond, rhs,
                                                   None, cpcm_cg_thresh)

        return self.get_contribution_to_Fock(molecule, basis, self._cpcm_grid,
                                             cpcm_rsp_q)

    def compute_gradient(self, molecule, basis, density):
        """
        Compute CPCM gradient contribution.

        :param molecule:
            The molecule.
        :param basis:
            The atomic basis.
        :param density:
            The density matrix.

        :return:
            The CPCM gradient contribution.
        """

        gradA = self.grad_Aij(molecule, self._cpcm_grid, self._cpcm_q,
                              self.epsilon, self.x)

        gradA += self.grad_Aii(molecule, self._cpcm_grid, self._cpcm_sw_func,
                               self._cpcm_q, self.epsilon, self.x)

        gradB = self.grad_B(molecule, basis, self._cpcm_grid, self._cpcm_q)

        gradC = self.grad_C(molecule, basis, self._cpcm_grid, self._cpcm_q,
                            density)

        return gradA + gradB + gradC

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

    def compute_solv_energy(self, Bzvec, Cvec, q):
        """
        Computes (electrostatic component of) C-PCM energy.
        TODO: add other components of the energy

        :param Bzvec:
            The nuclear potential on the grid.
        :param Cvec:
            The electronic potential on the grid.
        :param q:
            The surface charges.

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
                'C-PCM: expecting even number of entries for user-defined C-PCM radii'
            )

            keys = self.custom_vdw_radii[0::2]
            vals = self.custom_vdw_radii[1::2]

            for key, val in zip(keys, vals):
                val_au = float(val) / bohr_in_angstrom()
                try:
                    idx = int(key) - 1
                    assert_msg_critical(
                        0 <= idx and idx < molecule.number_of_atoms(),
                        'C-PCM: invalid atom index for user-defined C-PCM radii'
                    )
                    atom_radii[idx] = val_au
                    if self.custom_vdw_radii_verbose:
                        self.ostream.print_info(
                            f'Applying user-defined C-PCM radius {val} for atom {key}'
                        )

                except ValueError:
                    elem_found = False
                    for idx, label in enumerate(molecule.get_labels()):
                        if label.upper() == key.upper():
                            atom_radii[idx] = val_au
                            elem_found = True

                    if elem_found:
                        if self.custom_vdw_radii_verbose:
                            self.ostream.print_info(
                                f'Applying user-defined C-PCM radius {val} for atom {key}'
                            )

            if self.custom_vdw_radii_verbose:
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

        # increase radii (default radii_scaling is 1.2)
        atom_radii = self.get_cpcm_vdw_radii(molecule) * self.radii_scaling
        atom_coords = molecule.get_coordinates_in_bohr()
        identifiers = molecule.get_identifiers()

        cpcm_grid_raw = np.zeros((0, 6))

        natoms = molecule.number_of_atoms()
        ave, rem = divmod(natoms, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        atom_start = sum(counts[:self.rank])
        atom_end = sum(counts[:self.rank + 1])

        for i in range(atom_start, atom_end):
            if identifiers[i] == 0:
                # dummy atom
                continue

            elif identifiers[i] == 1:
                # scale and shift unit grid of hydrogen atom
                atom_grid_coords = (unit_hydrogen_grid_coords * atom_radii[i] +
                                    atom_coords[i])
                grid_zeta = hydrogen_zeta / (
                    atom_radii[i] * np.sqrt(unit_hydrogen_grid_weights))
                atom_idx = np.full_like(grid_zeta, i)
                atom_grid = np.hstack(
                    (atom_grid_coords, unit_hydrogen_grid_weights, grid_zeta,
                     atom_idx))

            else:
                # scale and shift unit grid of non-hydrogen atom
                atom_grid_coords = (unit_grid_coords * atom_radii[i] +
                                    atom_coords[i])
                grid_zeta = zeta / (atom_radii[i] * np.sqrt(unit_grid_weights))
                atom_idx = np.full_like(grid_zeta, i)
                atom_grid = np.hstack(
                    (atom_grid_coords, unit_grid_weights, grid_zeta, atom_idx))

            cpcm_grid_raw = np.vstack((cpcm_grid_raw, atom_grid))

        gathered_cpcm_grid_raw = self.comm.allgather(cpcm_grid_raw)
        cpcm_grid_raw = np.vstack(gathered_cpcm_grid_raw)

        sw_func_raw = self.get_switching_function(atom_coords, atom_radii,
                                                  cpcm_grid_raw)

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

    def form_vector_Bz(self, grid, molecule, basis):
        """
        Forms the nuclear-cavity interaction vector.

        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and indices for which atom they belong to.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The (nuclear) electrostatic potential at each grid point due to
            each nucleus (not weighted by nuclear charge).
        """

        Bzvec = np.zeros(grid.shape[0])

        atom_coords = molecule.get_coordinates_in_bohr()
        elem_ids = molecule.get_effective_nuclear_charges(basis)

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
            The density matrix.

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
        """
        Gets CPCM contribution to Fock matrix.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param grid:
            The grid.
        :param q:
            The surface charges.

        :return:
            The CPCM contribution to Fock matrix.
        """

        ave, res = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        grid_coords = grid[start:end, :3].copy()
        zeta = grid[start:end, 4].copy()

        nerf_drv = NuclearPotentialErfDriver()

        local_V_es = -1.0 * nerf_drv.compute(molecule, basis, q[start:end],
                                             grid_coords, zeta).to_numpy()

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

    def visualize_cpcm_charges(self, molecule, colorbar=True):
        """
        Visualizes CPCM charges.

        :param molecule:
            The molecule.
        """

        try:
            import py3Dmol as p3d
        except ImportError:
            raise ImportError("Unable to import py3Dmol.")
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("Unable to import matplotlib.")
        try:
            import matplotlib.colors as mcolors
        except ImportError:
            raise ImportError("Unable to import matplotlib.colors.")

        assert_msg_critical(
            self._cpcm_grid is not None,
            "C-PCM grid not available. Driver not initialized.")

        assert_msg_critical(
            self._cpcm_q is not None, "C-PCM charges not available.\n"
            "Run SCF with C-PCM first, then call this function")

        grid = self._cpcm_grid
        q = self._cpcm_q

        assert_msg_critical(grid.shape[1] == 6,
                            "CpcmDriver.visualize_grid: Invalid grid size")
        assert_msg_critical(
            len(q) == len(grid),
            "CpcmDriver.visualize_grid: Invalid q-vector size")

        grid_in_angstrom = grid[:, :3] * bohr_in_angstrom()

        grid_xyz_string = f"{grid_in_angstrom.shape[0]}\n\n"

        for i in range(grid_in_angstrom.shape[0]):
            x, y, z = grid_in_angstrom[i]
            grid_xyz_string += f"He {x} {y} {z}\n"

        q_absmax = np.max(np.abs(q))
        norm = mcolors.TwoSlopeNorm(vmin=-q_absmax, vcenter=0.0, vmax=q_absmax)
        v = p3d.view(width=600, height=600)

        v.addModel(molecule.get_xyz_string(), "xyz")
        v.setStyle({"stick": {}})
        cmap = plt.get_cmap("coolwarm_r")

        for k in range(grid_in_angstrom.shape[0]):
            x, y, z = grid_in_angstrom[k, :3]

            color = mcolors.to_hex(cmap(norm(q[k])))
            v.addSphere({
                "center": {
                    "x": x,
                    "y": y,
                    "z": z
                },
                "radius": 0.05,
                "color": color,
                "opacity": 0.9,
            })

        v.zoomTo()
        v.show()
        if colorbar:
            fig, ax = plt.subplots(figsize=(6, 1))
            fig.subplots_adjust(bottom=0.5)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            cb = plt.colorbar(sm, cax=ax, orientation="horizontal")
            cb.set_label("Charge (a.u)")

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

        grad_Aij = cpcm_comp_grad_Aij(grid_coords, zeta, atom_indices, q, start,
                                      end, natoms)
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

        scale_f = -(eps - 1) / (eps + x)

        grid_coords = np.copy(grid[:, :3])
        zeta_i = np.copy(grid[:, 4])
        atom_idx = np.copy(grid[:, 5])

        atom_coords = molecule.get_coordinates_in_bohr()
        atom_radii = self.get_cpcm_vdw_radii(molecule) * self.radii_scaling

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

    def grad_B(self, molecule, basis, grid, q):
        """
        Calculates the nuclear-cavity gradient contribution.

        :param molecule:
            The molecule.
        :param basis:
            The atomic basis.
        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and indices for which atom they belong to.
        :param q:
            The grid point charges.

        :return:
            The gradient array of each cartesian component -- of shape (nAtoms, 3).
        """

        atom_coords = molecule.get_coordinates_in_bohr()
        natoms = molecule.number_of_atoms()

        grid_coords = np.copy(grid[:, :3])
        zeta_i = np.copy(grid[:, 4])
        atom_idx = np.copy(grid[:, 5]).astype(int)

        two_sqrt_invpi = 2.0 / np.sqrt(np.pi)

        elem_ids = molecule.get_effective_nuclear_charges(basis)
        gradB_vec = np.zeros((natoms, 3))

        # np.einsum('ia,bia,iac,i,a->bc', dB_dr, factor, dr_iA, q, Z)

        for a in list(range(natoms))[self.rank::self.nodes]:

            r_iA = grid_coords - atom_coords[a]
            r_iA_2 = np.sum(r_iA**2, axis=1)
            d_iA = np.sqrt(r_iA_2)

            dr_iA = r_iA / d_iA.reshape(-1, 1)

            zeta_r = zeta_i * d_iA
            erf_term = self.erf_array(zeta_r)
            exp_term = np.exp(-1.0 * (zeta_i**2) * r_iA_2)
            dB_dr = -1.0 * (erf_term -
                            two_sqrt_invpi * zeta_r * exp_term) / r_iA_2

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

        grid_coords = np.copy(grid[start:end, :3])
        zeta = np.copy(grid[start:end, 4])
        atom_indices = np.copy(grid[start:end, 5].astype(int))

        return compute_nuclear_potential_erf_gradient(molecule, basis,
                                                      grid_coords, q[start:end],
                                                      DM, zeta, atom_indices)

    def hess_B(self, molecule, grid, q):
        """
        Calculates the nuclear-cavity hessian contribution.

        :param molecule:
            The molecule.
        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and the indices for which atom the grid
            point belongs to.
        :param q:
            The c-pcm charges.

        :return:
            The hessian array of each cartesian pair, shape (nAtoms*3, nAtoms*3)
        """

        atom_coords = molecule.get_coordinates_in_bohr()
        natms = molecule.number_of_atoms()
        two_sqrt_invpi = 2 / math.sqrt(math.pi)

        d2B_mat = np.zeros((grid.shape[0], natms, natms, natms, 3, 3))
        z = molecule.get_element_ids()

        for m, (xk, yk, zk, wk, zeta_k, atom_idx) in enumerate(grid):
            for a in range(natms):
                R_A = atom_coords[a]
                rk = np.array([xk, yk, zk])
                r_kA_vec = rk - R_A
                r_kA = np.linalg.norm(r_kA_vec)
                unit_kA = r_kA_vec / r_kA

                exp_term = np.exp(-(zeta_k**2) * r_kA**2)
                erf_term = math.erf(zeta_k * r_kA)

                T1 = -(erf_term -
                       two_sqrt_invpi * zeta_k * r_kA * exp_term) / r_kA**2
                T2 = (2 * erf_term / r_kA**3 -
                      4 * zeta_k / np.sqrt(np.pi) * exp_term / r_kA**2 -
                      4 * zeta_k**3 / np.sqrt(np.pi) * exp_term)

                for L in range(natms):
                    factor1 = int(L == atom_idx) - int(L == a)
                    for alpha in range(3):
                        for J in range(natms):
                            factor2 = int(J == atom_idx) - int(J == a)
                            for beta in range(3):

                                delta = 1.0 if alpha == beta else 0.0

                                unitf = T2 * unit_kA[alpha] * unit_kA[beta] + (
                                    T1 / r_kA) * (
                                        delta - unit_kA[alpha] * unit_kA[beta])

                                d2B_mat[m, a, L, J, alpha,
                                        beta] += (factor1 * factor2 * unitf)

        return (np.einsum("m,maijxy,a->ijxy", q, d2B_mat,
                          z).transpose(0, 2, 1,
                                       3).reshape(3 * natms,
                                                  3 * natms))  # change this

    def hess_C(self, molecule, basis, grid, q, D):
        """
        Calculates the electron-cavity and electron-nuclear hessian contribution.

        :param molecule:
            The molecule object.
        :param basis:
            The basis.
        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and indices for which atom the grid point belongs to.
        :param q:
            The grid point charges.
        :param D:
            The converged density matrix.

        :return:
            The hessian array of each cartesian component, shape (nAtoms *3, nAtoms*3).
        """

        point_coords = np.copy(grid[:, :3])
        pq = np.copy(q)
        omega = np.copy(grid[:, 4])
        idx = np.copy(grid[:, -1].astype(int))

        ana_200 = compute_nuclear_potential_erf_hessian_200(
            molecule, basis, point_coords, pq, D, omega)
        ana_101 = compute_nuclear_potential_erf_hessian_101(
            molecule, basis, point_coords, pq, D, omega)
        ana_020 = compute_nuclear_potential_erf_hessian_020(
            molecule, basis, point_coords, pq, D, omega, idx)
        ana_110 = compute_nuclear_potential_erf_hessian_110(
            molecule, basis, point_coords, pq, D, omega, idx)

        return ana_200 + ana_101 + ana_020 + ana_110

    def hessA_ii(self, molecule, grid, swf, q, eps=78.39):
        """
        Calculates the diagonal cavity-cavity hessian contribution.

        :param molecule:
            The molecule object.
        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and indices for which atom the grid point belongs to.
        :param sw_f:
            The switching function.
        :param q:
            The grid point charges.
        :param eps:
            Dielectric constant.

        :return:
            The hessian array of each cartesian component, shape (nAtoms*3, nAtoms*3).
        """

        natms = molecule.number_of_atoms()
        atom_coords = molecule.get_coordinates_in_bohr()
        grid_coords = grid[:, :3]
        vdw_R = molecule.vdw_radii_to_numpy() * 1.2
        hessian = np.zeros((natms, 3, natms, 3))
        scale_f = -(eps - 1) / (eps)

        for k in range(grid.shape[0]):
            atm = int(grid[k, -1])
            z = grid[k, 4]
            F_k = swf[k]

            sum_df_over_f = np.zeros((natms, 3))
            sum_d2f = np.zeros((natms, natms, 3, 3))

            for L in range(natms):

                V_L = vdw_R[L]
                vec_kL = grid_coords[k] - atom_coords[L]
                r_kL = np.linalg.norm(vec_kL)
                unitvec_kL = vec_kL / r_kL

                f_kL = 1 - 1 / 2 * (math.erf(z * (V_L + r_kL)) +
                                    math.erf(z * (V_L - r_kL)))

                df_kL_pref = (
                    z / np.sqrt(np.pi) *
                    (np.exp(-(z**2) *
                            (V_L - r_kL)**2) - np.exp(-(z**2) *
                                                      (V_L + r_kL)**2)))
                df_kL = df_kL_pref * unitvec_kL

                d2f_pref1 = (z**3 / np.sqrt(np.pi) *
                             (2 * (V_L + r_kL) * np.exp(-(z**2) *
                                                        (V_L + r_kL)**2) - 2 *
                              (V_L - r_kL) * np.exp(-(z**2) * (V_L - r_kL)**2)))

                d2f_pref2 = (
                    z / np.sqrt(np.pi) *
                    (np.exp(-(z**2) *
                            (V_L + r_kL)**2) - np.exp(-(z**2) *
                                                      (V_L - r_kL)**2)) / r_kL)

                d2f_unitvec = np.outer(unitvec_kL, unitvec_kL)

                d2f_delta_unitvec = np.eye(3) - d2f_unitvec

                d2f_kL = d2f_pref1 * d2f_unitvec + d2f_pref2 * d2f_delta_unitvec

                for A in range(natms):

                    deltaA = int(A == atm) - int(A == L)  # delta_AI - delta_AL

                    sum_df_over_f[A] += deltaA * df_kL / f_kL

                    for B in range(natms):
                        deltaB = int(B == atm) - int(
                            B == L)  # delta_BI - delta_BL

                        sum_d2f[A, B] += (
                            deltaA * deltaB *
                            (-np.outer(df_kL, df_kL) / f_kL**2 + d2f_kL / f_kL))

            dF = np.zeros((natms, 3))
            d2F = np.zeros((natms, natms, 3, 3))

            for A in range(natms):
                dF[A] = -F_k * sum_df_over_f[A]

            for A in range(natms):
                for B in range(natms):

                    d2F[A, B] = -F_k * (np.outer(
                        sum_df_over_f[A], sum_df_over_f[B]) + sum_d2f[A, B])

            for A in range(natms):
                for B in range(natms):

                    A_contrib = (z * np.sqrt(2 / np.pi) *
                                 (2 * np.outer(dF[A], dF[B]) / F_k**3 -
                                  d2F[A, B] / F_k**2))

                    hessian[A, :,
                            B, :] += (-0.5 / scale_f) * q[k]**2 * A_contrib

        return hessian.reshape(3 * natms, 3 * natms)

    def hessA_ij(self, molecule, grid, q, eps=78.39):
        """
        Calculates the off-diagonal cavity-cavity hessian contribution.

        :param molecule:
            The molecule object.
        :param grid:
            The grid object containing the grid positions, weights,
            the Gaussian exponents, and indices for which atom they belong to.
        :param q:
            The grid point charges.
        :param eps:
            Dielectric constant.

        :return:
            The hessian array of each cartesian component, shape (nAtoms*3, nAtoms*3).
        """

        natms = molecule.number_of_atoms()
        hess = np.zeros((3 * natms, 3 * natms))
        scale_f = -(eps - 1) / (eps)

        for m in range(grid.shape[0]):
            grid_coords_m = grid[m, :3].squeeze()
            z_m = grid[m, 4]
            z_m2 = np.dot(z_m, z_m)
            idx_m = int(grid[m, -1])
            q_m = q[m]

            for n in range(m + 1, grid.shape[0]):
                grid_coords_n = grid[n, :3]
                z_n = grid[n, 4]
                z_n2 = np.dot(z_n, z_n)
                idx_n = int(grid[n, -1])
                q_n = q[n]

                z_mn = (z_m * z_n) / np.sqrt(z_n2 + z_m2)
                r_mn_vec = grid_coords_m - grid_coords_n
                r_mn = np.linalg.norm(r_mn_vec)
                r_mn2 = np.dot(r_mn_vec, r_mn_vec)
                rmn_unit = r_mn_vec / r_mn

                O_mn = (-(math.erf(z_mn * r_mn) - 2 / np.sqrt(np.pi) * z_mn *
                          r_mn * math.exp(-(z_mn**2) * r_mn**2)) / r_mn2)
                P_mn = -2 * O_mn / r_mn - 4 * z_mn**3 / np.sqrt(np.pi) * np.exp(
                    -(z_mn**2) * r_mn**2)

                for L in range(natms):
                    fL = int(idx_m == L) - int(idx_n == L)
                    for x in range(3):
                        drL = fL * r_mn_vec[x] / r_mn

                        for J in range(natms):
                            fJ = int(idx_m == J) - int(idx_n == J)
                            for y in range(3):
                                drJ = fJ * r_mn_vec[y] / r_mn

                                delta = 1.0 if x == y else 0.0

                                d2r = (-(fL * fJ) *
                                       (delta - rmn_unit[x] * rmn_unit[y]) /
                                       r_mn)
                                d2A_kl = P_mn * drL * drJ - O_mn * d2r
                                hess[3 * L + x, 3 * J + y] += (
                                    (-0.5 / scale_f) * q_m * d2A_kl * q_n) * 2

        return hess

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

            local_v = cpcm_local_matrix_A_dot_vector(grid, sw_func, v, start,
                                                     end)

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
                                             tol=(cg_thresh *
                                                  np.linalg.norm(b)),
                                             atol=0)

        assert_msg_critical(
            cg_conv == 0, 'C-PCM: conjugate gradient solver did not converge')

        cg_solution = self.comm.bcast(cg_solution, root=mpi_master())

        return cg_solution

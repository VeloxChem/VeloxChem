#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
import sys

from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import mpi_master
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .inputparser import parse_input


class RespChargesDriver:
    """
    Implements ESP and RESP charges.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
        - n_layers: The number of layers of scaled van der Waals surfaces.
        - density: The density of grid points in points per square Angstrom.
        - constraints: The constraints for charges to be equal.
        - restrained_hydrogens: If hydrogens should be restrained or not.
        - weak_restraint: The strength of the restraint in first stage of RESP fit.
        - strong_restraint: The strength of the restraint in second stage of RESP fit.
        - max_iter: The maximum number of iterations of the RESP fit.
        - threshold: The convergence threshold of the RESP fit.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes RESP charges driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        # grid information
        self.n_layers = 4
        self.density = 1.0

        # resp fitting
        self.constraints = None
        self.restrained_hydrogens = False
        self.weak_restraint = 0.0005
        self.strong_restraint = 0.001
        self.max_iter = 50
        self.threshold = 1.0e-6

    def update_settings(self, resp_dict):
        """
        Updates settings in RESP charges driver.

        :param resp_dict:
            The input dictionary of RESP charges group.
        """

        # TODO: add method_dict

        resp_keywords = {
            'n_layers': 'int',
            'density': 'float',
            'restrained_hydrogens': 'bool',
            'weak_restraint': 'float',
            'strong_restraint': 'float',
            'max_iter': 'int',
            'threshold': 'float',
        }

        parse_input(self, resp_keywords, resp_dict)

        # TODO: improve the format of constraints

        if 'constraints' in resp_dict:
            self.constraints = []
            for i in list(resp_dict['constraints'].split(', ')):
                self.constraints.append(list(map(int, list(i.split(' ')))))

    def compute(self, molecule, basis, scf_tensors):
        """
        Computes RESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :return:
            The RESP charges.
        """

        # generate grid points and calculate QM ESP for all conformers
        grid = self.get_grid_points(molecule)
        esp = self.get_electrostatic_potential(grid, molecule, basis,
                                               scf_tensors)

        if self.rank == mpi_master():
            self.print_header(esp.size)
            self.print_resp_header()

            constr_1, constr_2 = self.generate_constraints(molecule)

            # first stage of resp fit
            self.print_resp_stage_header('first')
            q0 = np.zeros(molecule.number_of_atoms())
            q = self.optimize_charges(grid, esp, q0, constr_1,
                                      self.weak_restraint, molecule)

            # quality of fit
            cur_str = 'RRMS          :   {:9.6f}  '.format(
                self.get_rrms(grid, esp, q, molecule))
            self.ostream.print_header(cur_str)
            self.ostream.print_blank()

            # second stage of RESP fit
            if constr_2 == [-1] * molecule.number_of_atoms():
                cur_str = '*** No refitting in second stage needed.'
                self.ostream.print_header(cur_str.ljust(40))
            else:
                self.print_resp_stage_header('second')
                q = self.optimize_charges(grid, esp, q, constr_2,
                                          self.strong_restraint, molecule)

                # quality of fit
                cur_str = 'RRMS          :   {:9.6f}  '.format(
                    self.get_rrms(grid, esp, q, molecule))
                self.ostream.print_header(cur_str)

            self.ostream.print_blank()
            self.ostream.flush()

            return q
        else:
            return None

    def compute_esp_charges(self, molecule, basis, scf_tensors):
        """
        Computes Merz-Kollman ESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :return:
            The ESP charges.
        """

        # generate grid points and calculate QM ESP for all conformers
        grid = self.get_grid_points(molecule)
        esp = self.get_electrostatic_potential(grid, molecule, basis,
                                               scf_tensors)

        if self.rank == mpi_master():
            self.print_header(esp.size)
            self.print_esp_header()

            # generate and solve equation system (ESP fit)
            a, b = self.generate_equation_system(grid, esp, molecule)
            q = np.linalg.solve(a, b)

            n_atoms = molecule.number_of_atoms()

            for i in range(n_atoms):
                cur_str = '{:3d}     {:2s}     {:11.6f}    '.format(
                    i + 1,
                    molecule.get_labels()[i], q[i])
                self.ostream.print_header(cur_str)

            self.ostream.print_header(31 * '-')
            cur_str = 'Total Charge  : {:9.6f}   '.format(np.sum(q[:n_atoms]))
            self.ostream.print_header(cur_str)
            cur_str = 'RRMS          : {:9.6f}   '.format(
                self.get_rrms(grid, esp, q, molecule))
            self.ostream.print_header(cur_str)

            self.ostream.print_blank()
            self.ostream.flush()

            return q[:n_atoms]
        else:
            return None

    def optimize_charges(self, grid, esp, q0, constr, restraint_strength,
                         molecule):
        """
        Iterative procedure to optimize charges for a single stage of the RESP fit.

        :param grid:
            The grid points.
        :param esp:
            The QM ESP on grid points.
        :param q0:
            The intial charges.
        :param constr:
            The constraints of the charges.
        :param restraint_strength:
            The strength of the hyperbolic penalty function.
        :param molecule:
            The molecule.

        :return:
            The optimized charges.
        """

        n_atoms = molecule.number_of_atoms()

        # initializes equation system
        a, b = self.generate_equation_system(grid, esp, molecule, q0, constr)
        q_new = np.concatenate((q0, np.zeros(b.size - n_atoms)), axis=None)

        dq_norm = 1.0
        it = 0

        # iterative RESP fitting procedure
        while dq_norm > self.threshold and it < self.max_iter:

            it += 1
            q_old = q_new

            # add restraint to matrix a
            rstr = np.zeros(b.size)
            for i in range(n_atoms):
                if (not self.restrained_hydrogens and
                        molecule.get_labels()[i] == 'H'):
                    continue
                elif constr[i] == -1:
                    continue
                else:
                    # TODO: add documentation for 0.1
                    rstr[i] += restraint_strength / np.sqrt(q_old[i]**2 +
                                                            0.1**2)
            a_tot = a + np.diag(rstr)

            q_new = np.linalg.solve(a_tot, b)
            dq_norm = np.linalg.norm(q_new - q_old)

        if dq_norm > self.threshold:
            # TODO: terminate calculation if not converged
            cur_str = '*** Warning: Charge fitting is not converged!'
            self.ostream.print_header(cur_str.ljust(40))
        else:
            cur_str = f'*** Charge fitting converged in {it} iterations.'
            self.ostream.print_header(cur_str.ljust(40))
        self.ostream.print_blank()

        cur_str = '{} | {} | {} | {} | {}'.format('No.', 'Atom', 'Variation',
                                                  'Initial Charge',
                                                  'Final Charge')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(56 * '-')

        for i in range(n_atoms):
            cur_str = '{:3d}     {:2s}     {:3d}      {:15.6f}   {:12.6f}'.format(
                i + 1,
                molecule.get_labels()[i], constr[i], q0[i], q_new[i])
            self.ostream.print_header(cur_str)

        self.ostream.print_header(56 * '-')
        cur_str = 'Total Charge  :   {:9.6f}   '.format(np.sum(q_new[:n_atoms]))
        self.ostream.print_header(cur_str)

        return q_new[:n_atoms]

    def generate_equation_system(self,
                                 grid,
                                 esp,
                                 molecule,
                                 q0=None,
                                 constr=None):
        """
        Generates the equation system for charge fitting to QM ESP.

        :param grid:
            The grid points.
        :param esp:
            The QM ESP on grid points.
        :param molecule:
            The molecule.
        :param constr:
            The constraints of the charges.

        :return:
            The matrix a and vector b.
        """

        n_atoms = molecule.number_of_atoms()

        if q0 is None:
            q0 = np.zeros(n_atoms)
        if constr is None:
            constr = [0] * n_atoms

        # number of Lagrangian constraints
        n_c = 0
        for i in range(n_atoms):
            if constr[i] > 0:
                n_c += 1

        # calculate unrestrained matrix a and vector b
        a = np.zeros((n_atoms + 1 + n_c, n_atoms + 1 + n_c))
        b = np.zeros(n_atoms + 1 + n_c)

        # total charge constraint
        a[:n_atoms, n_atoms] = 1
        a[n_atoms, :n_atoms] = 1
        b[n_atoms] = molecule.get_charge()

        coords = molecule.get_coordinates()

        for i in range(n_atoms):
            if constr[i] != -1:
                for p in range(esp.size):
                    r_pi = np.linalg.norm(grid[p] - coords[i])
                    b[i] += esp[p] / r_pi
                    for j in range(n_atoms):
                        a[i, j] += 1.0 / (r_pi *
                                          np.linalg.norm(grid[p] - coords[j]))

        # equivalent charge constraints and frozen/constant charges
        j = 0
        for i in range(n_atoms):
            if constr[i] == -1:
                a[i, i] = 1
                a[i, n_atoms] = 0
                b[i] = q0[i]
            elif constr[i] > 0:
                a[n_atoms + 1 + j, constr[i] - 1] = 1
                a[n_atoms + 1 + j, i] = -1
                a[constr[i] - 1, n_atoms + 1 + j] = 1
                a[i, n_atoms + 1 + j] = -1
                j += 1

        return a, b

    def generate_constraints(self, molecule):
        """
        Determines Lagrangian constraints and frozen charges for both RESP fit stages.

        :param molecule:
            The molecule.

        return:
            Two lists (constr_1 and constr_2) for variations in the two stages of
            the RESP fitting procedure for which -1 tells a charge to be frozen,
            and x > 0 that the atom has the same charge like atom x.
        """

        n_atoms = molecule.number_of_atoms()
        coords = molecule.get_coordinates()
        labels = molecule.get_labels()
        covalent_radii = molecule.covalent_radii_to_numpy()

        connectivity = np.full((n_atoms, n_atoms), False, dtype=bool)

        # connect atoms with distances close to sum of covalent radii with
        # tolerance of 0.4 Angstrom
        tolerance = 0.4 / bohr_in_angstroms()
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r_ij = np.linalg.norm(coords[i] - coords[j])
                if abs(r_ij - covalent_radii[i] -
                       covalent_radii[j]) < tolerance:
                    connectivity[i, j] = True
                    connectivity[j, i] = True

        # identify C(sp3)-H bonds from connectivity
        # c_nh: number of H bonded to atom if atom is C(sp3), else 0
        # h_cbond: atom number of C(sp3) to which H is bonded
        c_nh = [0] * n_atoms
        h_cbond = [0] * n_atoms

        for i in range(n_atoms):
            if labels[i] == 'C' and np.sum(connectivity[i]) == 4:
                for j in range(n_atoms):
                    if labels[j] == 'H' and connectivity[i, j]:
                        c_nh[i] += 1
                        h_cbond[j] = i + 1

        constr_1 = [0] * n_atoms
        constr_2 = [-1] * n_atoms

        # automatic constraints (refitting methyl and methylene in the 2nd stage)
        if self.constraints is None:
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    constr_2[i] = 0
                    j = 0
                    for k in range(n_atoms):
                        if h_cbond[k] == (i + 1) and j == 0:
                            constr_2[k] = 0
                            j = k + 1
                        elif h_cbond[k] == (i + 1) and j > 0:
                            constr_2[k] = j
                            j = k + 1

        # use constraints from input
        else:
            # put all constraints to constr_1
            for equal_atoms in self.constraints:
                for i in range(len(equal_atoms) - 1):
                    constr_1[equal_atoms[i + 1] - 1] = equal_atoms[i]

            # move methyl and methylene constraints to constr_2
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    constr_2[i] = constr_1[i]
                    constr_1[i] = 0
                    for j in range(n_atoms):
                        if h_cbond[j] == (i + 1):
                            constr_2[j] = constr_1[j]
                            constr_1[j] = 0

        return constr_1, constr_2

    def get_grid_points(self, molecule):
        """
        Gets grid points in the solvent-accessible region.

        :param molecule:
            The molecule.

        :return:
            The coordinates of each grid point listed in an array.
        """

        grid = []
        coords = molecule.get_coordinates()

        for layer in range(self.n_layers):

            # MK radii with layer-dependent scaling factor
            scaling_factor = 1.4 + layer * 0.4 / np.sqrt(self.n_layers)
            r = scaling_factor * molecule.mk_radii_to_numpy()

            for atom in range(molecule.number_of_atoms()):
                # number of points fitting on the equator
                n_eq = int(2.0 * np.pi * r[atom] * np.sqrt(self.density) *
                           bohr_in_angstroms())

                # number of latitudes
                n_i = n_eq // 2

                for i in range(n_i + 1):
                    phi = np.pi * i / n_i

                    # number of points on the latitude
                    if i == 0 or i == n_i:
                        n_j = 1
                    else:
                        n_j = int(n_eq * np.sin(phi) + 1e-10)

                    # generate new point and check for overlap
                    for j in range(n_j):
                        theta = 2 * np.pi * j / n_j
                        point = coords[atom] + r[atom] * np.array([
                            np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)
                        ])

                        save = True
                        for atom2 in range(molecule.number_of_atoms()):
                            if atom != atom2 and np.linalg.norm(
                                    point - coords[atom2]) < r[atom2]:
                                save = False
                                break
                        if save:
                            grid.append(point)

        return np.array(grid)

    def get_electrostatic_potential(self, grid, molecule, basis, scf_tensors):
        """
        Gets the QM ESP on the grid points.

        :param grid:
            The grid points.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :return:
            The ESP at each grid point listed in an array.
        """

        if self.rank == mpi_master():
            D = scf_tensors['D_alpha'] + scf_tensors['D_beta']
        else:
            D = None
        D = self.comm.bcast(D, root=mpi_master())

        coords = molecule.get_coordinates()
        elem_ids = molecule.elem_ids_to_numpy()

        esp = np.zeros(grid.shape[0])

        # classical electrostatic potential
        if self.rank == mpi_master():
            for i in range(esp.size):
                for j in range(molecule.number_of_atoms()):
                    esp[i] += elem_ids[j] / np.linalg.norm(coords[j] - grid[i])

        # electrostatic potential integrals

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(esp.size, self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        epi_drv = NuclearPotentialIntegralsDriver(local_comm)

        for i in range(start, end):
            epi_matrix = epi_drv.compute(molecule, basis, np.array([1.0]),
                                         np.array([grid[i]]))
            if local_comm.Get_rank() == mpi_master():
                esp[i] -= np.sum(epi_matrix.to_numpy() * D)

        if local_comm.Get_rank() == mpi_master():
            esp = cross_comm.reduce(esp, op=MPI.SUM, root=mpi_master())

        if self.rank == mpi_master():
            return esp
        else:
            return None

    def get_rrms(self, grid, esp, q, molecule):
        """
        Gets the relative root-mean-square (RRMS) error.

        :param grid:
            The grid points.
        :param esp:
            The QM ESP on grid points.
        :param q:
            The atom-centered charges.
        :param molecule:
            The molecule.

        :return:
            The RRMS.
        """

        # sum of the squared error between QM ESP and partial charge ESP
        chi_square = 0.0

        # sum of the QM ESP
        norm = 0.0

        coords = molecule.get_coordinates()

        for i in range(esp.size):

            # partial charge ESP
            pot = 0.0

            for j in range(molecule.number_of_atoms()):
                pot += q[j] / np.linalg.norm(grid[i] - coords[j])

            chi_square += (esp[i] - pot)**2
            norm += esp[i]**2

        return np.sqrt(chi_square / norm)

    def print_header(self, n_points):
        """
        Prints header for the RESP charges driver.

        :param n_points:
            The number of grid points.
        """

        self.ostream.print_blank()
        self.ostream.print_header('RESP Charges Driver Setup')
        self.ostream.print_header(29 * '=')
        self.ostream.print_blank()

        str_width = 40
        cur_str = 'Number of Layers             :  ' + str(self.n_layers)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Points per square Angstrom   :  ' + str(self.density)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Total Number of Grid Points  :  ' + str(n_points)
        self.ostream.print_header(cur_str.ljust(str_width))

    def print_resp_header(self):
        """
        Prints header for RESP charges calculation.

        """

        self.ostream.print_blank()
        self.ostream.print_header('Restrained ESP Charges')
        self.ostream.print_header(24 * '-')
        self.ostream.print_blank()

        if (self.n_layers != 4 or self.density != 1.0 or
                self.weak_restraint != 0.0005 or
                self.strong_restraint != 0.001):
            cur_str = '*** Warning: Parameters for RESP fitting differ from '
            cur_str += 'optimal choice!'
            self.ostream.print_header(cur_str.ljust(40))
            self.ostream.print_blank()

    def print_resp_stage_header(self, stage):
        """
        Prints header for a stage of the RESP fit.

        :param stage:
            The stage of the RESP fit.
        """

        str_width = 40

        if stage.lower() == 'first':
            self.ostream.print_header('First Stage Fit')
            self.ostream.print_header(17 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength         :  ' + str(
                self.weak_restraint)
        elif stage.lower() == 'second':
            self.ostream.print_header('Second Stage Fit')
            self.ostream.print_header(18 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength         :  ' + str(
                self.strong_restraint)

        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Restrained Hydrogens       :  ' + str(
            self.restrained_hydrogens)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Iterations  :  ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold      :  ' + str(self.threshold)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

    def print_esp_header(self):
        """
        Prints header for ESP charges calculation.

        """

        self.ostream.print_blank()
        self.ostream.print_header('Merz-Kollman ESP Charges')
        self.ostream.print_header(26 * '-')
        self.ostream.print_blank()
        cur_str = '{}   {}         {}'.format('No.', 'Atom', 'Charge    ')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(31 * '-')

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
from .outputstream import OutputStream


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
        - ostream: The output stream.
        - number_layers: The number of layers of scaled van der Waals surfaces.
        - density: The density of grid points in points per Å².
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

        # outputstream
        self.ostream = ostream

        # grid information
        self.number_layers = 4
        self.density = 1.0

        # resp fitting
        self.constraints = None
        self.restrained_hydrogens = False
        self.weak_restraint = 0.0005
        self.strong_restraint = 0.001
        self.max_iter = 50
        self.threshold = 1.0e-6

    def update_settings(self, resp_charges_dict):
        """
        Updates settings in RESP charges driver.

        :param resp_charges_dict:
            The input dictionary of RESP charges group.
        """

        if 'number_layers' in resp_charges_dict:
            self.number_layers = int(resp_charges_dict['number_layers'])
        if 'density' in resp_charges_dict:
            self.density = float(resp_charges_dict['density'])

        if 'constraints' in resp_charges_dict:
            self.constraints = []
            for i in list(resp_charges_dict['constraints'].split(', ')):
                self.constraints.append(list(map(int, list(i.split(' ')))))
        if 'restrained_hydrogens' in resp_charges_dict:
            self.restrained_hydrogens = bool(resp_charges_dict['restrained_hydrogens'])
        if 'weak_restraint' in resp_charges_dict:
            self.weak_restraint = float(resp_charges_dict['weak_restraint'])
        if 'strong_restraint' in resp_charges_dict:
            self.strong_restraint = float(resp_charges_dict['strong_restraint'])
        if 'max_iter' in resp_charges_dict:
            self.max_iter = int(resp_charges_dict['max_iter'])
        if 'threshold' in resp_charges_dict:
            self.threshold = float(resp_charges_dict['threshold'])

    def compute(self, molecule, basis, scf_tensors, weights=[1]):
        """
        Computes RESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param weights:
            The weight factors of different conformers.
        """

        # for one conformation
        if not isinstance(molecule, list):
            molecule = [molecule]
            basis = [basis]
            scf_tensors = [scf_tensors]

        grid = []
        esp = []
        n_points = 0

        # generate grid points and calculate QM ESP for all conformers
        for i in range(len(molecule)):
            grid.append(self.get_grid_points(molecule[i]))
            esp.append(self.get_electrostatic_potential(
                grid[i], molecule[i], basis[i], scf_tensors[i]))
            n_points += esp[i].size

        self.print_header(len(molecule), n_points)

        vary1, vary2 = self.generate_constraints(molecule[0])

        # first stage of resp fit
        self.print_resp_stage_header(True) 
        q0 = np.zeros(molecule[0].number_of_atoms())
        q = self.optimize_charges(grid, esp, q0, vary1, self.weak_restraint, molecule, weights)
        self.print_resp_stage_results(True, molecule[0], q, vary1)
        self.print_fit_quality(self.get_rrms(grid, esp, q, molecule, weights))

        self.ostream.print_blank()

        if vary2 == [-1] * molecule[0].number_of_atoms():
            cur_str = '*** No refitting in second stage needed.'
            self.ostream.print_header(cur_str.ljust(40))
        else:
            # second stage of RESP fit
            self.print_resp_stage_header(False)
            q = self.optimize_charges(grid, esp, q, vary2, self.strong_restraint, molecule, weights)
            self.print_resp_stage_results(False, molecule[0], q, vary2)
            self.print_fit_quality(self.get_rrms(grid, esp, q, molecule, weights))

        self.ostream.print_blank()
        self.ostream.flush()

        return q

    def compute_esp_charges(self, molecule, basis, scf_tensors, weights=[1]):
        """
        Computes Merz-Kollman ESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param weights:
            The weight factors of different conformers.

        :return:
            The ESP charges.
        """

        # for one conformation
        if not isinstance(molecule, list):
            molecule = [molecule]
            basis = [basis]
            scf_tensors = [scf_tensors]

        grid = []
        esp = []
        n_points = 0

        # generate grid points and calculate QM ESP for all conformers
        for i in range(len(molecule)):
            grid.append(self.get_grid_points(molecule[i]))
            esp.append(self.get_electrostatic_potential(
                grid[i], molecule[i], basis[i], scf_tensors[i]))
            n_points += esp[i].size

        self.print_header(len(molecule), n_points)
        self.print_esp_header()

        # generate and solve equation system (ESP fit)
        a, b = self.generate_equation_system(grid, esp, molecule, weights)
        q = np.linalg.solve(a, b)

        n_atoms = molecule[0].number_of_atoms()

        # print results
        for i in range(n_atoms):
            cur_str = '{:3d}     {:2s}     {:11.6f}    '.format(
                i + 1, molecule[0].get_labels()[i], q[i])
            self.ostream.print_header(cur_str)

        self.ostream.print_header(31 * '-')

        q_tot = np.sum(q[:n_atoms])        
        cur_str = 'Total Charge  : {:9.6f}   '.format(
             abs(q_tot) if q_tot >= -5e-7 else q_tot)   # no -0.000000
        self.ostream.print_header(cur_str)
        self.print_fit_quality(self.get_rrms(grid, esp, q, molecule, weights))
        self.ostream.print_blank()
        self.ostream.flush()

        return q[:n_atoms]

    def optimize_charges(self, grid, esp, q0, vary, restraint_strength, molecule, weights=[1]):
        """
        Iterative procedure to optimize charges for a single stage of the RESP fit.

        :param grid:
            The grid points.
        :param esp:
            The QM ESP on grid points.
        :param q0:
            The intial charges.
        :param vary:
            The variations of the charges (constraints and constants).
        :param restraint_strength:
            The strength of the hyperbolic penalty function.
        :param molecule:
            The molecule.
        :param weights:
            The weight factors of different conformers.

        :return:
            The optimized charges.
        """

        n_atoms = molecule[0].number_of_atoms()

        # initializes equation system
        a, b = self.generate_equation_system(grid, esp, molecule, weights, q0, vary)
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
                if not self.restrained_hydrogens and molecule[0].get_labels()[i] == 'H':
                    continue
                elif vary[i] == -1:
                    continue
                else:
                    rstr[i] += restraint_strength / np.sqrt(q_old[i]**2 + 0.1**2)
            a_tot = a + np.diag(rstr)

            q_new = np.linalg.solve(a_tot, b)
            dq_norm = np.linalg.norm(q_new - q_old)

        self.ostream.print_blank()
        if dq_norm > self.threshold:
            cur_str = '*** Warning: Charge fitting is not converged!'
            self.ostream.print_header(cur_str.ljust(40))
        else:
            cur_str = '*** Charge fitting converged in {} iterations.'.format(it)
            self.ostream.print_header(cur_str.ljust(40))

        return q_new[:n_atoms]

    def generate_equation_system(self, grid, esp, molecule, weights=[1], q0=None, vary=None):
        """
        Generates the equation system for charge fitting to QM ESP.

        :param grid:
            The grid points.
        :param esp:
            The QM ESP on grid points.
        :param molecule:
            The molecule.
        :param weights:
            The weight factors of different conformers.
        :param vary:
            The variations of the charges (constraints and constants).

        :return:
            The matrix a and vector b.
        """

        n_atoms = molecule[0].number_of_atoms()

        if vary is None:
            vary = [0] * n_atoms

        # number of Lagrangian constraints
        n_c = 0
        for i in range(n_atoms):
            if vary[i] > 0:
                n_c += 1

        # calculate unrestrained matrix a and vector b
        a = np.zeros((n_atoms + 1 + n_c, n_atoms + 1 + n_c))
        b = np.zeros(n_atoms + 1 + n_c)

        # total charge constraint
        a[:n_atoms, n_atoms] = 1
        a[n_atoms, :n_atoms] = 1
        b[n_atoms] = molecule[0].get_charge()

        for m in range(len(molecule)):
            a_m = np.zeros((n_atoms + 1 + n_c, n_atoms + 1 + n_c))
            b_m = np.zeros(n_atoms + 1 + n_c)

            coords = molecule[m].get_coordinates()

            for i in range(n_atoms):
                if vary[i] != -1:
                    for p in range(esp[m].size):
                        r_pi = np.linalg.norm(grid[m][p] - coords[i])
                        b_m[i] += esp[m][p] / r_pi
                        for j in range(n_atoms):
                            a_m[i, j] += 1 / (r_pi * np.linalg.norm(grid[m][p] - coords[j]))

            a += weights[m] * a_m
            b += weights[m] * b_m

        # equivalent charge constraints and frozen/constant charges
        j = 0
        for i in range(n_atoms):
            if vary[i] == -1:
                a[i, i] = 1
                a[i, n_atoms] = 0
                b[i] = q0[i]
            elif vary[i] > 0:
                a[n_atoms + 1 + j, int(vary[i]-1)] = 1
                a[n_atoms + 1 + j, i] = -1
                a[int(vary[i]-1), n_atoms + 1 + j] = 1
                a[i, n_atoms + 1 + j] = -1
                j += 1

        return a, b

    def generate_constraints(self, molecule):
        """
        Determines Lagrangian constraints and frozen charges for both RESP fit stages.

        :param molecule:
            The molecule.

        return:
            Two lists (vary1 and vary2) for variations in the two stages of
            the RESP fitting procedure for which -1 tells a charge to be frozen,
            and x > 0 that the atom has the same charge like atom x.
        """

        n_atoms = molecule.number_of_atoms()
        coords = molecule.get_coordinates()
        labels = molecule.get_labels()
        covalent_radii = molecule.covalent_radii_to_numpy()

        connectivity = np.zeros((n_atoms, n_atoms), dtype=bool)

        # connect atoms with distances close to sum of covalent radii with tolerance of 0.4 Å
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r_ij = np.linalg.norm(coords[i] - coords[j])
                if abs(r_ij - covalent_radii[i] - covalent_radii[j]) < (0.4 / bohr_in_angstroms()):
                    connectivity[i, j] = connectivity[j, i] = True

        # identify C(sp3)-H bonds from connectivity
        c_nh = [0] * n_atoms   # number of H bonded to atom if atom is C(sp3), else 0
        h_cbond = [0] * n_atoms   # atom number of C(sp3) to which H is bonded

        for i in range(n_atoms):
            if labels[i] == 'C' and np.sum(connectivity[i]) == 4:
                for j in range(n_atoms):
                    if labels[j] == 'H' and connectivity[i, j]:
                        c_nh[i] += 1
                        h_cbond[j] = i + 1

        vary1 = [0] * n_atoms
        vary2 = [-1] * n_atoms

        # automatic constraints (refitting methyl and methylene in the 2nd stage)
        if self.constraints is None:
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    vary2[i] = 0
                    j = 0
                    for k in range(n_atoms):
                        if h_cbond[k] == (i + 1) and j == 0:
                            vary2[k] = 0
                            j = k + 1
                        elif h_cbond[k] == (i + 1) and j > 0:
                            vary2[k] = j
                            j = k + 1

        # use constraints from input
        else:
            # put all constraints to vary1
            for equal_atoms in self.constraints:
                for i in range(len(equal_atoms) - 1):
                    vary1[equal_atoms[i + 1] - 1] = equal_atoms[i]

            # move methyl and methylene constraints to vary2
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    vary2[i] = vary1[i]
                    vary1[i] = 0
                    for j in range(n_atoms):
                        if h_cbond[j] == (i + 1):
                            vary2[j] = vary1[j]
                            vary1[j] = 0

        return vary1, vary2

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

        for layer in range(self.number_layers):

            # MK radii with layer-dependent scaling factor
            r = (1.4 + layer * 0.4 / np.sqrt(self.number_layers)) * molecule.mk_radii_to_numpy()

            for atom in range(molecule.number_of_atoms()):
                # number of points fitting on the equator
                n_eq = int(2 * np.pi * r[atom] * np.sqrt(self.density) * bohr_in_angstroms())

                # number of latitudes
                n_i = int(n_eq / 2)

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
                        point = coords[atom] + r[atom] * np.array(
                            [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])

                        save = True
                        for atom2 in range(molecule.number_of_atoms()):
                            if atom != atom2 and np.linalg.norm(point - coords[atom2]) < r[atom2]:
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

        D = scf_tensors['D']
        coords = molecule.get_coordinates()
        elem_ids = molecule.elem_ids_to_numpy()

        esp = np.zeros(grid.shape[0])
        for i in range(esp.size):

            # classical electrostatic potential
            for j in range(molecule.number_of_atoms()):
                esp[i] += elem_ids[j] / np.linalg.norm(coords[j] - grid[i])

            # electrostatic potential integrals
            epi_drv = NuclearPotentialIntegralsDriver()
            epi_matrix = epi_drv.compute(
                molecule, basis, np.array([[1.0]]), np.array([grid[i]]))

            esp[i] -= np.sum(epi_matrix.to_numpy() * D)

        return esp

    def get_rrms(self, grid, esp, q, molecule, weights=[1]):
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
        :param weights:
            The weight factors of different conformers.

        :return:
            The RRMS.
        """

        # sum of the squared error between QM ESP and partial charge ESP
        chi_square = 0.0

        # sum of the QM ESP
        norm = 0.0

        for m in range(len(molecule)):

            coords = molecule[m].get_coordinates()

            for i in range(esp[m].size):

                # partial charge ESP
                pot = 0.0

                for j in range(molecule[m].number_of_atoms()):
                    pot += q[j] / np.linalg.norm(grid[m][i] - coords[j])

                chi_square += weights[m] * (esp[m][i] - pot)**2
                norm += weights[m] * esp[m][i]**2

        return np.sqrt(chi_square / norm)

    def print_header(self, n_conf, n_points):
        """
        Prints header for the RESP charges driver.

        :param n_conf:
            The number of conformers.
        :param n_points:
            The number of grid points.
        """

        self.ostream.print_blank()
        self.ostream.print_header('RESP Charges Driver Setup')
        self.ostream.print_header(29 * '=')
        self.ostream.print_blank()

        str_width = 40
        cur_str = 'Number of Conformers         :  ' + str(n_conf)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Number of Layers             :  ' + str(self.number_layers)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Points per Square Angstrom   :  ' + str(self.density)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Total Number of Grid Points  :  ' + str(n_points)
        self.ostream.print_header(cur_str.ljust(str_width))

    def print_resp_stage_header(self, stage):
        """
        Prints header for a stage of the RESP fit.

        :param stage:
            The stage of the RESP fit.
        """

        str_width = 40
        self.ostream.print_blank()

        if stage:
            if (self.number_layers != 4 or self.density != 1.0 or
                    self.weak_restraint != 0.0005 or self.strong_restraint != 0.001):
                cur_str = '*** Warning: Parameters for RESP fitting differ from recommended choice!'
                self.ostream.print_header(cur_str.ljust(str_width))
                self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.print_header('First Stage Fit')
            self.ostream.print_header(17 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength           :  ' + str(self.weak_restraint)
        else:
            self.ostream.print_header('Second Stage Fit')
            self.ostream.print_header(18 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength           :  ' + str(self.strong_restraint)

        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Restrained Hydrogens         :  ' + str(self.restrained_hydrogens)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Iterations    :  ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold (a.u.) :  ' + str(self.threshold)
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
        cur_str = '{}   {}      {}'.format(
            'No.', 'Atom', 'Charge (a.u.)')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(31 * '-')

    def print_resp_stage_results(self, stage, molecule, q, vary):
        """
        Prints results for a stage of the RESP fit.

        :param stage:
            The stage of the RESP fit.
        :param q:
            The final charges of the stage.
        :param vary:
            The variations and constraints for the stage.
        """

        if stage:
            str_variation = ''
            width = 44
        else:
            str_variation = 'Frozen |'
            width = 52

        self.ostream.print_blank()
        cur_str = '{} | {} | {} {} | {}'.format(
            'No.', 'Atom', str_variation, 'Constraints', 'Charges (a.u.)')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(width * '-')

        for i in range(q.size):
            constraint = ''
            if vary[i] < 0:
                str_variation = 'Yes     '
            elif vary[i] == 0:
                str_variation = ' No     '
            else:
                str_variation = ' No     '
                constraint = str(vary[i])

            if stage:
                str_variation = ''

            cur_str = '{:3d}     {:2s}     {}    {:3.6s}     {:12.6f}   '.format(
                i + 1, molecule.get_labels()[i], str_variation, constraint, q[i])
            self.ostream.print_header(cur_str)

        self.ostream.print_header(width * '-')

        q_tot = np.sum(q)
        cur_str = 'Total Charge  : {:9.6f}   '.format(
             abs(q_tot) if q_tot >= -5e-7 else q_tot)   # no -0.000000
        self.ostream.print_header(cur_str)

    def print_fit_quality(self, rrms):
        """
        Prints fit quality.

        :param rrms:
            The relative root-mean-square error.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Fit Quality')
        self.ostream.print_header('-' * 13)
        cur_str = 'Relative Root-Mean-Square Error  : {:9.6f}'.format(rrms)
        self.ostream.print_header(cur_str)

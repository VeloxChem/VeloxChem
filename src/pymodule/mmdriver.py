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
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom, hartree_in_kcalpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical, safe_arccos


class MMDriver:
    """
    Implements molecular mechanics driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes molecular mechanics driver.
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
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def generate_force_field(self, molecule):
        """
        Generate the force field for the given molecule.

        :param molecule:
            The molecule.
        """

        # NOTE: Never use this method inside ForceFieldGenerator.
        # Otherwise it will be a cyclic import.

        from .forcefieldgenerator import ForceFieldGenerator

        ff_generator = ForceFieldGenerator()
        ff_generator.ostream.mute()
        ff_generator.create_topology(molecule)

        self.load_force_field(ff_generator)

    def load_force_field(self, ff_generator):
        """
        Load the force field from a force field generator.

        :param ff_generator:
            The force field generator.
        """

        self.force_field = {}

        self.force_field['atoms'] = dict(ff_generator.atoms)
        self.force_field['bonds'] = dict(ff_generator.bonds)
        self.force_field['pairs'] = dict(ff_generator.pairs)

        self.force_field['fudgeLJ'] = ff_generator.fudgeLJ
        self.force_field['fudgeQQ'] = ff_generator.fudgeQQ

        self.force_field['angles'] = dict(ff_generator.angles)
        self.force_field['dihedrals'] = dict(ff_generator.dihedrals)
        self.force_field['impropers'] = dict(ff_generator.impropers)

    def compute_bonds(self, coordinates):
        """
        Compute the potential energy and forces for the bonds.

        :param coordinates:
            The coordinates of the atoms.
        """

        for (i, j), bond in self.force_field['bonds'].items():

            equilibrium = bond['equilibrium']
            force_constant = bond['force_constant']

            r_ij = coordinates[j] - coordinates[i]
            distance = np.linalg.norm(r_ij)
            n_ij = r_ij / distance

            self.stretch_potential += 0.5 * force_constant * (distance -
                                                              equilibrium)**2

            g = force_constant * (distance - equilibrium)

            self.gradient[i] += -g * n_ij
            self.gradient[j] += g * n_ij

            self.exclusions.append((i, j))

    def compute_angles(self, coordinates):
        """
        Compute the potential energy and forces for the angles.

        :param coordinates:
            The coordinates of the atoms.
        """

        deg_to_rad = np.pi / 180.0

        for (i, j, k), angle in self.force_field['angles'].items():

            # Note: equilibrium_angle in radian
            equilibrium_angle = angle['equilibrium'] * deg_to_rad
            force_constant = angle['force_constant']

            # Vectors from j to i and j to k
            r_ji = coordinates[i] - coordinates[j]
            r_kj = coordinates[j] - coordinates[k]

            # Distances
            distance_ji = np.linalg.norm(r_ji)
            distance_kj = np.linalg.norm(r_kj)

            # Normalized vectors
            n_ji = r_ji / distance_ji
            n_kj = r_kj / distance_kj

            # Cosine and sine of the current angle
            cos_theta = -np.dot(n_ji, n_kj)
            theta_ijk = safe_arccos(cos_theta)
            sin_theta = np.sin(theta_ijk)

            # simple workaround for straight angles
            sin_theta = max(sin_theta, 1e-8)

            self.bend_potential += 0.5 * force_constant * (theta_ijk -
                                                           equilibrium_angle)**2

            g = force_constant * (theta_ijk - equilibrium_angle)

            grad_i = (cos_theta * n_ji + n_kj) / (sin_theta * distance_ji) * g
            grad_k = -(cos_theta * n_kj + n_ji) / (sin_theta * distance_kj) * g

            self.gradient[i] += grad_i
            self.gradient[j] += -(grad_i + grad_k)
            self.gradient[k] += grad_k

            self.exclusions.append((i, k))

    def compute_dihedrals(self, coordinates):
        """
        Compute the potential energy and forces for the dihedrals and impropers.

        :param coordinates:
            The coordinates of the atoms.
        """

        for dihedral_key in ['dihedrals', 'impropers']:

            for (i, j, k, l), dih in self.force_field[dihedral_key].items():

                assert_msg_critical(
                    dih['type'] in ['Fourier', 'RB'],
                    'MMDriver: Invalid type for dihedral potential')

                if dih['type'] == 'RB':
                    RB_coefs = dih['RB_coefficients']
                elif dih['type'] == 'Fourier':
                    RB_coefs = self.get_RB_coefficients(dih['barrier'],
                                                        dih['phase'],
                                                        dih['periodicity'])

                (potential_energy, grad_i, grad_j, grad_k,
                 grad_l) = MMDriver.compute_Ryckaert_Bellemans(
                     coordinates[i], coordinates[j], coordinates[k],
                     coordinates[l], RB_coefs)

                if dihedral_key == 'dihedrals':
                    self.torsion_potential += potential_energy
                elif dihedral_key == 'impropers':
                    self.improper_potential += potential_energy

                self.gradient[i] += grad_i
                self.gradient[j] += grad_j
                self.gradient[k] += grad_k
                self.gradient[l] += grad_l

                self.exclusions.append((i, l))

    def compute_nonbonded(self, coordinates):
        """
        Compute the potential energy and forces for the non-bonded interactions.

        :param coordinates:
            The coordinates of the atoms.
        """

        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):

                if (i, j) in self.force_field['pairs']:
                    nonbonded_type = 'pair'
                    prefac_LJ = self.force_field['fudgeLJ']
                    prefac_QQ = self.force_field['fudgeQQ']
                elif (i, j) not in self.exclusions:
                    nonbonded_type = 'full'
                    prefac_LJ = 1.0
                    prefac_QQ = 1.0
                else:
                    nonbonded_type = None
                    prefac_LJ = None
                    prefac_QQ = None

                if nonbonded_type in ['pair', 'full']:

                    charge_i = self.force_field['atoms'][i]['charge']
                    charge_j = self.force_field['atoms'][j]['charge']

                    sigma_i = self.force_field['atoms'][i]['sigma']
                    sigma_j = self.force_field['atoms'][j]['sigma']

                    epsilon_i = self.force_field['atoms'][i]['epsilon']
                    epsilon_j = self.force_field['atoms'][j]['epsilon']

                    r_ij = coordinates[j] - coordinates[i]
                    distance_ij = np.linalg.norm(r_ij)
                    n_ij = r_ij / distance_ij

                    # Coulomb

                    # Coulomb's constant in MD units: 138.9354576...
                    ke = (4.184 * hartree_in_kcalpermol() * 0.1 *
                          bohr_in_angstrom())

                    self.coulomb_potential += prefac_QQ * ke * (
                        charge_i * charge_j) / distance_ij

                    g = -prefac_QQ * ke * (charge_i * charge_j) / distance_ij**2

                    self.gradient[i] += -g * n_ij
                    self.gradient[j] += g * n_ij

                    # VDW

                    epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                    sigma_ij = 0.5 * (sigma_i + sigma_j)

                    sigma_r_6 = (sigma_ij / distance_ij)**6
                    sigma_r_12 = sigma_r_6**2

                    self.vdw_potential += prefac_LJ * 4.0 * epsilon_ij * (
                        sigma_r_12 - sigma_r_6)

                    g = -prefac_LJ * 24.0 * epsilon_ij * (
                        2.0 * sigma_r_12 / distance_ij -
                        sigma_r_6 / distance_ij)

                    self.gradient[i] += -g * n_ij
                    self.gradient[j] += g * n_ij

    @staticmethod
    def get_RB_coefficients(barrier, phase, periodicity):
        """
        Get Ryckaert-Bellemans coefficients.

        :param barrier:
            The barrier of dihedral potential.
        :param phase:
            The phase of dihedral potential.
        :param periodicity:
            The periodicity of dihedral potential.

        :return:
            The Ryckaert-Bellemans coefficients.
        """

        assert_msg_critical(
            abs(phase) < 1.0e-6 or abs(phase - 180.0) < 1.0e-6,
            'MMDriver: Invalid phase for dihedral potential')

        assert_msg_critical(
            abs(periodicity) in [1, 2, 3, 4],
            'MMDriver: Invalid periodicity for dihedral potential')

        F_coefs = {x: 0.0 for x in [1, 2, 3, 4]}
        E_shift = 0.0

        if abs(phase) < 1.0e-6:
            # phase == 0 degree
            F_coefs[abs(periodicity)] += barrier
        else:
            # phase == 180 degree
            F_coefs[abs(periodicity)] -= barrier
            E_shift += 2.0 * barrier

        # JPCA 2021, 125, 2673-2681
        # Note that we also take into account the phases in Fourier series
        c0 = (F_coefs[1] + F_coefs[3] + 2.0 * F_coefs[4] + E_shift)
        c1 = -1.0 * F_coefs[1] + 3.0 * F_coefs[3]
        c2 = 2.0 * F_coefs[2] - 8.0 * F_coefs[4]
        c3 = -4.0 * F_coefs[3]
        c4 = 8.0 * F_coefs[4]
        c5 = 0.0

        return (c0, c1, c2, c3, c4, c5)

    def compute(self, molecule):
        """
        Compute the potential energy and gradient for the given molecule.

        :param molecule:
            The molecule.
        """

        # Convert coordinates to nm
        bohr_in_nm = bohr_in_angstrom() * 0.1

        # Get coordinates from the molecule in a.u.
        coordinates = molecule.get_coordinates_in_bohr() * bohr_in_nm

        # Initialization
        self.gradient = np.zeros((len(coordinates), 3))
        self.exclusions = []

        # For reporting components puposes only
        # TODO: use dictionary for energy components
        self.stretch_potential = 0.0
        self.bend_potential = 0.0
        self.torsion_potential = 0.0
        self.improper_potential = 0.0
        self.coulomb_potential = 0.0
        self.vdw_potential = 0.0
        self.energy = 0.0

        # Compute potential energy components and update atom forces
        self.compute_bonds(coordinates)
        self.compute_angles(coordinates)
        self.compute_dihedrals(coordinates)
        self.compute_nonbonded(coordinates)

        # energy
        self.energy = (self.stretch_potential + self.bend_potential +
                       self.torsion_potential + self.improper_potential +
                       self.coulomb_potential + self.vdw_potential)
        self.energy /= (4.184 * hartree_in_kcalpermol())

        # gradient
        self.gradient /= (4.184 * hartree_in_kcalpermol() * 10.0 /
                          bohr_in_angstrom())

    def get_energy(self):
        """
        Get the potential energy.

        :return:
            The potential energy.
        """

        return self.energy

    def get_gradient(self):
        """
        Get the gradient.

        :return:
            The gradient.
        """

        return self.gradient

    @staticmethod
    def compute_Ryckaert_Bellemans(coord_i, coord_j, coord_k, coord_l,
                                   RB_coefs):
        """
        Compute the dihedral potential and forces using the Ryckaert-Bellemans
        potential.

        :param coord_i:
            The coordinates of atom i.
        :param coord_j:
            The coordinates of atom j.
        :param coord_k:
            The coordinates of atom k.
        :param coord_l:
            The coordinates of atom l.
        :param RB_coefs:
            The Ryckaert-Bellemans coefficients.
        """

        # J. Comput. Chem. 2000, 21, 553-561

        r21 = coord_i - coord_j
        r32 = coord_j - coord_k
        r43 = coord_k - coord_l

        r21_d = np.linalg.norm(r21)
        r32_d = np.linalg.norm(r32)
        r43_d = np.linalg.norm(r43)

        r21_unit = r21 / r21_d
        r32_unit = r32 / r32_d
        r43_unit = r43 / r43_d

        cos123 = -np.dot(r21_unit, r32_unit)
        sin123 = np.sqrt(1.0 - cos123**2)

        cos234 = -np.dot(r32_unit, r43_unit)
        sin234 = np.sqrt(1.0 - cos234**2)

        cos134 = -np.dot(r21_unit, r43_unit)

        cos_phi = (cos123 * cos234 + cos134) / (sin123 * sin234)

        c0, c1, c2, c3, c4, c5 = RB_coefs

        potential_energy = (c0 - c1 * cos_phi + c2 * cos_phi**2 -
                            c3 * cos_phi**3 + c4 * cos_phi**4 - c5 * cos_phi**5)

        g = (-c1 + 2.0 * c2 * cos_phi - 3.0 * c3 * cos_phi**2 +
             4.0 * c4 * cos_phi**3 - 5.0 * c5 * cos_phi**4)

        r_21_32 = np.cross(r21_unit, r32_unit)
        r_43_32 = np.cross(r43_unit, r32_unit)

        const_1 = np.dot(r43_unit, r_21_32) / (r21_d * sin123**3 * sin234)
        const_4 = np.dot(r43_unit, r_21_32) / (r43_d * sin234**3 * sin123)

        c123 = r21_d * cos123 / r32_d - 1.0
        b432 = r43_d * cos234 / r32_d

        grad_i = -const_1 * r_21_32 * g
        grad_l = -const_4 * r_43_32 * g

        grad_j = c123 * grad_i - b432 * grad_l
        grad_k = -(grad_i + grad_j + grad_l)

        return (potential_energy, grad_i, grad_j, grad_k, grad_l)

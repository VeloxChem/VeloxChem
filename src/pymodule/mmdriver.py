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
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol
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

        # TODO: separate the basic part of MMDriver into another class to avoid
        # cyclic import

        # NOTE: Never use this method inside MMForceFieldGenerator.
        # Otherwise it will be a cyclic import.

        from .mmforcefieldgenerator import MMForceFieldGenerator

        ff_generator = MMForceFieldGenerator()
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

                # Assert msg critical if the dihedral potential is not Fourier
                assert_msg_critical(
                    dih['type'] == 'Fourier',
                    'MMDriver: Only Fourier dihedral potential is supported')

                # Multiple dihedral potentials can be true or false
                # Use get method with default value False
                if dih.get('multiple', False):
                    # Barrier, phase, and periodicity are in lists
                    for barrier, phase, periodicity in zip(
                            dih['barrier'], dih['phase'], dih['periodicity']):

                        (potential_energy, grad_i, grad_j, grad_k,
                         grad_l) = self.compute_Ryckaert_Bellemans(
                             coordinates[i], coordinates[j], coordinates[k],
                             coordinates[l], barrier, phase, periodicity)

                        if dihedral_key == 'dihedrals':
                            self.torsion_potential += potential_energy
                        elif dihedral_key == 'impropers':
                            self.improper_potential += potential_energy

                        self.gradient[i] += grad_i
                        self.gradient[j] += grad_j
                        self.gradient[k] += grad_k
                        self.gradient[l] += grad_l

                else:
                    (potential_energy, grad_i, grad_j, grad_k,
                     grad_l) = self.compute_Ryckaert_Bellemans(
                         coordinates[i], coordinates[j], coordinates[k],
                         coordinates[l], dih['barrier'], dih['phase'],
                         dih['periodicity'])

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
                    ke = (hartree_in_kjpermol() * 0.1 * bohr_in_angstrom())

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
            abs(periodicity) in [1, 2, 3, 4, 5, 6],
            'MMDriver: Invalid periodicity for dihedral potential')

        F_coefs = {x: 0.0 for x in [1, 2, 3, 4, 5, 6]}
        E_shift = 0.0

        if abs(phase) < 1.0e-6:
            # phase == 0 degree
            F_coefs[abs(periodicity)] += barrier
        else:
            # phase == 180 degree
            F_coefs[abs(periodicity)] -= barrier
            E_shift += 2.0 * barrier

        # See e.g. JPCA 2021, 125, 2673-2681 for periodicity 1-4.
        # F_coefs[5] and F_coefs[6] were added for periodicity up to 6,
        # based on Chebyshev polynomials of the first kind
        c0 = (F_coefs[1] + F_coefs[3] + 2.0 * F_coefs[4] + F_coefs[5] + E_shift)
        c1 = -1.0 * F_coefs[1] + 3.0 * F_coefs[3] - 5.0 * F_coefs[5]
        c2 = 2.0 * F_coefs[2] - 8.0 * F_coefs[4] + 18.0 * F_coefs[6]
        c3 = -4.0 * F_coefs[3] + 20.0 * F_coefs[5]
        c4 = 8.0 * F_coefs[4] - 48.0 * F_coefs[6]
        c5 = -16.0 * F_coefs[5]
        c6 = 32.0 * F_coefs[6]

        return (c0, c1, c2, c3, c4, c5, c6)

    def compute(self, molecule):
        """
        Compute the potential energy and gradient for the given molecule.

        :param molecule:
            The molecule.
        """

        # TODO: automatically generate and load force field if self.force_field
        # is empty

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
        self.energy /= hartree_in_kjpermol()

        # gradient
        self.gradient /= (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())

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

    def compute_Ryckaert_Bellemans(self, coord_i, coord_j, coord_k, coord_l,
                                   barrier, phase, periodicity):
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
        :param barrier:
            The barrier.
        :param phase:
            The phase.
        :param periodicity:
            The periodicity.

        :return:
            The dihedral potential energy and gradients.
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

        RB_coefs = self.get_RB_coefficients(barrier, phase, periodicity)

        c0, c1, c2, c3, c4, c5, c6 = RB_coefs

        potential_energy = (c0 - c1 * cos_phi + c2 * cos_phi**2 -
                            c3 * cos_phi**3 + c4 * cos_phi**4 -
                            c5 * cos_phi**5 + c6 * cos_phi**6)

        g = (-c1 + 2.0 * c2 * cos_phi - 3.0 * c3 * cos_phi**2 +
             4.0 * c4 * cos_phi**3 - 5.0 * c5 * cos_phi**4 +
             6.0 * c6 * cos_phi**5)

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

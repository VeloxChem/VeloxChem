#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

import h5py
import numpy as np
import py3Dmol as p3d
from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import Molecule
from .veloxchemlib import MolecularBasis
from .errorhandler import assert_msg_critical
from .visualizationdriver import VisualizationDriver


class ExcitedStateAnalysisDriver:
    """
    Implements excited state descriptors at TDDFT level of theory.
    """

    def __init__(self):
        """
        Initialize excited state analysis driver.
        """
        self.fragment_dict = None

        self.include_transition_density_matrix = False
        self.include_particle_hole_density_matrices = False
        self.include_ct_matrix = False
        self.include_participation_ratios = False
        self.include_average_particle_hole_position = False

    def update_settings(self, fragment_dict=None):
        """
        Update fragment dictionary and descriptor settings.

        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.
        """

        if fragment_dict is not None:
            self.fragment_dict = fragment_dict

    def compute(self, state_index, molecule, basis, scf_tensors, rsp_tensors):
        """
        Computes dictionary containing excited state descriptors for
        excited state nstate.

        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.
        :param nstate:
            The excited state for which the descriptors are computed.
        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.
        :param scf_checkpoint_file:
            The checkpoint file from the converged SCF calculation.
        :param rsp_checkpoint_file:
            The checkpoint file for the rsp calculation.

        :return:
            A dictionary containing the transition, hole and particle density matrices
            in MO and AO basis, the charge transfer matrix,
            the hole/particle/average participation ratios,
            the average hole/particle position and their difference vector,
            and the relative charge transfer length.
        """

        assert_msg_critical(self.fragment_dict is not None,
                            'fragment_dict not defined ')

        num_exc_states = sum([1 for key in rsp_tensors if key.startswith('S')])

        assert_msg_critical(state_index <= num_exc_states,
                            'excited state not included in rsp calculation')
        ret_dict = {}

        # necessary objects for some descriptor calculations
        trans_dens_mo, trans_dens_ao = self.compute_transition_density_matrix(
            molecule, scf_tensors, rsp_tensors, state_index)
        ct_matrix = self.compute_ct_matrix(molecule, basis, scf_tensors,
                                           trans_dens_ao, self.fragment_dict)

        if self.include_transition_density_matrix is True:
            ret_dict['Transition_density_matrix_MO'] = trans_dens_mo
            ret_dict['Transition_density_matrix_AO'] = trans_dens_ao

        if self.include_particle_hole_density_matrices is True:
            ret_dict['Hole_density_matrix_MO'], ret_dict[
                'hole_density_matrix_AO'] = self.compute_hole_density_matrix(
                    molecule, scf_tensors, rsp_tensors, state_index)
            ret_dict['Particle_density_matrix_MO'], ret_dict[
                'Particle_density_matrix_AO'] = self.compute_particle_density_matrix(
                    molecule, scf_tensors, rsp_tensors, state_index)

        if self.include_ct_matrix is True:
            ret_dict['CT_matrix'] = ct_matrix

        if self.include_participation_ratios is True:
            ret_dict['Hole_participation_ratio'], ret_dict[
                'Particle_participation_ratio'], ret_dict[
                    'Average_participation_ratio'] = self.compute_participation_ratios(
                        ct_matrix)

        if self.include_average_particle_hole_position is True:
            ret_dict['Average_particle_position'], ret_dict[
                'Average_hole_position'], ret_dict[
                    'Average_difference_vector'] = self.compute_avg_position(
                        molecule, ct_matrix, self.fragment_dict)
            ret_dict['Relative_CT_length'] = self.compute_relative_ct_length(
                molecule, ret_dict['Average_difference_vector'])

        return ret_dict

    def read_scf_checkpoint_file(self, fname):
        """
        Reads data from hdf5 scf checkpoint file
        and returns it as a dictionary.

        :param filename:
            The name of the scf checkpoint file.

        :return:
            The scf checkpoint file dictionary.
        """

        dict = {}
        h5f = h5py.File(fname, "r")

        for key in h5f.keys():
            data = np.array(h5f.get(key))
            dict[key] = data
        h5f.close()

        assert_msg_critical(
            'atom_coordinates' in dict,
            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: ' +
            'atom_coordinates not found')

        assert_msg_critical(
            'nuclear_charges' in dict,
            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: ' +
            'nuclear_charges not found')

        assert_msg_critical(
            'basis_set' in dict,
            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: ' +
            'basis_set not found')

        assert_msg_critical(
            'C_alpha' in dict,
            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: ' +
            'C_alpha not found')

        assert_msg_critical(
            'S' in dict,
            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: ' +
            'S not found')

        return dict

    def read_rsp_checkpoint_file(self, fname):
        """
        Read data from hdf5 rsp checkpoint file
        and return it as a dictionary.

        :param filename:
            The name of the rsp checkpoint file.

        :return:
            The rsp checkpoint file dictionary.
        """

        dict = {}
        h5f = h5py.File(fname, "r")

        for key in h5f.keys():
            data = np.array(h5f.get(key))
            dict[key] = data
        h5f.close()

        return dict

    def create_molecule_basis(self, scf_tensors):
        """
        Creates molecule and basis objects from scf dictionary.

        :param scf_tensors:
            The dictionary containing the scf tensors.

        :return:
            A tuple containing the Molecule and Basis objects.
        """

        coordinates = scf_tensors['atom_coordinates']
        nuclear_charges = np.array(scf_tensors['nuclear_charges'])
        nuclear_charges = nuclear_charges.astype(int)
        basis_set_label = scf_tensors['basis_set'][0].decode("utf-8")
        molecule = Molecule(nuclear_charges, coordinates, units="au")
        basis = MolecularBasis.read(molecule, basis_set_label)

        return molecule, basis

    def compute_transition_density_matrix(self, molecule, scf_tensors,
                                          rsp_tensors, state_index):
        """
        Computes the transition density matrix (TDM) for state
        nstate in MO and AO basis.

        :param molecule:
            The Molecule object
        :param scf_tensors:
            The dictionary containing the scf tensors.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.
        :param nstate:
            The excited state for which the TDM is computed.

        :return:
            A tuple containing the TDM in MO and AO basis.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = scf_tensors["C_alpha"].shape[1]
        nvirt = norb - nocc

        excstate = "S" + str(state_index)
        x_f = rsp_tensors[excstate]

        nexc = nocc * nvirt
        z_f = x_f[:nexc]
        if x_f.shape[0] == nexc:
            trans_dens_mo = np.reshape(z_f, (nocc, nvirt))
        else:
            y_f = x_f[nexc:]
            trans_dens_mo = np.reshape(z_f - y_f, (nocc, nvirt))

        trans_dens_ao = np.linalg.multi_dot([
            scf_tensors["C_alpha"][:, :nocc], trans_dens_mo,
            scf_tensors["C_alpha"][:, nocc:].T
        ])

        return trans_dens_mo, trans_dens_ao

    def compute_hole_density_matrix(self, molecule, scf_tensors, rsp_tensors,
                                    state_index):
        """
        Computes the hole density matrix for state nstate in MO and AO basis.

        :param molecule:
            The Molecule object
        :param scf_tensors:
            The dictionary containing the scf tensors.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.
        :param nstate:
            The excited state for which the hole density matrix is computed.

        :return:
            A tuple containing the hole density matrix in MO and AO basis.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = scf_tensors["C_alpha"].shape[1]
        nvirt = norb - nocc

        excstate = "S" + str(state_index)
        x_f = rsp_tensors[excstate]

        nexc = nocc * nvirt
        z_f = x_f[:nexc]
        if x_f.shape[0] == nexc:
            trans_dens_mo = np.reshape(z_f, (nocc, nvirt))
        else:
            y_f = x_f[nexc:]
            trans_dens_mo = np.reshape(z_f - y_f, (nocc, nvirt))

        hole_dens_mo = np.matmul(trans_dens_mo, trans_dens_mo.T)
        hole_dens_ao = np.linalg.multi_dot([
            scf_tensors["C_alpha"][:, :nocc], hole_dens_mo,
            scf_tensors["C_alpha"][:, :nocc].T
        ])

        return hole_dens_mo, hole_dens_ao

    def compute_particle_density_matrix(self, molecule, scf_tensors,
                                        rsp_tensors, state_index):
        """
        Computes the particle density matrix for state nstate in MO and
        AO basis.

        :param molecule:
            The Molecule object
        :param scf_tensors:
            The dictionary containing the scf tensors.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.
        :param nstate:
            The excited state for which the particle density matrix is
            computed.

        :return:
            A tuple containing the particle density matrix in MO and AO basis.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = scf_tensors["C_alpha"].shape[1]
        nvirt = norb - nocc

        excstate = "S" + str(state_index)
        x_f = rsp_tensors[excstate]

        nexc = nocc * nvirt
        z_f = x_f[:nexc]
        if x_f.shape[0] == nexc:
            trans_dens_mo = np.reshape(z_f, (nocc, nvirt))
        else:
            y_f = x_f[nexc:]
            trans_dens_mo = np.reshape(z_f - y_f, (nocc, nvirt))

        part_dens_mo = np.matmul(trans_dens_mo.T, trans_dens_mo)
        part_dens_ao = np.linalg.multi_dot([
            scf_tensors["C_alpha"][:, nocc:], part_dens_mo,
            scf_tensors["C_alpha"][:, nocc:].T
        ])

        return part_dens_mo, part_dens_ao

    def map_fragments_to_atomic_orbitals(self, molecule, basis, fragment_dict):
        """
        Maps the atom indices in the fragment dictionary to indices of atomic
        orbitals centered on the corresponding atom.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param fragment_dict:
            The dictionary containing indices of atoms in each fragment.

        :return:
            The dictionary containing the indices of the atomic orbitals
            centered in each fragment.
        """

        vis_drv = VisualizationDriver()
        ao_list = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)

        ao_map_fragments = {}
        for fragment in fragment_dict:
            ao_map_fragments[fragment] = []
            for i in fragment_dict[fragment]:
                ao_map_fragments[fragment] += ao_list[i -
                                                      1]  # i is 1-based index
        return ao_map_fragments

    def compute_ct_matrix(self, molecule, basis, scf_tensors, T, fragment_dict):
        """
        Computes the charge transfer (CT) matrix from the
        transition density matrix (TDM).

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param T:
            The TDM in AO basis.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            The CT matrix.
        """

        ao_map_fragments = self.map_fragments_to_atomic_orbitals(
            molecule, basis, fragment_dict)
        S = scf_tensors["S"]

        TS = np.matmul(T, S)
        ST = np.matmul(S, T)
        STS = np.linalg.multi_dot([S, T, S])
        TS_ST = TS * ST  # element by element
        T_STS = T * STS

        n_fragments = len(fragment_dict.keys())
        ct_matrix = np.zeros((n_fragments, n_fragments))

        for a, A in enumerate(fragment_dict.keys()):
            for b, B in enumerate(fragment_dict.keys()):
                TSST_AB = TS_ST[ao_map_fragments[A]][:, ao_map_fragments[B]]
                T_STS_AB = T_STS[ao_map_fragments[A]][:, ao_map_fragments[B]]
                ct_matrix[a, b] = np.einsum("mn->", TSST_AB + T_STS_AB)

        return (1 / 2) * ct_matrix

    def compute_participation_ratios(self, ct_matrix):
        """
        Computes the particle, hole, and average participation ratios
        from the charge transfer (CT) matrix.

        :param ct_matrix:
            The CT matrix.

        :return:
            A tuple containing the hole, particle, and average participation ratios.
        """

        omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / omega
        particle_weight = np.sum(ct_matrix, axis=0) / omega
        hole_participation_ratio = 1 / np.sum(hole_weight**2)
        particle_participation_ratio = 1 / np.sum(particle_weight**2)
        avg_participation_ratio = (hole_participation_ratio +
                                   particle_participation_ratio) / 2

        return hole_participation_ratio, particle_participation_ratio, avg_participation_ratio

    def compute_fragment_coordinates(self, molecule, fragment_dict):
        """
        Computes average coordinates of fragments in fragment_dict.

        :param molecule:
            The molecule.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            A list containing the average coordinates of each fragment.
        """

        coordinates = molecule.get_coordinates_in_bohr()
        avg_coord_list = []
        for i, A in enumerate(fragment_dict):
            coord_array = coordinates[np.array(fragment_dict[A]) - 1, :]
            avg_coord_list.append(
                np.sum(coord_array, axis=0) / coord_array.shape[0])
        return avg_coord_list

    def max_dist_from_mol_center(self, molecule):
        """
        Computes the maximum distance from an atom in the molecule to the
        center of the molecule.

        :param molecule:
            The molecule.

        :return:
            The maximum distance from an atom in the molecule to the
            center of the molecule.
        """

        coordinates = molecule.get_coordinates_in_angstrom()
        center = np.sum(coordinates, axis=0) / coordinates.shape[0]
        print(center)
        atom_center_distance = []
        for i in range(coordinates.shape[0]):
            atom_center_distance.append(np.linalg.norm(coordinates[i] - center))
        return np.max(np.array(atom_center_distance))

    def compute_avg_position(self, molecule, ct_matrix, fragment_dict):
        """
        Computes the average particle position, average hole position, and the
        the difference vector between them.

        :param molecule:
            The molecule.
        :param ct_matrix:
            The charge transfer matrix.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            A tuple containing the average particle position,
            average hole position, and the difference vector
            from average hole to average particle.
        """

        avg_coord_list = self.compute_fragment_coordinates(
            molecule, fragment_dict)

        omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / omega
        particle_weight = np.sum(ct_matrix, axis=0) / omega

        avg_hole_position = [0, 0, 0]
        avg_particle_position = [0, 0, 0]

        for j in range(len(avg_coord_list)):
            avg_hole_position += (avg_coord_list[j] * hole_weight[j] *
                                  bohr_in_angstrom())
            avg_particle_position += (avg_coord_list[j] * particle_weight[j] *
                                      bohr_in_angstrom())
        avg_diff_vec = avg_particle_position - avg_hole_position
        return avg_particle_position, avg_hole_position, avg_diff_vec

    def compute_relative_ct_length(self, molecule, avg_diff_vec):
        """
        Compute the relative charge transfer (CT) length.

        :param molecule:
            The molecule.
        :param ET_net:
            The difference vector from average position of hole to average
            position of particle.

        :return:
            The relative charge transfer length.
        """
        return np.linalg.norm(avg_diff_vec) / self.max_dist_from_mol_center(
            molecule)

    def show_avg_ph_position(self,
                             molecule,
                             avg_hole_position,
                             avg_particle_position,
                             width=400,
                             height=300):
        """
        Displays the molecule and the average positions of the particle and hole.

        :param molecule:
            The molecule.
        :param avg_hole_pos:
            The average position of the hole.
        :param avg_particle_pos:
            The average positions of the particle.
        :param width:
            The width of the plot.
        :param height:
            The height of the plot.
        """
        r_H = avg_hole_position
        r_E = avg_particle_position
        viewer = p3d.view(width=width, height=height)
        viewer.addModel(molecule.get_xyz_string(), "xyz")
        viewer.setStyle({"line": {}, "sphere": {"scale": 0.05}})
        viewer.addArrow({
            'start': {
                'x': r_H[0],
                'y': r_H[1],
                'z': r_H[2]
            },
            'end': {
                'x': r_E[0],
                'y': r_E[1],
                'z': r_E[2]
            },
            'radius': 0.06,
            'color': 'green',
        })
        viewer.addSphere({
            'center': {
                'x': r_H[0],
                'y': r_H[1],
                'z': r_H[2]
            },
            'radius': 0.12,
            'color': 'blue',
        })
        viewer.addSphere({
            'center': {
                'x': r_E[0],
                'y': r_E[1],
                'z': r_E[2]
            },
            'radius': 0.12,
            'color': 'red',
        })
        viewer.show()

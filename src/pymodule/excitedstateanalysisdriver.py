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
import sys
import time as tm

from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import Molecule
from .veloxchemlib import MolecularBasis
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .visualizationdriver import VisualizationDriver
from .lreigensolver import LinearResponseEigenSolver


class ExcitedStateAnalysisDriver:
    """
    Implements excited state descriptors at TDDFT level of theory.

    :param ostream:
        The output stream.

    Instance variables
        - fragment_dict: dictionary of indices for atoms in each fragment.
        - include_particle_hole_density_matrices: The flag for particle and
          hole density matrices.
        - include_participation_ratios: The flag for participation ratios.
        - include_avg_particle_hole_position: The flag for average particle.
          and hole positions.
        - include_relative_ct_length: The flag for relative CT length.
        - _tda: The flag for the Tamm-Dankoff Approximation.
    """

    def __init__(self, ostream=None):
        """
        Initialize excited state analysis driver to default setup.
        """

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.fragment_dict = None

        self.include_particle_hole_density_matrices = True
        self.include_participation_ratios = True
        self.include_avg_particle_hole_position = True
        self.include_relative_ct_length = True
        self._tda = False

        # output stream
        self.ostream = ostream

    def compute(self, state_index, molecule, basis, scf_tensors, rsp_tensors):
        """
        Computes dictionary containing excited state descriptors for
        excited state nstate.

        :param state_index:
            The excited state for which the descriptors are computed.
        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.

        :return:
            A dictionary containing the transition, hole and particle density matrices
            in MO and AO basis, the charge transfer matrix,
            the hole/particle/average participation ratios,
            the average hole/particle position and their difference vector,
            and the relative charge transfer length.
        """

        assert_msg_critical(self.fragment_dict is not None,
                            'fragment_dict not defined ')

        num_exc_states = len(rsp_tensors['eigenvalues'])

        assert_msg_critical(state_index <= num_exc_states,
                            'excited state not included in rsp calculation')

        ret_dict = {}

        self.print_header('Excited State Analysis.', state_index)

        start_time = tm.time()

        # necessary objects for some descriptor calculations
        tdens_mo, tdens_ao = self.compute_transition_density_matrix(
            molecule, scf_tensors, rsp_tensors, state_index)

        self.ostream.print_info("Transition density matrix computed in " +
                                "{:.3f} sec.".format(tm.time() - start_time))
        self.ostream.print_blank()

        tmp_start_time = tm.time()
        ct_matrix = self.compute_ct_matrix(molecule, basis, scf_tensors,
                                           tdens_ao, self.fragment_dict)
        self.ostream.print_info("Charge transfer matrix computed in " +
                                "{:.3f} sec.".format(tm.time() -
                                                     tmp_start_time))
        self.ostream.print_blank()

        tmp_start_time = tm.time()
        ret_dict['transition_density_matrix_MO'] = tdens_mo
        ret_dict['transition_density_matrix_AO'] = tdens_ao

        ret_dict['ct_matrix'] = ct_matrix

        if self.include_particle_hole_density_matrices:

            ret_dict['hole_density_matrix_MO'], ret_dict[
                'hole_density_matrix_AO'] = self.compute_hole_density_matrix(
                    molecule, scf_tensors, rsp_tensors, state_index)
            ret_dict['particle_density_matrix_MO'], ret_dict[
                'particle_density_matrix_AO'] = self.compute_particle_density_matrix(
                    molecule, scf_tensors, rsp_tensors, state_index)
            self.ostream.print_info(
                "Particle and hole density matrices computed in " +
                "{:.3f} sec.".format(tm.time() - tmp_start_time))
            self.ostream.print_blank()
            tmp_start_time = tm.time()

        if self.include_participation_ratios:

            ret_dict['hole_participation_ratio'], ret_dict[
                'particle_participation_ratio'], ret_dict[
                    'avg_participation_ratio'] = self.compute_participation_ratios(
                        ct_matrix)
            self.ostream.print_info("Participation ratios computed in " +
                                    "{:.3f} sec.".format(tm.time() -
                                                         tmp_start_time))
            self.ostream.print_blank()
            tmp_start_time = tm.time()

        if self.include_avg_particle_hole_position or self.include_relative_ct_length:

            ret_dict['avg_particle_position'], ret_dict[
                'avg_hole_position'], ret_dict[
                    'avg_difference_vector'] = self.compute_avg_position(
                        molecule, ct_matrix, self.fragment_dict)
            ret_dict['relative_ct_length'] = self.compute_relative_ct_length(
                molecule, ret_dict['avg_difference_vector'])
            self.ostream.print_info(
                "Average positions and relative ct length computed in " +
                "{:.3f} sec.".format(tm.time() - tmp_start_time))
            self.ostream.print_blank()

        self.ostream.print_info(
            "Total time spent performing excited state analysis: " +
            "{:.3f} sec.".format(tm.time() - start_time))

        self.ostream.print_blank()
        self.ostream.print_blank()

        self.print_results(ret_dict)

        self.ostream.flush()

        return ret_dict
    
    def read_from_h5(self, filename):
        """
        Reads data from hdf5 checkpoint file and 
        returns it as tuple containing the scf and rsp dictionaries.

        :param filename:
            The name of the checkpoint file.

        :return:
            A tuple containing the scf and rsp dictionaries.
        """
        h5f = h5py.File(filename, "r")

        scf_tensors = {}
        rsp_tensors = {}
        for key in h5f:
            if key == "atom_coordinates":
                scf_tensors[key] = np.array(h5f[key])
            if key == "nuclear_charges":
                scf_tensors[key] = np.array(h5f[key])
            if key == "basis_set":
                scf_tensors[key] = np.array(h5f[key])
            if key == "scf":
                scf_tensors_dict = dict(h5f.get(key))
                for scf_key in scf_tensors_dict:
                    scf_tensors[scf_key] = np.array(scf_tensors_dict[scf_key])
            elif key == "rsp":
                rsp_tensors_dict = dict(h5f.get(key))
                for rsp_key in rsp_tensors_dict:
                    rsp_tensors[rsp_key] = np.array(rsp_tensors_dict[rsp_key]) 
        h5f.close()

        return scf_tensors, rsp_tensors
    
    def format_rsp_tensors(self, rsp_tensors):
        """
        Reads the eigenvectors in rsp_tensors and adds them
        as numpy arrays with individual keys to the dictionary.

        :param rsp_tensors:
            The dictionary containing the results from a 
            rsp calculation.

        :return:
            The rsp results dictionary in the format needed for
            the excited state analysis driver class.
        """

        if 'eigenvectors_distributed' in rsp_tensors.keys():
            rsp_drv = LinearResponseEigenSolver()
            num_states = len(rsp_tensors['eigenvectors_distributed'])
            for i in range(num_states):
                name = 'S'+str(i+1)
                rsp_tensors[name] = rsp_drv.get_full_solution_vector(rsp_tensors['eigenvectors_distributed'][i])
        else:
            num_states = len(rsp_tensors['eigenvectors'][0, :])
            for i in range(num_states):
                name = 'S'+str(i+1)
                rsp_tensors[name] = rsp_tensors['eigenvectors'][:, i]
        rsp_tensors['formatted'] = True

        return rsp_tensors

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
        :param state_index:
            The excited state for which the TDM is computed.

        :return:
            A tuple containing the TDM in MO and AO basis.
        """

        if any(key.startswith('eigenvector') for key in rsp_tensors.keys()) and 'formatted' not in rsp_tensors.keys():
            rsp_tensors = self.format_rsp_tensors(rsp_tensors)

        nocc = molecule.number_of_alpha_electrons()
        norb = scf_tensors["C_alpha"].shape[1]
        nvirt = norb - nocc

        
        excstate = "S" + str(state_index)
        eigvec = rsp_tensors[excstate]

        nexc = nocc * nvirt
        z_mat = eigvec[:nexc]
        if eigvec.shape[0] == nexc:
            tdens_mo = np.reshape(z_mat, (nocc, nvirt))
        else:
            y_mat = eigvec[nexc:]
            tdens_mo = np.reshape(z_mat - y_mat, (nocc, nvirt))

        tdens_ao = np.linalg.multi_dot([
            scf_tensors["C_alpha"][:, :nocc], tdens_mo,
            scf_tensors["C_alpha"][:, nocc:].T
        ])

        return tdens_mo, tdens_ao

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
        :param state_index:
            The excited state for which the hole density matrix is computed.

        :return:
            A tuple containing the hole density matrix in MO and AO basis.
        """

        if any(key.startswith('eigenvector') for key in rsp_tensors.keys()) and 'formatted' not in rsp_tensors.keys():
            rsp_tensors = self.format_rsp_tensors(rsp_tensors)

        nocc = molecule.number_of_alpha_electrons()

        tdens_mo, tdens_ao = self.compute_transition_density_matrix(
            molecule, scf_tensors, rsp_tensors, state_index)

        hole_dens_mo = -np.matmul(tdens_mo, tdens_mo.T)
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
        :param state_index:
            The excited state for which the particle density matrix is
            computed.

        :return:
            A tuple containing the particle density matrix in MO and AO basis.
        """

        if any(key.startswith('eigenvector') for key in rsp_tensors.keys()) and 'formatted' not in rsp_tensors.keys():
            rsp_tensors = self.format_rsp_tensors(rsp_tensors)

        nocc = molecule.number_of_alpha_electrons()

        tdens_mo, tdens_ao = self.compute_transition_density_matrix(
            molecule, scf_tensors, rsp_tensors, state_index)

        part_dens_mo = np.matmul(tdens_mo.T, tdens_mo)
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
        hole_pr = 1 / np.sum(hole_weight**2)
        particle_pr = 1 / np.sum(particle_weight**2)
        avg_pr = (hole_pr + particle_pr) / 2

        return hole_pr, particle_pr, avg_pr

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
        :param avg_diff_vec:
            The difference vector from average position of hole to average
            position of particle.

        :return:
            The relative charge transfer length.
        """
        return np.linalg.norm(avg_diff_vec) / self.max_dist_from_mol_center(
            molecule)

    def show_avg_ph_position(self,
                             molecule,
                             descriptor_dict,
                             width=400,
                             height=300):
        """
        Displays the molecule and the average positions of the particle and hole.

        :param molecule:
            The molecule.
        :param descriptor_dict:
            The dictionary containing the excited state descriptors.
        :param width:
            The width of the plot.
        :param height:
            The height of the plot.
        """

        avg_hole_position = descriptor_dict['avg_hole_position']
        avg_particle_position = descriptor_dict['avg_particle_position']

        viewer = p3d.view(width=width, height=height)
        viewer.addModel(molecule.get_xyz_string(), "xyz")
        viewer.setStyle({"line": {}, "sphere": {"scale": 0.05}})
        viewer.addCylinder({
            'start': {
                'x': avg_hole_position[0],
                'y': avg_hole_position[1],
                'z': avg_hole_position[2]
            },
            'end': {
                'x': avg_particle_position[0],
                'y': avg_particle_position[1],
                'z': avg_particle_position[2]
            },
            'radius': 0.06,
            'color': 'darkcyan',
        })
        viewer.addSphere({
            'center': {
                'x': avg_hole_position[0],
                'y': avg_hole_position[1],
                'z': avg_hole_position[2]
            },
            'radius': 0.12,
            'color': 'red',
        })
        viewer.addSphere({
            'center': {
                'x': avg_particle_position[0],
                'y': avg_particle_position[1],
                'z': avg_particle_position[2]
            },
            'radius': 0.12,
            'color': 'blue',
        })
        viewer.zoomTo()
        viewer.show()

    def print_header(self, title, state_index):
        """
        Print header in output stream.

        :param title:
            The name of the driver.
        :param state_index:
            The excited state index
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s}'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 60

        cur_str = 'Excited state index             : ' + str(state_index)
        self.ostream.print_header(cur_str.ljust(str_width))

        for key in self.fragment_dict:
            if len(self.fragment_dict[key]
                  ) >= 7:  # Splice atomic index list if longer than 7 indices
                fragment_slice = [
                    self.fragment_dict[key][i:i + 7]
                    for i in range(0, len(self.fragment_dict[key]), 7)
                ]
                cur_str = f'Fragment {key:15s}        : ' + str(
                    fragment_slice[0])[1:-1]
                self.ostream.print_header(cur_str.ljust(str_width))
                for i in range(1, len(fragment_slice)):
                    cur_str = '                                  ' + str(
                        fragment_slice[i])[1:-1]
                    self.ostream.print_header(cur_str.ljust(str_width))
            else:
                cur_str = f'Fragment {key:15s}        : ' + str(
                    self.fragment_dict[key])[1:-1]
                self.ostream.print_header(cur_str.ljust(str_width))

        # is this part necessary and if so, what is the best way to add it?
        #if self._tda:
        #    cur_str = 'Tamm-Dancoff Approximation      : yes'
        #    self.ostream.print_header(cur_str.ljust(str_width))
        #else:
        #    cur_str = 'Tamm-Dancoff Approximation      : no'
        #    self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_avg_particle_hole_positions(self, descriptor_dict):
        """
        Prints average particle and hole positions.

        :param descriptor_dict:
            The dictionary of descriptors.
        """

        r_p = descriptor_dict['avg_particle_position']
        r_h = descriptor_dict['avg_hole_position']

        valstr = "Average Particle and Hole Positions (Angstrom)"
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = '                     '
        valstr += '{:>13s}{:>13s}{:>13s}'.format('X', 'Y', 'Z')
        self.ostream.print_header(valstr.ljust(92))

        valstr = 'Particle:            '
        valstr += '{:13.6f}{:13.6f}{:13.6f}'.format(r_p[0], r_p[1], r_p[2])
        self.ostream.print_header(valstr.ljust(92))

        valstr = 'Hole:                '
        valstr += '{:13.6f}{:13.6f}{:13.6f}'.format(r_h[0], r_h[1], r_h[2])
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def print_participation_ratios(self, descriptor_dict):
        """
        Prints participation ratios.

        :param descriptor_dict:
            The dictionary of descriptors.
        """
        pr = [None] * 3

        pr[0] = descriptor_dict['particle_participation_ratio']
        pr[1] = descriptor_dict['hole_participation_ratio']
        pr[2] = descriptor_dict['avg_participation_ratio']

        valstr = "Participation Ratios"
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = '                     '
        valstr += '{:>13s}{:>13s}{:>13s}'.format('Particle', 'Hole', 'Average')
        self.ostream.print_header(valstr.ljust(92))

        valstr = '                    '
        valstr += '{:13.6f}{:13.6f}{:13.6f}'.format(pr[0], pr[1], pr[2])
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def print_relative_charge_transfer_length(self, descriptor_dict):
        """
        Prints relative charge transfer length.

        :param descriptor_dict:
            The dictionary of descriptors.
        """
        rel_ct_length = descriptor_dict['relative_ct_length']

        valstr = "Relative Charge Transfer Length"
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = 'Distance:                 '
        valstr += '{:.6f}'.format(rel_ct_length)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def print_results(self, descriptor_dict):
        """
        Prints results of analysis.

        :param descriptor_dict:
            The dictionary of descriptors.
        """

        if self.include_avg_particle_hole_position or self.include_relative_ct_length:
            self.print_avg_particle_hole_positions(descriptor_dict)

        if self.include_participation_ratios:
            self.print_participation_ratios(descriptor_dict)

        if self.include_avg_particle_hole_position or self.include_relative_ct_length:
            self.print_relative_charge_transfer_length(descriptor_dict)

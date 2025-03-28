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
from .densityviewer import DensityViewer


class ExcitedStateAnalysisDriver:
    """
    Implements excited state descriptors at TDDFT level of theory.

    :param ostream:
        The output stream.

    Instance variables
        - fragment_dict: dictionary of indices for atoms in each fragment.
    """

    def __init__(self, ostream=None):
        """
        Initialize excited state analysis driver to default setup.
        """

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.fragment_dict = None

        # output stream
        self.ostream = ostream

    def compute(self, molecule, basis, scf_results, rsp_results, state_index=1):
        """
        Computes dictionary containing excited state descriptors for
        excited state nstate.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param state_index:
            The excited state for which the descriptors are
            computed (indexing starts at 1).

        :return:
            A dictionary containing the transition, hole and particle density matrices
            in MO and AO basis, the charge transfer matrix,
            the hole/particle/average participation ratios,
            the average hole/particle position and their difference vector,
            and the charge transfer length.
        """

        assert_msg_critical(self.fragment_dict is not None,
                            'fragment_dict not defined ')

        num_exc_states = len(rsp_results['eigenvalues'])

        assert_msg_critical(state_index <= num_exc_states,
                            'excited state not included in rsp calculation')

        ret_dict = {}

        self.print_header('Excited State Analysis.', state_index)

        start_time = tm.time()
        dens_dict = self.compute_density_matrices(molecule, scf_results,
                                                  rsp_results, state_index)
        tdens_ao = dens_dict['transition_density_matrix_AO']
        ret_dict.update(dens_dict)
        self.ostream.print_info("Density matrices computed in " +
                                "{:.3f} sec.".format(tm.time() - start_time))

        self.ostream.print_blank()

        tmp_start_time = tm.time()
        ct_matrix = self.compute_ct_matrix(molecule, basis, scf_results,
                                           tdens_ao, self.fragment_dict)
        ret_dict['ct_matrix'] = ct_matrix
        self.ostream.print_info("Charge transfer matrix computed in " +
                                "{:.3f} sec.".format(tm.time() -
                                                     tmp_start_time))

        self.ostream.print_blank()

        tmp_start_time = tm.time()
        pr_dict = self.compute_participation_ratios(ct_matrix)
        ret_dict.update(pr_dict)
        self.ostream.print_info("Participation ratios computed in " +
                                "{:.3f} sec.".format(tm.time() -
                                                     tmp_start_time))

        self.ostream.print_blank()

        tmp_start_time = tm.time()
        avg_pos_dict = self.compute_avg_position(molecule, ct_matrix,
                                                 self.fragment_dict)
        ret_dict.update(avg_pos_dict)
        self.ostream.print_info("Average positions and ct length computed in " +
                                "{:.3f} sec.".format(tm.time() -
                                                     tmp_start_time))

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

        scf_results = {}
        rsp_results = {}
        for key in h5f:
            if key == "atom_coordinates":
                scf_results[key] = np.array(h5f[key])
            if key == "nuclear_charges":
                scf_results[key] = np.array(h5f[key])
            if key == "basis_set":
                scf_results[key] = np.array(h5f[key])
            if key == "scf":
                scf_results_dict = dict(h5f.get(key))
                for scf_key in scf_results_dict:
                    scf_results[scf_key] = np.array(scf_results_dict[scf_key])
            elif key == "rsp":
                rsp_results_dict = dict(h5f.get(key))
                for rsp_key in rsp_results_dict:
                    rsp_results[rsp_key] = np.array(rsp_results_dict[rsp_key])
        h5f.close()

        return scf_results, rsp_results

    def format_rsp_results(self, rsp_results):
        """
        Reads the eigenvectors in rsp_results and adds them
        as numpy arrays with individual keys to the dictionary.

        :param rsp_results:
            The dictionary containing the results from a
            rsp calculation.

        :return:
            The rsp results dictionary in the format needed for
            the excited state analysis driver class.
        """

        if 'eigenvectors_distributed' in rsp_results.keys():
            rsp_drv = LinearResponseEigenSolver()
            num_states = len(rsp_results['eigenvectors_distributed'])
            for i in range(num_states):
                name = 'S' + str(i + 1)
                rsp_results[name] = rsp_drv.get_full_solution_vector(
                    rsp_results['eigenvectors_distributed'][i])
        else:
            num_states = len(rsp_results['eigenvectors'][0, :])
            for i in range(num_states):
                name = 'S' + str(i + 1)
                rsp_results[name] = rsp_results['eigenvectors'][:, i]
        rsp_results['formatted'] = True

        return rsp_results

    def create_molecule_and_basis(self, scf_results):
        """
        Creates molecule and basis objects from scf dictionary.

        :param scf_results:
            The dictionary containing the scf results dictionary.

        :return:
            A tuple containing the Molecule and Basis objects.
        """

        coordinates = scf_results['atom_coordinates']
        nuclear_charges = np.array(scf_results['nuclear_charges'])
        nuclear_charges = nuclear_charges.astype(int)
        basis_set_label = scf_results['basis_set'][0].decode("utf-8")
        molecule = Molecule(nuclear_charges, coordinates, units="au")
        basis = MolecularBasis.read(molecule, basis_set_label)

        return molecule, basis

    def compute_density_matrices(self, molecule, scf_results, rsp_results,
                                 state_index):
        """
        Computes the transition density matrix (TDM) for state
        nstate in MO and AO basis.

        :param molecule:
            The Molecule object
        :param scf_results:
            The dictionary containing the scf results.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param state_index:
            The excited state for which the TDM is computed.

        :return:
            A tuple containing the transition, particle, and hole
            density matrices in MO and AO basis.
        """

        if any(key.startswith('eigenvector') for key in rsp_results.keys()
              ) and 'formatted' not in rsp_results.keys():
            rsp_results = self.format_rsp_results(rsp_results)

        mo = scf_results["C_alpha"]
        nocc = molecule.number_of_alpha_electrons()
        norb = mo.shape[1]
        nvirt = norb - nocc
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()

        excstate = "S" + str(state_index)
        eigvec = rsp_results[excstate]

        nexc = nocc * nvirt
        z_mat = eigvec[:nexc]
        if eigvec.shape[0] == nexc:
            tdens_mo = np.reshape(z_mat, (nocc, nvirt))
        else:
            y_mat = eigvec[nexc:]
            tdens_mo = np.reshape(z_mat - y_mat, (nocc, nvirt))

        tdens_ao = np.linalg.multi_dot([mo_occ, tdens_mo, mo_vir.T])

        hole_dens_mo = -np.matmul(tdens_mo, tdens_mo.T)
        hole_dens_ao = np.linalg.multi_dot([mo_occ, hole_dens_mo, mo_occ.T])

        part_dens_mo = np.matmul(tdens_mo.T, tdens_mo)
        part_dens_ao = np.linalg.multi_dot([mo_vir, part_dens_mo, mo_vir.T])

        return {
            'transition_density_matrix_MO': tdens_mo,
            'transition_density_matrix_AO': tdens_ao,
            'particle_density_matrix_MO': part_dens_mo,
            'particle_density_matrix_AO': part_dens_ao,
            'hole_density_matrix_MO': hole_dens_mo,
            'hole_density_matrix_AO': hole_dens_ao,
        }

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
                # i is 1-based index
                ao_map_fragments[fragment] += ao_list[i - 1]
        return ao_map_fragments

    def compute_ct_matrix(self, molecule, basis, scf_results, T, fragment_dict):
        """
        Computes the charge transfer (CT) matrix from the
        transition density matrix (TDM).

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param T:
            The transition density matrix in AO basis.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            The CT matrix.
        """

        ao_map_fragments = self.map_fragments_to_atomic_orbitals(
            molecule, basis, fragment_dict)
        S = scf_results["S"]

        TS = np.matmul(T, S)
        ST = np.matmul(S, T)
        STS = np.matmul(S, TS)
        TS_ST = TS * ST  # element by element
        T_STS = T * STS

        n_fragments = len(fragment_dict.keys())
        ct_matrix = np.zeros((n_fragments, n_fragments))

        for a, A in enumerate(fragment_dict.keys()):
            for b, B in enumerate(fragment_dict.keys()):
                TSST_AB = TS_ST[ao_map_fragments[A]][:, ao_map_fragments[B]]
                T_STS_AB = T_STS[ao_map_fragments[A]][:, ao_map_fragments[B]]
                ct_matrix[a, b] = np.einsum("mn->", TSST_AB + T_STS_AB)

        return 0.5 * ct_matrix

    def compute_participation_ratios(self, ct_matrix):
        """
        Computes the particle, hole, and average participation ratios
        from the charge transfer matrix.

        :param ct_matrix:
            The charge transfer matrix.

        :return:
            A dictionary containing the hole, particle,
            and average participation ratios.
        """

        omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / omega
        particle_weight = np.sum(ct_matrix, axis=0) / omega
        hole_pr = 1 / np.sum(hole_weight**2)
        particle_pr = 1 / np.sum(particle_weight**2)
        avg_pr = (hole_pr + particle_pr) / 2

        return {
            'hole_participation_ratio': hole_pr,
            'particle_participation_ratio': particle_pr,
            'avg_participation_ratio': avg_pr,
        }

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
            A dictionary containing the average particle position,
            average hole position, the difference vector
            from average hole to average particle,
            and the charge transfer length.
        """

        avg_coord_list = np.array(
            self.compute_fragment_coordinates(molecule, fragment_dict))

        omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / omega
        particle_weight = np.sum(ct_matrix, axis=0) / omega

        avg_hole_position = np.matmul(hole_weight,
                                      avg_coord_list) * bohr_in_angstrom()
        avg_particle_position = np.matmul(particle_weight,
                                          avg_coord_list) * bohr_in_angstrom()
        avg_diff_vec = avg_particle_position - avg_hole_position
        ct_length = np.linalg.norm(avg_diff_vec)

        return {
            'avg_hole_position': avg_hole_position,
            'avg_particle_position': avg_particle_position,
            'avg_difference_vector': avg_diff_vec,
            'ct_length': ct_length,
        }

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

    def show_density(self, molecule, basis, descriptors):
        """Displays the particle and hole densities of molecule.

        :param molecule:
            The molecule object.
        :param basis:
            The basis object.
        :param descriptors:
            Either descriptor dictionary for one excited state or
            dictionary containing the descriptor dictionaries of
            multiple excited states.
        """

        dens_viewer = DensityViewer()
        dens_viewer_dict = {}
        if 'hole_density_matrix_AO' in descriptors.keys():
            dens_viewer_dict['particle'] = descriptors[
                'particle_density_matrix_AO']
            dens_viewer_dict['hole'] = descriptors['hole_density_matrix_AO']
        else:
            for key in descriptors.keys():
                dens_viewer_dict[key + ' particle'] = descriptors[key][
                    'particle_density_matrix_AO']
                dens_viewer_dict[
                    key + ' hole'] = descriptors[key]['hole_density_matrix_AO']
        dens_viewer.plot(molecule, basis, dens_viewer_dict)

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
            if len(self.fragment_dict[key]) >= 7:
                # Splice atomic index list if longer than 7 indices
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

        self.ostream.print_blank()

        self.ostream.print_reference(self.get_reference())
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

    def print_charge_transfer_length(self, descriptor_dict):
        """
        Prints charge transfer length.

        :param descriptor_dict:
            The dictionary of descriptors.
        """
        rel_ct_length = descriptor_dict['ct_length']

        valstr = "Charge Transfer Length (Angstrom)"
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = 'Distance:            '
        valstr += '{:13.6f}'.format(rel_ct_length)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def print_results(self, descriptor_dict):
        """
        Prints results of analysis.

        :param descriptor_dict:
            The dictionary of descriptors.
        """

        self.print_avg_particle_hole_positions(descriptor_dict)
        self.print_participation_ratios(descriptor_dict)
        self.print_charge_transfer_length(descriptor_dict)

    def get_reference(self):
        """
        Gets reference string for excited state descriptors.
        """

        ref_str = "F. Plasser, M. Wormit, A. Dreuw, "
        ref_str += "J. Chem. Phys., 2014, 141, 024106"

        return ref_str

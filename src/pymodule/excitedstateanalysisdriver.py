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
from matplotlib import colormaps as cmaps
from matplotlib import pyplot as plt
from .errorhandler import assert_msg_critical
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check)
from .visualizationdriver import VisualizationDriver

class ExcitedStateAnalysisDriver:
    """
    Implements excited state descriptors at TDDFT level of theory.
    """

    def __init__(self):
        """
        Initialize excited state analysis driver.
        """
        pass

    def compute(self, fragment_dict, nstate, molecule=None, basis=None, scf_tensors=None, 
                rsp_tensors=None, scf_checkpoint_file = None,
                rsp_checkpoint_file=None):
        """
        Computes dictionary containing excited state descriptors for
        excited state nstate.

        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.
        :param nstate:
            The excited state for which the descriptors are computed.
        :param molecule:
            The Molecule object
        :param basis:
            The Basis object
        :param scf_tensors:
            The dictionary containing the scf tensors.
        :param rsp_tensors:
            The dictionary containing the rsp tensors.
        :param scf_checkpoint_file:
            The checkpoint file for the scf calculation.
        :param rsp_checkpoint_file:
            The checkpoint file for the rsp calculation.
        
        :return:
            The dictionary containing the descriptors.
        """
        
        if scf_checkpoint_file is not None and rsp_checkpoint_file is not None:
            checkpoint_files = True
            assert_msg_critical(molecule is None and basis is None and
                                scf_tensors is None and rsp_tensors is None,
                                'molecule/basis/scf_tensors/rsp_tensors not '+
                                'from checkpoint file found')
            
        elif scf_checkpoint_file is None and rsp_checkpoint_file is None:
            checkpoint_files = False
            assert_msg_critical(molecule is not None and basis is not None and
                                scf_tensors is not None and rsp_tensors is not None,
                                'molecule/basis/scf_tensors/rsp_tensors not found')

        if checkpoint_files == True:
            scf_tensors = self.read_scf_checkpoint_file(scf_checkpoint_file)
            rsp_tensors = self.read_rsp_checkpoint_file(rsp_checkpoint_file)
            molecule, basis = self.create_molecule_basis(scf_tensors)


        T_MO, T_AO = self.compute_T(molecule, scf_tensors, rsp_tensors, nstate)

        T_hole_MO, T_hole_AO = self.compute_T_hole(molecule, scf_tensors,
                                                        rsp_tensors, nstate)

        T_particle_MO, T_particle_AO = self.compute_T_particle(molecule, scf_tensors, 
                                                               rsp_tensors, nstate)

        ct_matrix = self.compute_ct_matrix(molecule, basis, scf_tensors,
                                                       T_AO, fragment_dict)

        hole_participation_ratio, particle_participation_ratio, average_participation_ratio = self.compute_participation_ratios(ct_matrix)

        average_particle_position, average_hole_position, average_difference_vector = self.compute_avg_position(molecule, basis, ct_matrix,
                                                                                                                    fragment_dict)
        relative_ct_length = self.compute_relative_ct_length(molecule, average_difference_vector)

        ret_dict = {
                    'T_MO': T_MO,
                    'T_AO': T_AO,
                    'T_hole_MO': T_hole_MO,
                    'T_hole_AO': T_hole_AO,
                    'T_particle_MO': T_particle_MO,
                    'T_particle_AO': T_particle_AO,
                    'CT_matrix': ct_matrix,
                    'Hole_participation_ratio': hole_participation_ratio,
                    'Particle_participation_ratio': particle_participation_ratio,
                    'Average_participation_ratio': average_participation_ratio,
                    'Average_hole_position': average_hole_position,
                    'Average_particle_position': average_particle_position,
                    'Average_difference_vector': average_difference_vector,
                    'Relative_ct_length': relative_ct_length,
                    }
        
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
        
        assert_msg_critical('atom_coordinates' in dict, 
                            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: '
                            + 'atom_coordinates not found')
        
        assert_msg_critical('nuclear_charges' in dict, 
                            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: '
                            + 'nuclear_charges not found')
        
        assert_msg_critical('basis_set' in dict, 
                            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: '
                            + 'basis_set not found')
        
        assert_msg_critical('C_alpha' in dict, 
                            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: '
                            + 'C_alpha not found')
        
        assert_msg_critical('S' in dict, 
                            'ExcitedStateAnalysisDriver.read_scf_checkpoint_file: '
                            + 'S not found')

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
    
    def compute_T(self, molecule, scf_tensors, rsp_tensors, nstate):
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

        C = scf_tensors["C_alpha"]
        nocc = molecule.number_of_alpha_electrons()
        norb = C.shape[1]
        nvirt = norb-nocc
    
        excstate = "S"+str(nstate)
        Xf = rsp_tensors[excstate]

        nexc = nocc * nvirt
        Zf = Xf[:nexc]
        if Xf.shape[0] == nexc:
            T_MO = np.reshape(Zf, (nocc, nvirt))
        else:
            Yf = Xf[nexc:]
            T_MO = np.reshape(Zf - Yf, (nocc, nvirt))

        T_AO = np.linalg.multi_dot([C[:, :nocc], T_MO, C[:, nocc:].T])

        return T_MO, T_AO
    
    def compute_T_hole(self, molecule, scf_tensors, rsp_tensors, nstate):
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

        C = scf_tensors["C_alpha"]
        nocc = molecule.number_of_alpha_electrons()
        norb = C.shape[1]
        nvirt = norb-nocc
    
        excstate = "S"+str(nstate)
        Xf = rsp_tensors[excstate]

        nexc = nocc * nvirt
        Zf = Xf[:nexc]
        if Xf.shape[0] == nexc:
            T_MO = np.reshape(Zf, (nocc, nvirt))
        else:
            Yf = Xf[nexc:]
            T_MO = np.reshape(Zf - Yf, (nocc, nvirt))

        T_hole_MO = np.matmul(T_MO, T_MO.T)
        T_hole_AO = np.linalg.multi_dot([C[:, :nocc], T_hole_MO, C[:, :nocc].T])

        return T_hole_MO, T_hole_AO

    def compute_T_particle(self, molecule, scf_tensors, rsp_tensors, nstate):
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

        C = scf_tensors["C_alpha"]
        nocc = molecule.number_of_alpha_electrons()
        norb = C.shape[1]
        nvirt = norb-nocc
    
        excstate = "S"+str(nstate)
        Xf = rsp_tensors[excstate]

        nexc = nocc * nvirt
        Zf = Xf[:nexc]
        if Xf.shape[0] == nexc:
            T_MO = np.reshape(Zf, (nocc, nvirt))
        else:
            Yf = Xf[nexc:]
            T_MO = np.reshape(Zf - Yf, (nocc, nvirt))

        T_particle_MO = np.matmul(T_MO.T, T_MO)
        T_particle_AO = np.linalg.multi_dot([C[:, nocc:], T_particle_MO, C[:, nocc:].T])

        return T_particle_MO, T_particle_AO
    
    def map_fragments_to_atomic_orbitals(self, molecule, basis, fragment_dict):
        """
        Maps the atom indices in the fragment dictionary to indices of atomic
        orbitals centered on the corresponding atom.

        :param molecule:
            The Molecule object.
        :param basis:
            The Basis object.
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
                ao_map_fragments[fragment] += ao_list[i - 1]  # i is 1-based index
        return ao_map_fragments
    
    def compute_ct_matrix(self, molecule, basis, scf_tensors, T, fragment_dict):
        """
        Computes the charge transfer (CT) matrix from the 
        transition density matrix (TDM).

        :param molecule:
            The Molecule object.
        :param basis:
            The Basis object.
        :param scf_tensors:
            The dictionary containing the scf tensors.
        :param T: 
            The TDM in AO basis.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            The CT matrix.
        """

        ao_map_fragments = self.map_fragments_to_atomic_orbitals(molecule, basis, fragment_dict)
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

        Omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / Omega
        particle_weight = np.sum(ct_matrix, axis=0) / Omega
        hole_PR = 1 / np.sum(hole_weight ** 2)
        particle_PR = 1 / np.sum(particle_weight ** 2)
        avg_PR = (hole_PR + particle_PR) / 2

        return hole_PR, particle_PR, avg_PR
    
    def compute_fragment_coordinates(self, molecule, fragment_dict):
        """
        Computes average coordinates of fragments in fragment_dict.

        :param molecule:
            The Molecule object.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            A list containing the average coordinates of each fragment.
        """

        coordinates = molecule.get_coordinates_in_bohr()
        avg_coord_list = []
        for i, A in enumerate(fragment_dict):
            coord_array = coordinates[np.array(fragment_dict[A]) - 1, :]
            avg_coord_list.append(np.sum(coord_array, axis=0) / coord_array.shape[0])
        return avg_coord_list
    
    def max_dist_from_mol_center(self, molecule):
        """
        Computes the maximum distance from an atom in the molecule to the 
        center of the molecule.

        :param molecule:
            The Molecule object.

        :return:
            The maximum distance from an atom in the molecule to the
            center of the molecule.
        """

        coordinates = molecule.get_coordinates_in_angstrom()
        center = np.sum(coordinates,axis=0)/coordinates.shape[0]
        print(center)
        atom_center_distance = []
        for i in range(coordinates.shape[0]):
            atom_center_distance.append(np.linalg.norm(coordinates[i]-center))
        return np.max(np.array(atom_center_distance))
    
    def compute_avg_position(self, molecule, ct_matrix, fragment_dict):
        """
        Computes the average particle position, average hole position, and the 
        the difference vector between them.

        :param molecule:
            The Molecule object.
        :param ct_matrix:
            The charge transfer matrix.
        :param fragment_dict:
            The dictionary containing the indices of atoms in each fragment.

        :return:
            A tuple containing the average particle position, 
            average hole position, and the difference vector 
            from average hole to average particle.
        """

        avg_coord_list = self.compute_fragment_coordinates(molecule, fragment_dict)
    
        Omega = np.sum(ct_matrix)
        hole_weight = np.sum(ct_matrix, axis=1) / Omega
        particle_weight = np.sum(ct_matrix, axis=0) / Omega
    
        avg_hole_pos = [0, 0, 0]
        avg_particle_pos = [0, 0, 0]
    
        for j in range(len(avg_coord_list)):
            avg_hole_pos += (
                avg_coord_list[j]
                * hole_weight[j]
                * bohr_in_angstrom()
            )
            avg_particle_pos += (
                avg_coord_list[j]
                * particle_weight[j]
                * bohr_in_angstrom()
            )
        ET_net = avg_particle_pos - avg_hole_pos
        return avg_particle_pos, avg_hole_pos, ET_net
    
    def compute_relative_ct_length(self, molecule, ET_net):
        """
        Compute the relative charge transfer (CT) length.

        :param molecule:
            The Molecule object.
        :param ET_net:
            The difference vector from average position of hole to average
            position of particle.

        :return:
            The relative charge transfer length.
        """
        return np.linalg.norm(ET_net)/self.max_dist_from_mol_center(molecule)
    
    def compute_avg_ph_viewer(self, molecule, avg_hole_pos, avg_particle_pos, width=400, height=300):
        """
        Displays the molecule and the average positions of the particle and hole.

        :param molecule:
            The Molecule object.
        :param avg_hole_pos:
            The average position of the hole.
        :param avg_particle_pos:
            The average positions of the particle.
        :param width:
            The width of the plot.
        :param height:
            The height of the plot.
        """
        r_H = avg_hole_pos
        r_E = avg_particle_pos
        viewer = p3d.view(width=width, height=height)
        viewer.addModel(molecule.get_xyz_string(), "xyz")
        viewer.setStyle({"line": {}, "sphere": {"scale": 0.05}})
        viewer.addArrow({
                'start': {'x': r_H[0], 'y': r_H[1], 'z': r_H[2]},
                'end': {'x': r_E[0], 'y': r_E[1], 'z': r_E[2]},
                'radius': 0.06,
                'color': 'green',
            })
        viewer.addSphere({
            'center': {'x': r_H[0], 'y': r_H[1], 'z': r_H[2]},
            'radius': 0.12,
            'color': 'blue',
            })
        viewer.addSphere({
            'center': {'x': r_E[0], 'y': r_E[1], 'z': r_E[2]},
            'radius': 0.12,
            'color': 'red',
            })
        viewer.show()
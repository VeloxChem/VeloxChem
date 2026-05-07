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

import sys
from pathlib import Path

import numpy as np
import networkx as nx
from mpi4py import MPI
import h5py
import re

try:
    from scipy.optimize import minimize
except ImportError:
    pass

from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...molecule import Molecule

from ..io.basic import nn, nl, pname, is_list_A_in_B
from ..utils.geometry import (unit_cell_to_cartesian_matrix,
                              fractional_to_cartesian, cartesian_to_fractional,
                              locate_min_idx, reorthogonalize_matrix,
                              find_optimal_pairings, find_edge_pairings)
from .other import fetch_X_atoms_ind_array
from .superimpose import superimpose_rotation_only


class NetOptimizer:
    """Optimizes node rotations and cell parameters so linkers fit the net, then places edges.

    Uses OptimizationDriver for two-stage rotation optimization and cell scaling.
    Requires G, V_data, V_X_data, E_data, E_X_data, sorted_nodes, sorted_edges,
    cell_info, and linker_frag_length (and optionally EC_data for multitopic). Sets
    opt_rots, sc_unit_cell, sc_rot_node_X_pos, sG, optimized_pair, etc.
    """

    def __init__(self, comm=None, ostream=None):
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        #NEED to be set before use
        self.G = None
        self.V_data = None
        self.V_X_data = None
        self.EC_data = None
        self.EC_X_data = None
        self.E_data = None
        self.E_X_data = None
        self.sorted_nodes = None
        self.sorted_edges = None
        self.cell_info = None

        #will be generated
        self.sorted_edges_of_sortednodeidx = None
        self.pname_set_dict = None
        self.pname_set = None
        self.opt_rots = None
        self.opt_params = None
        self.new_edge_length = None
        self.optimized_pair = None
        self.sc_node_pos_dict = None
        self.sc_rot_node_X_pos = None
        self.sc_unit_cell = None
        self.sc_unit_cell_inv = None
        self.fake_edge = False
        #self.constant_length = 1.54  #default C-C single bond length
        self.linker_frag_length = None

        self.load_optimized_rotations = None
        self.skip_rotation_optimization = False
        self.rotation_filename = None

        # Optimization parameters
        self.opt_drv = OptimizationDriver(comm=self.comm, ostream=self.ostream)
        self.opt_drv.pname_set_dict = None  #to be set later
        #self.opt_drv.opt_method = 'L-BFGS-B'
        #self.opt_drv.maxfun = 15000
        #self.opt_drv.maxiter = 15000
        #self.opt_drv.display = True
        #self.opt_drv.eps = 1e-8

        #other parameters
        self._debug = False
        self.opt_drv._debug = self._debug

        self.optimized_cell_info = None
        """
        - node_target_type (str):metal atom type of the node
        - node_unit_cell (array):unit cell of the node
        - node_atom (array):2 columns, atom_name, atom_type of the node
        - node_x_fcoords (array):fractional coordinates of the X connected atoms of node
        - node_fcoords (array):fractional coordinates of the whole node
        - node_x_ccoords (array):cartesian coordinates of the X connected atoms of node
        - node_coords (array):cartesian coordinates of the whole node
        - linker_unit_cell (array):unit cell of the ditopic linker or branch of multitopic linker
        - linker_atom (array):2 columns, atom_name, atom_type of the ditopic linker or branch of multitopic linker
        - linker_x_fcoords (array):fractional coordinates of the X connected atoms of ditopic linker or branch of multitopic linker
        - linker_fcoords (array):fractional coordinates of the whole ditopic linker or branch of multitopic linker
        - linker_x_ccoords (array):cartesian coordinates of the X connected atoms of ditopic linker or branch of multitopic linker
        - linker_frag_length (float):distance between two X-X connected atoms of the ditopic linker or branch of multitopic linker
        - linker_ccoords (array):cartesian coordinates of the whole ditopic linker or branch of multitopic linker
        - linker_center_cif (str):cif file of the center of the multitopic linker
        - ec_unit_cell (array):unit cell of the center of the multitopic linker
        - ec_atom (array):2 columns, atom_name, atom_type of the center of the multitopic linker
        - ec_x_vecs (array):fractional coordinates of the X connected atoms of the center of the multitopic linker
        - ec_fcoords (array):fractional coordinates of the whole center of the multitopic linker
        - ec_xcoords (array):cartesian coordinates of the X connected atoms of the center of the multitopic linker
        - eccoords (array):cartesian coordinates of the whole center of the multitopic linker
        - constant_length (float):constant length to add to the linker length, normally 1.54 for single bond of C-C,
        because C is always used as the connecting atom in the builder
        - maxfun (int):maximum number of function evaluations for the node rotation optimization
        - opt_method (str):optimization method for the node rotation optimization
        - G (networkx graph):graph of the template
        - node_max_degree (int):maximum degree of the node in the template, should be the same as the node topic
        - sorted_nodes (list):sorted nodes in the template by connectivity
        - sorted_edges (list):sorted edges in the template by connectivity
        - sorted_edges_of_sortednodeidx (list):sorted edges in the template by connectivity with the index of the sorted nodes
        - optimized_rotations (dict):optimized rotations for the nodes in the template
        - optimized_params (array):optimized cell parameters for the template topology to fit the target MOF cell
        - new_edge_length (float):new edge length of the ditopic linker or branch of multitopic linker, 2*constant_length+linker_frag_length
        - optimized_pair (dict): pair of connected nodes in the template with the index of the X connected atoms, used for the edge placement
        - scaled_rotated_node_positions (dict):scaled and rotated node positions in the target MOF cell
        - scaled_rotated_Xatoms_positions (dict):scaled and rotated X connected atom positions of nodes in the target MOF cell
        - sc_unit_cell (array):(scaled) unit cell of the target MOF cell
        - sc_unit_cell_inv (array):inverse of the (scaled) unit cell of the target MOF cell
        - sG_node (networkx graph):graph of the target MOF cell
        - nodes_atom (dict):atom name and atom type of the nodes
        - rotated_node_positions (dict):rotated node positions
        - supercell (array): supercell set by user, along x,y,z direction
        - multiedge_bundlings (dict):multiedge bundlings of center and branches of the multitopic linker, used for the supercell
        construction and merging of center and branches to form one EDGE
        - prim_multiedge_bundlings (dict):multiedge bundlings in primitive cell, used for the supercell construction
        - super_multiedge_bundlings (dict):multiedge bundlings in the supercell, used for the supercell construction
        - dv_v_pairs (dict):DV and V pairs in the template, used for the supercell construction
        - super_multiedge_bundlings (dict):multiedge bundlings in the supercell, used for the supercell construction
        - superG (networkx graph):graph of the supercell
        - add_virtual_edge (bool): add virtual edge to the target MOF cell
        - vir_edge_range (float): range to search the virtual edge between two Vnodes directly, should <= 0.5,
        used for the virtual edge addition of bridge type nodes: nodes and nodes can connect directly without linker
        - vir_edge_max_neighbor (int): maximum number of neighbors of the node with virtual edge, used for the virtual edge addition of bridge type nodes
        - remove_node_list (list):list of nodes to remove in the target MOF cell
        - remove_edge_list (list):list of edges to remove in the target MOF cell
        - eG (networkx graph):graph of the target MOF cell with only EDGE and V nodes
        - node_topic (int):maximum degree of the node in the template, should be the same as the node_max_degree
        - unsaturated_node (list):unsaturated nodes in the target MOF cell
        - term_info (array):information of the node terminations
        - term_coords (array):coordinates of the node terminations
        - term_xoovecs (array):X and O vectors (usually carboxylate group) of the node terminations
        - unsaturated_vnode_xind_dict (dict):unsaturated node and the exposed X connected atom index
        - unsaturated_vnode_xoo_dict (dict):unsaturated node and the exposed X connected atom index and the corresponding O connected atoms
        """

    def rotation_and_cell_optimization(self):
        """
        two optimization steps:
        1. optimize the node rotation (vertex or edge center)
        2. optimize the cell parameters to fit the target MOF cell
        """
        #if self._debug:
        self.ostream.print_info(f"constant_length: {self.constant_length}")
        self.ostream.flush()

        G = self.G.copy()
        self.node_x_ccoords = self.V_X_data[:, 5:8].astype(float)
        self.node_ccoords = self.V_data[:, 5:8].astype(float)
        self.ec_x_ccoords = self.EC_X_data[:, 5:8].astype(
            float) if self.EC_X_data is not None else None
        self.ec_ccoords = self.EC_data[:, 5:8].astype(
            float) if self.EC_data is not None else None
        self.e_x_ccoords = self.E_X_data[:, 5:8].astype(float)
        self.e_ccoords = self.E_data[:, 5:8].astype(float)
        self.opt_drv.sorted_nodes = self.sorted_nodes
        self.opt_drv.pname_set_dict = self.pname_set_dict

        node_atom = self.V_data[:, 0:2]
        ec_atom = self.EC_data[:, 0:2] if self.EC_data is not None else None

        linker_frag_length = self.linker_frag_length
        constant_length = self.constant_length
        x_com_length = np.mean(
            [np.linalg.norm(i) for i in self.node_x_ccoords])
        # firstly, check if all V nodes have highest connectivity
        # secondly, sort all DV nodes by connectivity
        sorted_nodes = self.sorted_nodes
        sorted_edges = self.sorted_edges

        self.nodes_atom = {}
        for n in sorted_nodes:
            if "CV" in n:
                self.nodes_atom[n] = ec_atom
            else:
                self.nodes_atom[n] = node_atom

    # reindex the nodes in the Xatoms_positions with the index in the sorted_nodes, like G has 16 nodes[2,5,7], but the new dictionary should be [0,1,2]
        node_pos_dict, node_X_pos_dict = self._generate_pos_dict(G)

        # reindex the edges in the G with the index in the sorted_nodes
        self.sorted_edges_of_sortednodeidx = [(sorted_nodes.index(e[0]),
                                               sorted_nodes.index(e[1]))
                                              for e in self.sorted_edges]
        self.opt_drv.sorted_edges = self.sorted_edges_of_sortednodeidx
        # Optimize rotations
        num_nodes = G.number_of_nodes()
        pname_set, pname_set_dict = self._generate_pname_set(
            G, sorted_nodes, node_X_pos_dict)
        self.pname_set_dict = pname_set_dict
        self.opt_drv.pname_set_dict = pname_set_dict

        node_pos_dict, node_X_pos_dict = self._apply_rot_trans2dict(
            G, node_pos_dict, node_X_pos_dict)
        ###3D free rotation
        saved_optimized_rotations = None

        ini_rot = (np.eye(3, 3).reshape(1, 3, 3).repeat(len(pname_set),
                                                        axis=0))

        if self.load_optimized_rotations is not None and Path(
                self.load_optimized_rotations).is_file():
            #load the saved optimized rotations
            with h5py.File(self.load_optimized_rotations, 'r') as hf:
                saved_optimized_rotations = hf['optimized_rotations'][:]
                if self._debug:
                    self.ostream.print_info(
                        f"Loaded saved optimized rotations from {self.load_optimized_rotations}"
                    )
                    self.ostream.print_info(
                        f"Loaded saved optimized rotations shape: {saved_optimized_rotations.shape}"
                    )
            if not self.skip_rotation_optimization:  #load but not skip
                self.ostream.print_info(
                    "use the loaded optimized_rotations from the previous optimization as initial guess"
                )
                ini_rot = saved_optimized_rotations.reshape(-1, 3, 3)

        if saved_optimized_rotations is None:
            self.skip_rotation_optimization = False

        if not self.skip_rotation_optimization:
            ####TODO: modified for mil53
            opt_rot_pre, _ = self.opt_drv._optimize_rotations_pre(
                num_nodes, G, node_X_pos_dict, ini_rot)
            opt_rot_aft, _ = self.opt_drv._optimize_rotations_after(
                num_nodes, G, node_X_pos_dict, opt_rot_pre)
        else:
            opt_rot_aft = saved_optimized_rotations.reshape(-1, 3, 3)

        if self.rotation_filename is not None:
            #save it as h5 file
            with h5py.File(self.rotation_filename, 'w') as hf:
                hf.create_dataset('optimized_rotations', data=opt_rot_aft)

        opt_rots = expand_set_rots(pname_set_dict, opt_rot_aft, sorted_nodes)
        # Apply rotations
        rot_node_pos, _ = self._apply_rot2atoms_pos(opt_rots, G, node_pos_dict)

        if self._debug:
            self.ostream.print_info(
                f"Optimized Rotations (after optimization): {opt_rots}")

            def temp_save_xyz(filename, rot_node_pos):
                with open(filename, "w") as file:
                    num_atoms = sum(
                        len(positions) for positions in rot_node_pos.values())
                    file.write(f"{num_atoms}\n")
                    file.write("Optimized structure\n")
                    for node, positions in rot_node_pos.items():
                        for pos in positions:
                            file.write(
                                f"X{node}   {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n"
                            )

            temp_save_xyz("optimized_nodesstructure.xyz", rot_node_pos)

        rot_node_X_pos_dict, _ = self._apply_rot2atoms_pos(
            opt_rots, G, node_X_pos_dict)

        start_node = self.sorted_edges[0][
            0]  # find_nearest_node_to_beginning_point(G)

        # loop all of the edges in G and get the lengths of the edges, length is the distance between the two nodes ccoords
        edge_lengths, lengths = self._get_edge_lengths(G)

        x_com_length = np.mean(
            [np.linalg.norm(i) for i in self.node_x_ccoords])
        if self.fake_edge:
            new_edge_length = self.constant_length + 2 * x_com_length
        else:
            new_edge_length = self.linker_frag_length + 2 * self.constant_length + 2 * x_com_length

        # update the node ccoords in G by loop edge, start from the start_node, and then update the connected node ccoords by the edge length, and update the next node ccords from the updated node

        new_ccoords, old_ccoords = self._update_node_ccoords(
            G, edge_lengths, start_node, new_edge_length)
        # exclude the start_node in updated_ccoords and original_ccoords
        new_ccoords = {k: v for k, v in new_ccoords.items() if k != start_node}
        old_ccoords = {k: v for k, v in old_ccoords.items() if k != start_node}

        # use optimized_params to update all of nodes ccoords in G, according to the fccoords
        if self.optimized_cell_info is None:
            self.ostream.print_info("-" * 80)
            self.ostream.print_info(
                "Start to optimize the cell parameters to fit the target MOF cell"
            )
            self.ostream.print_info("-" * 80)

            optimized_cell_info = self.opt_drv._optimize_cell_params(
                self.cell_info, old_ccoords, new_ccoords)

            self.ostream.print_info("-" * 80)
            self.ostream.print_info(
                "Finished the cell parameters optimization")
            self.ostream.print_info("-" * 80)
        else:
            self.ostream.print_info(
                "use the optimized_cell_info from the previous optimization")
            optimized_cell_info = self.optimized_cell_info

        # get scaled unit cell and inverse and sG
        a, b, c, alpha, beta, gamma = optimized_cell_info
        sc_unit_cell = unit_cell_to_cartesian_matrix(a, b, c, alpha, beta,
                                                     gamma)
        sc_unit_cell_inv = np.linalg.inv(sc_unit_cell)
        sG, scaled_ccoords = self._update_ccoords_by_optimized_cell_params(
            G, optimized_cell_info)

        # update ccoords in sG
        sc_node_pos_dict, sc_node_X_pos_dict = self._generate_pos_dict(sG)

        # Apply rotations and translations to the scaled positions in sG
        sc_node_pos_dict, sc_node_X_pos_dict = self._apply_rot_trans2dict(
            sG, sc_node_pos_dict, sc_node_X_pos_dict)
        sc_rot_node_pos, _ = self._apply_rot2atoms_pos(opt_rots, sG,
                                                       sc_node_pos_dict)
        sc_rot_node_X_pos, self.optimized_pair = self._apply_rot2atoms_pos(
            opt_rots, sG, sc_node_X_pos_dict)

        # Save results to XYZ
        if self._debug:
            temp_save_xyz("scaled_optimized_nodesstructure.xyz",
                          sc_rot_node_pos)

        self.opt_rots = opt_rots
        self.optimized_cell_info = optimized_cell_info
        self.new_edge_length = new_edge_length

        self.sG = sG
        self.sc_node_pos_dict = sc_node_pos_dict
        self.sc_rot_node_X_pos = sc_rot_node_X_pos
        self.sc_rot_node_pos = sc_rot_node_pos
        self.sc_unit_cell = sc_unit_cell
        self.sc_unit_cell_inv = sc_unit_cell_inv

        self.rotated_node_positions = rot_node_pos
        self.Xatoms_positions_dict = node_X_pos_dict
        self.node_pos_dict = node_pos_dict

    def place_edge_in_net(self):
        """
        based on the optimized rotations and cell parameters, use optimized pair to find connected X-X pair in optimized cell,
        and place the edge in the target MOF cell

        return:
            sG (networkx graph):graph of the target MOF cell, with scaled and rotated node and edge positions
        """
        # linker_middle_point = np.mean(linker_x_vecs,axis=0)
        e_xx_vec = self.e_x_ccoords
        self.e_atom = self.E_data[:, 0:2]
        linker_frag_length = self.linker_frag_length
        optimized_pair = self.optimized_pair
        scaled_rotated_Xatoms_positions = self.sc_rot_node_X_pos
        scaled_rotated_node_positions = self.sc_rot_node_pos
        sorted_nodes = self.sorted_nodes
        sG = self.sG.copy()
        sc_unit_cell_inv = self.sc_unit_cell_inv
        nodes_atom = self.nodes_atom
        if linker_frag_length > 0.0:
            scalar = (linker_frag_length +
                      2 * self.constant_length) / linker_frag_length
        else:
            scalar = 1.0

        extended_e_xx_vec = [i * scalar for i in e_xx_vec]
        norm_xx_vector_record = []
        rot_record = []

        def correct_rot_direction(rot, source_vec, target_vec):
            if (np.linalg.norm(source_vec) == 0
                    or np.linalg.norm(target_vec) == 0):
                return rot

            rotated_source_vec = np.dot(source_vec, rot)
            if np.dot(rotated_source_vec, target_vec) >= 0:
                return rot

            axis = np.cross(rotated_source_vec, target_vec)
            if np.linalg.norm(axis) == 0:
                axis = np.cross(rotated_source_vec, [1, 0, 0])
                if np.linalg.norm(axis) == 0:
                    axis = np.cross(rotated_source_vec, [0, 1, 0])

            axis = axis / np.linalg.norm(axis)
            flip_matrix = -np.eye(3) + 2 * np.outer(axis, axis)
            return np.dot(rot, flip_matrix)

        # edges = {}
        for (i, j), pair in optimized_pair.items():
            x_idx_i, x_idx_j = pair
            reindex_i = sorted_nodes.index(i)
            reindex_j = sorted_nodes.index(j)
            x_i = scaled_rotated_Xatoms_positions[reindex_i][x_idx_i][1:]
            x_j = scaled_rotated_Xatoms_positions[reindex_j][x_idx_j][1:]
            x_i_x_j_middle_point = np.mean([x_i, x_j], axis=0)
            xx_vector = np.vstack(
                [x_i - x_i_x_j_middle_point, x_j - x_i_x_j_middle_point])
            norm_xx_vector = xx_vector / np.linalg.norm(xx_vector)

            if self._debug:
                self.ostream.print_info(
                    f"Placing edge between Node {i} and Node {j}")
                self.ostream.print_info(
                    f"  X atom index in Node {i}: {x_idx_i}, coordinates: {x_i}"
                )
                self.ostream.print_info(
                    f"  X atom index in Node {j}: {x_idx_j}, coordinates: {x_j}"
                )
                self.ostream.print_info(
                    f"  Middle point: {x_i_x_j_middle_point}")
                self.ostream.print_info(f"  Original XX vector: {xx_vector}")
                self.ostream.print_info(
                    f"  Normalized XX vector: {norm_xx_vector}")
                self.ostream.flush()
            # use superimpose to get the rotation matrix
            # use record to record the rotation matrix for get rid of the repeat calculation
            if self.linker_frag_length >= 0.0:
                # for normal linker, the direction is important
                indices = [
                    index for index, value in enumerate(norm_xx_vector_record)
                    if is_list_A_in_B(norm_xx_vector, value)
                ]
                if len(indices) == 1:
                    rot = rot_record[indices[0]]
                    # rot = reorthogonalize_matrix(rot)
                else:
                    _, rot, trans = superimpose_rotation_only(
                        extended_e_xx_vec, xx_vector)
                    # rot = reorthogonalize_matrix(rot)
                    norm_xx_vector_record.append(norm_xx_vector)
                    # the rot may be opposite, so we need to check the angle between the two vectors
                    # if the angle is larger than 90 degree, we need to reverse the rot
                    if self.EC_X_data is None:
                        roted_xx = np.dot(extended_e_xx_vec, rot)

                        if np.dot(roted_xx[1] - roted_xx[0],
                                  xx_vector[1] - xx_vector[0]) < 0:
                            ##rotate 180 around the axis of the cross product of the two vectors
                            axis = np.cross(roted_xx[1] - roted_xx[0],
                                            xx_vector[1] - xx_vector[0])
                            # if 001 not linear to the two vectors
                            if np.linalg.norm(axis) == 0:
                                check_z_axis = np.cross(
                                    roted_xx[1] - roted_xx[0], [0, 0, 1])
                                if np.linalg.norm(check_z_axis) == 0:
                                    axis = np.array([1, 0, 0])
                                else:
                                    axis = np.array([0, 0, 1])

                            axis = axis / np.linalg.norm(axis)
                            flip_matrix = np.eye(3) - 2 * np.outer(
                                axis, axis)  # Householder matrix for reflection
                            rot = np.dot(rot, flip_matrix)
                    # Flip the last column of the rotation matrix if the determinant is negative
                    rot_record.append(rot)
                if self.EC_X_data is not None:
                    source_vec = e_xx_vec[1] - e_xx_vec[0]
                    if "CV" in str(i):
                        target_vec = x_j - x_i
                    elif "CV" in str(j):
                        target_vec = x_i - x_j
                    else:
                        target_vec = xx_vector[1] - xx_vector[0]
                    rot = correct_rot_direction(rot, source_vec, target_vec)
            else:
                #get a random rotation matrix
                rot = np.eye(3)

            # use the rotation matrix to rotate the linker x coords
            placed_edge_ccoords = (np.dot(self.e_ccoords, rot) +
                                   x_i_x_j_middle_point)

            placed_edge = np.hstack(
                (np.asarray(self.e_atom), placed_edge_ccoords))
            sG.edges[(i, j)]["coords"] = x_i_x_j_middle_point
            sG.edges[(i, j)]["c_points"] = placed_edge

            sG.edges[(i, j)]["f_points"] = np.hstack((
                placed_edge[:, 0:2],
                cartesian_to_fractional(placed_edge[:, 2:5], sc_unit_cell_inv),
            ))

            _, sG.edges[(i, j)]["x_coords"] = fetch_X_atoms_ind_array(
                placed_edge, 0, "X")
        for i, v in scaled_rotated_node_positions.items():
            k = sorted_nodes[i]
            pos = v[:, 1:]  #cause first column is index added by addidx
            sG.nodes[k]["c_points"] = np.hstack((nodes_atom[k], pos))
            sG.nodes[k]["f_points"] = np.hstack(
                (nodes_atom[k], cartesian_to_fractional(pos,
                                                        sc_unit_cell_inv)))
            # find the atoms starts with "x" and extract the coordinates
            _, sG.nodes[k]["x_coords"] = fetch_X_atoms_ind_array(
                sG.nodes[k]["c_points"], 0, "X")
        self.sG = sG
        return sG

    def _get_edge_lengths(self, G):
        """Compute edge length (distance between node ccoords) for each edge. Returns (edge_lengths dict, set of lengths)."""
        edge_lengths = {}
        lengths = []
        for e in G.edges():
            i, j = e
            length = np.linalg.norm(G.nodes[i]["ccoords"] -
                                    G.nodes[j]["ccoords"])
            length = np.round(length, 3)
            edge_lengths[(i, j)] = length
            edge_lengths[(j, i)] = length
            lengths.append(length)
        if self._debug:
            self.ostream.print_info(f"Edge lengths: {edge_lengths}")
            self.ostream.print_info(
                f"Set of unique edge lengths: {set(lengths)}")
        if len(set(lengths)) != 1:
            self.ostream.print_warning(
                "Warning: more than one type of edge length")
            # if the length are close, which can be shown by std
            if np.std(lengths) < 1:  #1 Angstrom
                self.ostream.print_info("the edge lengths are close")
            else:
                self.ostream.print_info("the edge lengths are not close")
            self.ostream.print_info(str(set(lengths)))
        return edge_lengths, set(lengths)

    def _apply_rot2atoms_pos(self, optimized_rotations, G, node_X_pos_dict):
        """Apply optimized rotation matrices to node X positions and compute edge pairings.

        Args:
            optimized_rotations: Rotation matrix per node (by sorted_nodes index).
            G: Net graph with node "ccoords".
            node_X_pos_dict: Dict mapping node index to (N, 4) array [idx, x, y, z].

        Returns:
            Tuple (rotated_positions, optimized_pair): Rotated positions dict and edge (i,j) -> (x_idx_i, x_idx_j).
        """
        rotated_positions = node_X_pos_dict.copy()
        sorted_nodes = self.sorted_nodes
        sorted_edges_of_sortednodeidx = self.sorted_edges_of_sortednodeidx

        for i, node in enumerate(sorted_nodes):
            # if node type is V
            # if 'DV' in G.nodes[node]['type']:
            # continue
            R = optimized_rotations[i]

            original_positions = rotated_positions[i][:, 1:]
            com = G.nodes[node]["ccoords"]

            # Translate, rotate, and translate back to preserve the mass center
            translated_positions = original_positions - com
            rotated_translated_positions = np.dot(translated_positions, R.T)
            rotated_positions[i][:, 1:] = rotated_translated_positions + com
        edge_pair = find_edge_pairings(sorted_nodes,
                                       sorted_edges_of_sortednodeidx,
                                       rotated_positions)
        if self._debug:
            self.ostream.print_info(
                f"Optimized Pairings (after optimization): {edge_pair}")

        optimized_pair = {}
        for (i, j), pair in edge_pair.items():
            if self._debug:
                self.ostream.print_info(
                    f"Node {sorted_nodes[i]} and Node {sorted_nodes[j]}:")
            idx_i, idx_j = pair

            if self._debug:
                self.ostream.print_info(
                    f"Node{sorted_nodes[i]}_{int(idx_i)} -- Node{sorted_nodes[j]}_{int(idx_j)}"
                )
            optimized_pair[sorted_nodes[i],
                           sorted_nodes[j]] = (int(idx_i), int(idx_j))

        return rotated_positions, optimized_pair

    def _generate_pos_dict(self, sG):
        """Build dicts of node positions and X-atom positions per node index from sG and node/EC coords."""
        ec_x_ccoords = self.ec_x_ccoords
        ec_ccoords = self.ec_ccoords
        node_x_ccoords = self.node_x_ccoords
        node_ccoords = self.node_ccoords
        sorted_nodes = self.sorted_nodes

        def addidx(array):
            row_indices = np.arange(array.shape[0]).reshape(-1, 1).astype(int)
            new_array = np.hstack((row_indices, array))
            return new_array

        sc_node_pos_dict = {}
        sc_node_X_pos_dict = {}
        for n in sorted_nodes:
            if "CV" in n:
                sc_node_X_pos_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + ec_x_ccoords)
                sc_node_pos_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + ec_ccoords)
            else:
                sc_node_X_pos_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + node_x_ccoords)
                sc_node_pos_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + node_ccoords)
        return sc_node_pos_dict, sc_node_X_pos_dict

    def _apply_rot_trans2dict(self, sG, sc_node_pos_dict, sc_node_X_pos_dict):
        """Apply per-pname rotation and translation to node and X positions in place."""
        sorted_nodes = self.sorted_nodes
        pname_set_dict = self.pname_set_dict
        opt_rots = self.opt_rots
        for p_name in pname_set_dict:
            rot, trans = pname_set_dict[p_name]["rot_trans"]
            for k in pname_set_dict[p_name]["ind_ofsortednodes"]:
                node = sorted_nodes[k]
                sc_node_X_pos_dict[k][:, 1:] = (np.dot(
                    sc_node_X_pos_dict[k][:, 1:] - sG.nodes[node]["ccoords"],
                    rot,
                ) + trans + sG.nodes[node]["ccoords"])

                sc_node_pos_dict[k][:, 1:] = (np.dot(
                    sc_node_pos_dict[k][:, 1:] - sG.nodes[node]["ccoords"],
                    rot) + trans + sG.nodes[node]["ccoords"])

        return sc_node_pos_dict, sc_node_X_pos_dict

    def _generate_pname_set(self, G, sorted_nodes, node_X_pos_dict):
        """Build pname set and dict mapping pname to sorted node indices and initial rot_trans per pname."""
        pname_list = [pname(n) for n in sorted_nodes]
        pname_set = set(pname_list)
        pname_set_dict = {}
        for n in pname_set:
            pname_set_dict[n] = {
                "ind_ofsortednodes": [],
            }
        for i, node in enumerate(sorted_nodes):
            pname_set_dict[pname(node)]["ind_ofsortednodes"].append(i)
            if len(pname_set_dict[pname(node)]
                   ["ind_ofsortednodes"]) == 1:  # first node
                pname_set_dict[pname(
                    node)]["rot_trans"] = get_rot_trans_matrix(
                        node, G, sorted_nodes,
                        node_X_pos_dict)  # initial guess

        if self._debug:
            self.ostream.print_info(f"sorted_nodes: {sorted_nodes} ")
            self.ostream.print_info(f"pname_set: {pname_set}")
            self.ostream.print_info(f"pname_set_dict: {pname_set_dict}")
        return pname_set, pname_set_dict

    def _update_node_ccoords(self, G, edge_lengths, start_node,
                             new_edge_length):
        """Propagate node positions from start_node so edge lengths match new_edge_length."""
        updated_ccoords = {}
        original_ccoords = {}
        updated_ccoords[start_node] = G.nodes[start_node]["ccoords"]
        original_ccoords[start_node] = G.nodes[start_node]["ccoords"]
        updated_node = [start_node]
        # to update all the nodes start from the startnode and spread to neighbors
        for i in range(len(G.nodes()) - 1):
            for n in updated_node:
                for nn in G.neighbors(n):
                    if nn in updated_node:
                        continue
                    edge = (n, nn)
                    edge_length = edge_lengths[edge]
                    updated_ccoords[nn] = (
                        updated_ccoords[n] +
                        (G.nodes[nn]["ccoords"] - G.nodes[n]["ccoords"]) *
                        new_edge_length / edge_length)
                    original_ccoords[nn] = G.nodes[nn]["ccoords"]
                    updated_node.append(nn)

        return updated_ccoords, original_ccoords

    def _update_ccoords_by_optimized_cell_params(self, G, optimized_params):
        """Update node ccoords in G from fcoords using the optimized unit cell matrix."""
        sG = G.copy()
        a, b, c, alpha, beta, gamma = optimized_params
        T_unitcell = unit_cell_to_cartesian_matrix(a, b, c, alpha, beta, gamma)
        updated_ccoords = {}
        for n in sG.nodes():
            updated_ccoords[n] = fractional_to_cartesian(
                T_unitcell, sG.nodes[n]["fcoords"].T).T
            sG.nodes[n]["ccoords"] = updated_ccoords[n]
        return sG, updated_ccoords


class OptimizationDriver:
    """Driver for two-stage rotation optimization and cell-parameter optimization.

    Stage 1: minimize distance from rotated X positions to neighbor COMs.
    Stage 2: minimize pairwise distances between paired X atoms across edges.
    Cell optimization minimizes fractional-coordinate change when scaling the cell.
    """

    def __init__(self, comm=None, ostream=None):
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        # attributes to be set before use
        self.sorted_nodes = None
        self.sorted_edges = None
        self.pname_set_dict = None

        self.initial_rotations = None
        self.initial_set_rotations = None
        self.optimized_rotations = None
        self.optimized_set_rotations = None
        self.static_atom_positions = None

        # Optimization parameters
        self.opt_method = 'L-BFGS-B'
        self.maxfun = 15000
        self.maxiter = 15000
        self.display = True
        self.eps = 1e-8

        self.fixed_cell_shape = True

        self._debug = False

    def _objective_function_pre(self, params, G, static_atom_positions):
        """
        Objective function to minimize distances between paired node to paired node_com along edges.

        Parameters:
            params (numpy.ndarray): Flattened array of rotation matrices.
            G (networkx.Graph): Graph structure.
            atom_positions (dict): Original positions of X atoms for each node.


        Returns:
            float: Total distance metric to minimize.
        """
        # num_nodes = len(G.nodes())

        sorted_nodes = self.sorted_nodes
        sorted_edges = self.sorted_edges
        pname_set_dict = self.pname_set_dict
        set_rotation_matrices = params.reshape(len(pname_set_dict), 3, 3)
        rotation_matrices = expand_set_rots(pname_set_dict,
                                            set_rotation_matrices,
                                            sorted_nodes)
        total_distance = 0.0

        for i, j in sorted_edges:
            R_i = reorthogonalize_matrix(rotation_matrices[i])

            com_i = G.nodes[sorted_nodes[i]]["ccoords"]
            com_j = G.nodes[sorted_nodes[j]]["ccoords"]
            # Rotate positions around their mass center
            rotated_i_positions = (
                np.dot(static_atom_positions[i][:, 1:] - com_i, R_i.T) + com_i)

            #dist_matrix = np.empty((len(rotated_i_positions), 1))
            #for idx_i in range(len(rotated_i_positions)):
            #    dist = np.linalg.norm(rotated_i_positions[idx_i] - com_j)
            #    dist_matrix[idx_i, 0] = dist
            dist_matrix = np.linalg.norm(rotated_i_positions - com_j,
                                         axis=1,
                                         keepdims=True)
            # total_distance += dist ** 2
            if np.argmin(dist_matrix) > 1:
                total_distance += 1e4  # penalty for the distance difference
            total_distance += np.min(dist_matrix)**2

            total_distance += 1e3 / (np.max(dist_matrix) - np.min(dist_matrix)
                                     )  # reward for the distance difference

        return total_distance

    def _objective_function_after(self, params, G, static_atom_positions):
        """
        Objective function to minimize distances between paired atoms along edges. just use minimum distance

        Parameters:
            params (numpy.ndarray): Flattened array of rotation matrices.
            G (networkx.Graph): Graph structure.
            atom_positions (dict): Original positions of X atoms for each node.
            edge_pairings (dict): Precomputed pairings for each edge.

        Returns:
            float: Total distance metric to minimize.
        """
        # num_nodes = len(G.nodes())
        set_rotation_matrices = params.reshape(len(self.pname_set_dict), 3, 3)
        rotation_matrices = expand_set_rots(self.pname_set_dict,
                                            set_rotation_matrices,
                                            self.sorted_nodes)
        total_distance = 0.0

        #for i, j in self.sorted_edges:
        #    R_i = reorthogonalize_matrix(rotation_matrices[i])
        #    R_j = reorthogonalize_matrix(rotation_matrices[j])
        #
        #    com_i = G.nodes[self.sorted_nodes[i]]["ccoords"]
        #    com_j = G.nodes[self.sorted_nodes[j]]["ccoords"]

        # Rotate positions around their mass center
        #rotated_i_positions = (
        #    np.dot(static_atom_positions[i][:, 1:] - com_i, R_i.T) + com_i)
        #rotated_j_positions = (
        #    np.dot(static_atom_positions[j][:, 1:] - com_j, R_j.T) + com_j)
        #
        #dist_matrix = np.empty(
        #    (len(rotated_i_positions), len(rotated_j_positions)))
        #for idx_i in range(len(rotated_i_positions)):
        #    for idx_j in range(len(rotated_j_positions)):
        #        dist = np.linalg.norm(rotated_i_positions[idx_i] -
        #                              rotated_j_positions[idx_j])
        #        dist_matrix[idx_i, idx_j] = dist
        for i, j in self.sorted_edges:
            R_i = reorthogonalize_matrix(rotation_matrices[i])
            R_j = reorthogonalize_matrix(rotation_matrices[j])

            com_i = G.nodes[self.sorted_nodes[i]]["ccoords"]
            com_j = G.nodes[self.sorted_nodes[j]]["ccoords"]

            # Rotate positions around their mass center
            rotated_i_positions = (
                np.dot(static_atom_positions[i][:, 1:] - com_i, R_i.T) + com_i
            )  # shape (Ni, 3)
            rotated_j_positions = (
                np.dot(static_atom_positions[j][:, 1:] - com_j, R_j.T) + com_j
            )  # shape (Nj, 3)

            # Vectorized pairwise distance matrix
            diff = rotated_i_positions[:, None, :] - rotated_j_positions[
                None, :, :]
            # shape (Ni, Nj, 3)

            dist_matrix = np.linalg.norm(diff, axis=2)

            if np.argmin(dist_matrix) > 1:
                total_distance += 1e4  # penalty for the distance difference

            total_distance += np.min(dist_matrix)**2

        return total_distance

    def _optimize_rotations_pre(self, num_nodes, G, atom_positions,
                                initial_set_rotations):
        """
        Optimize rotations for all nodes in the graph.

        Parameters:
            G (networkx.Graph): Graph structure with edges between nodes.
            atom_positions (dict): Positions of X atoms for each node.

        Returns:
            list: Optimized rotation matrices for all nodes.
        """

        assert_msg_critical("scipy" in sys.modules,
                            "scipy is required for optimize_rotations_pre.")

        self.ostream.print_info(f"Rotations optimization information:")
        self.ostream.print_info(f"opt_method:, {self.opt_method}")
        self.ostream.print_info(f"maxfun:, {self.maxfun}")
        self.ostream.print_info(f"maxiter:, {self.maxiter}")
        self.ostream.print_info(f"display:, {self.display}")
        self.ostream.print_info(f"eps:, {self.eps}")
        self.ostream.print_info(f"Number of nodes to optimize:, {num_nodes}")
        self.ostream.print_info("\n")
        self.ostream.print_separator()
        self.ostream.print_info(f"Rotation Optimization (stage 1)")
        self.ostream.flush()

        # initial_rotations = np.tile(np.eye(3), (num_nodes, 1)).flatten()
        # get a better initial guess, use random rotation matrix combination
        # initial_rotations  = np.array([reorthogonalize_matrix(np.random.rand(3,3)) for i in range(num_nodes)]).flatten()
        static_atom_positions = atom_positions.copy()
        # Precompute edge-specific pairings
        # edge_pairings = find_edge_pairings(sorted_edges, atom_positions).

        result = minimize(
            self._objective_function_pre,
            initial_set_rotations.flatten(),
            args=(G, static_atom_positions),
            method=self.opt_method,
            options={
                "maxfun": self.maxfun,
                "maxiter": self.maxiter,
                "disp": self.display,
                "eps": self.eps,
                "maxls": 50,
            },
        )

        optimized_rotations = result.x

        return optimized_rotations, static_atom_positions

    def _optimize_rotations_after(self, num_nodes, G, atom_positions,
                                  initial_rotations):
        """
        Optimize rotations for all nodes in the graph.

        Parameters:
            G (networkx.Graph): Graph structure with edges between nodes.
            atom_positions (dict): Positions of X atoms for each node.

        Returns:
            list: Optimized rotation matrices for all nodes.
        """

        assert_msg_critical("scipy" in sys.modules,
                            "scipy is required for optimize_rotations_after.")
        self.ostream.print_info('-' * 20)
        self.ostream.print_info(f"Rotation Optimization (stage 2)")
        self.ostream.print_separator()
        self.ostream.flush()

        # get a better initial guess, use random rotation matrix combination
        # initial_rotations  = np.array([reorthogonalize_matrix(np.random.rand(3,3)) for i in range(num_nodes)]).flatten()
        static_atom_positions = atom_positions.copy()
        # Precompute edge-specific pairings
        # edge_pairings = find_edge_pairings(sorted_edges, atom_positions)

        result = minimize(
            self._objective_function_after,
            initial_rotations.flatten(),
            args=(G, static_atom_positions),
            method=self.opt_method,
            options={
                "maxfun": self.maxfun,
                "maxiter": self.maxiter,
                "disp": self.display,
                "eps": self.eps,
            },
        )

        optimized_rotations = result.x.reshape(-1, 3, 3)
        optimized_rotations = [
            reorthogonalize_matrix(R) for R in optimized_rotations
        ]
        optimized_rotations = np.array(optimized_rotations)

        return optimized_rotations, static_atom_positions

    def _scale_objective_function(self, params, old_cell_params,
                                  old_cartesian_coords, new_cartesian_coords,
                                  ratio_ba, ratio_ca):
        """Sum of squared differences between old fractional coords and new coords in new cell (optionally fixed shape)."""
        a_old, b_old, c_old, alpha_old, beta_old, gamma_old = old_cell_params
        a_new, b_new, c_new, _, _, _ = params
        #constrain the angles to be the same as old cell
        if self.fixed_cell_shape:
            b_new = a_new * ratio_ba
            c_new = a_new * ratio_ca

        # Compute transformation matrix for the old unit cell, T is the unit cell matrix
        T_old = unit_cell_to_cartesian_matrix(a_old, b_old, c_old, alpha_old,
                                              beta_old, gamma_old)
        T_old_inv = np.linalg.inv(T_old)
        old_fractional_coords = cartesian_to_fractional(
            old_cartesian_coords, T_old_inv)

        # backup
        # old_fractional_coords = cartesian_to_fractional(old_cartesian_coords,T_old_inv)

        # Compute transformation matrix for the new unit cell
        T_new = unit_cell_to_cartesian_matrix(a_new, b_new, c_new, alpha_old,
                                              beta_old, gamma_old)
        T_new_inv = np.linalg.inv(T_new)

        # Convert the new Cartesian coordinates to fractional coordinate using the old unit cell

        # Recalculate fractional coordinates from updated Cartesian coordinates
        new_fractional_coords = cartesian_to_fractional(
            new_cartesian_coords, T_new_inv)

        # Compute difference from original fractional coordinates
        diff = new_fractional_coords - old_fractional_coords
        return np.sum(diff**2)  # Sum of squared differences

    def _optimize_cell_params(self, cell_info, original_ccoords,
                              updated_ccoords):
        """Minimize fractional coordinate change when scaling cell to fit updated_ccoords; returns (a, b, c, alpha, beta, gamma)."""
        assert_msg_critical("scipy" in sys.modules,
                            "scipy is required for optimize_cell_parameters.")

        # Old cell parameters (example values)
        old_cell_params = cell_info  # [a, b, c, alpha, beta, gamma]

        # Old Cartesian coordinates of points (example values)
        old_cartesian_coords = np.vstack(list(
            original_ccoords.values()))  # original_ccoords

        # New Cartesian coordinates of the same points (example values)
        new_cartesian_coords = np.vstack(list(
            updated_ccoords.values()))  # updated_ccoords
        # Initial guess for new unit cell parameters (e.g., slightly modified cell)
        initial_params = cell_info

        # Bounds: a, b, c > 3; angles [0, 180]
        bounds = [(3, None), (3, None), (3, None)] + [(20, 180)] * 3

        ratio_ba = round(initial_params[1] / initial_params[0], 5)
        ratio_ca = round(initial_params[2] / initial_params[0], 5)

        # Optimize using L-BFGS-B to minimize the objective function
        result = minimize(
            self._scale_objective_function,
            x0=initial_params,
            args=(old_cell_params, old_cartesian_coords, new_cartesian_coords,
                  ratio_ba, ratio_ca),
            method="L-BFGS-B",
            bounds=bounds,
        )

        # Extract optimized parameters
        optimized_params = np.round(result.x, 5)
        self.ostream.print_info(
            f"Optimized New Cell Parameters: {optimized_params}\nTemplate Cell Parameters: {cell_info}"
        )
        if self.fixed_cell_shape:
            self.ostream.print_info(
                "Note: Cell shape is fixed during optimization.")
            optimized_params[1] = optimized_params[0] * (old_cell_params[1] /
                                                         old_cell_params[0])
            optimized_params[2] = optimized_params[0] * (old_cell_params[2] /
                                                         old_cell_params[0])
        return optimized_params

    def _update_ccoords_by_optimized_cell_params(self, G, optimized_params):
        sG = G.copy()
        a, b, c, alpha, beta, gamma = optimized_params
        T_unitcell = unit_cell_to_cartesian_matrix(a, b, c, alpha, beta, gamma)
        updated_ccoords = {}
        for n in sG.nodes():
            updated_ccoords[n] = fractional_to_cartesian(
                T_unitcell, sG.nodes[n]["fcoords"].T).T
            sG.nodes[n]["ccoords"] = updated_ccoords[n]
        return sG, updated_ccoords


def recenter_and_norm_vectors(vectors, extra_mass_center=None):
    """Center vectors (optionally at extra_mass_center) and normalize each row. Returns (normalized_vectors, mass_center)."""
    vectors = np.array(vectors)
    if extra_mass_center is not None:
        mass_center = extra_mass_center
    else:
        mass_center = np.mean(vectors, axis=0)
    vectors = vectors - mass_center
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    return vectors, mass_center


def get_connected_nodes_vectors(node, G):
    """Return list of neighbor ccoords and this node's ccoords from G."""
    vectors = []
    for i in list(G.neighbors(node)):
        vectors.append(G.nodes[i]["ccoords"])
    return vectors, G.nodes[node]["ccoords"]


def get_rot_trans_matrix(node, G, sorted_nodes, Xatoms_positions_dict):
    """Compute rotation and translation to align node X vectors to neighbor directions (for initial guess)."""
    node_id = sorted_nodes.index(node)
    node_xvecs = Xatoms_positions_dict[node_id][:, 1:]
    vecsA, _ = recenter_and_norm_vectors(node_xvecs, extra_mass_center=None)
    v2, node_center = get_connected_nodes_vectors(node, G)
    vecsB, _ = recenter_and_norm_vectors(v2, extra_mass_center=node_center)
    _, rot, tran = superimpose_rotation_only(vecsA, vecsB)
    return rot, tran


def expand_set_rots(pname_set_dict, set_rotations, sorted_nodes):
    """Expand one rotation per pname to a full list of rotations per sorted_nodes index."""
    set_rotations = set_rotations.reshape(len(pname_set_dict), 3, 3)
    rotations = np.empty((len(sorted_nodes), 3, 3))
    idx = 0
    for name in pname_set_dict:
        for k in pname_set_dict[name]["ind_ofsortednodes"]:
            rotations[k] = set_rotations[idx]
        idx += 1
    return rotations

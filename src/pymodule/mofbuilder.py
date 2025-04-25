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

from pathlib import Path
import numpy as np
import time
from copy import deepcopy
import re

from .mofpreparer import MofPreparer
from .mofutils import (NetOptimizer, extract_node_name_from_gro_resindex,
                       replace_edges_by_callname, make_dummy_split_node_dict,
                       rename_node_arr, merge_metal_list_to_node_array,
                       save_node_edge_term_gro, update_matched_nodes_xind,
                       gro_string_show)
from .mofmdprepare import (
    get_residues_forcefield,
    get_itps,
    genrate_top_file,
    copy_mdps,
)


class MofBuilder:

    def __init__(self):
        # call preparation driver
        # find database path which should be decided later
        # load the MOF_topology_dict flie in database folder
        self.preparation = MofPreparer()
        self.mof_family = None
        self.node_metal = None
        self.linker_xyz_file = None
        self.supercell = (1, 1, 1)
        self.dummy_node = False
        self.bridge_node = False
        self.bridge_node_search_range = 0.3

    def show_available_mof_families(self):
        self.preparation.list_mof_family()
        # print(
        #    "MOF builder is initialized, please prepare the building material by calling the preparation driver"
        # )
        # print("***check the preparation status by calling preparation_check()***")

    def preparation_check(self):
        # check the preparation status
        # necessary for MofBuilder to start building
        preparation = self.preparation
        if preparation.check_status():
            print("-" * 80)
            print(" " * 20, "Preparation is completed", " " * 20)
            print("-" * 80)

            print("MOF builder is ready to build")
            self.mof_family = preparation.mof_family
            self.template_cif = preparation.selected_template_cif_file
            self.node_pdb = preparation.selected_node_pdb_file
            self.node_metal = preparation.node_metal
            self.linker_pdb = preparation.selected_linker_edge_pdb
            self.linker_center_pdb = (preparation.selected_linker_center_pdb
                                      )  # could be None if ditopic linker
            self.linker_topic = preparation.linker_topic
            self.linker_xyz = preparation.linker_xyz
            return True
        else:
            print("Error: Could not find the required files")
            print("Please redo the preparation steps")
            return False

    def set_supercell(self, supercell):
        self.supercell = supercell

    # maxfun, maxiter, display,eps, iprint,method
    def set_rotation_optimizer_maxfun(self, maxfun):
        self.rotation_optimizer_maxfun = maxfun

    def set_rotation_optimizer_maxiter(self, maxiter):
        self.rotation_optimizer_maxiter = maxiter

    def set_rotation_optimizer_display(self, display):
        self.rotation_optimizer_display = display

    def set_rotation_optimizer_eps(self, eps):
        self.rotation_optimizer_eps = eps

    def set_rotation_optimizer_iprint(self, iprint):
        self.rotation_optimizer_iprint = iprint

    def set_rotation_optimizer_method(self, method):
        self.rotation_optimizer_method = method

    def set_node_topic(self, node_topic):
        self.node_topic = node_topic

    def set_linker_topic(self, linker_topic):
        self.linker_topic = linker_topic

    def set_connection_constant_length(self, connection_constant_length):
        self.connection_constant_length = connection_constant_length

    def set_node_termination(self, node_termination):
        self.node_termination = node_termination

    def set_vitual_edge_search(self, vitual_edge_search):
        self.vitual_edge = True

    def save_optimized_rotations(self, filename):
        self.optimized_rotations_filename = filename

    def set_use_saved_optimized_rotations_npy(self, saved_optimized_rotations):
        saved_optimized_rotations = saved_optimized_rotations + ".npy"

        if Path(saved_optimized_rotations).exists():
            self.saved_optimized_rotations = np.load(saved_optimized_rotations,
                                                     allow_pickle=True)
            print("Optimized rotations are loaded from: ",
                  saved_optimized_rotations)
        else:
            print(
                f"Could not find the saved optimized rotations:  {saved_optimized_rotations} will start the optimization from the beginning"
            )
            pass

    def set_use_saved_rotations_as_initial_guess(
            self, use_saved_rotations_as_initial_guess):
        """
        use the saved optimized rotations as initial guess
        """
        if hasattr(self, "saved_optimized_rotations"):
            self.use_saved_rotations_as_initial_guess = (
                use_saved_rotations_as_initial_guess)
        else:
            print(
                "saved_optimized_rotations is not found, will start the optimization from the beginning"
            )
            pass

    def set_supercell_cleaved_buffer_plus(self, buffer_plus_ratio):
        self.supercell_cleaved_buffer_plus = buffer_plus_ratio

    def set_supercell_cleaved_buffer_minus(self, buffer_minus_ratio):
        self.supercell_cleaved_buffer_minus = buffer_minus_ratio

    def build(self):
        self.preparation.select_mof_family(self.mof_family)
        self.preparation.select_node_metal(self.node_metal)
        if hasattr(self, "dummy_node"):
            self.preparation.use_dummy_node(self.dummy_node)
        self.preparation.fetch_node()
        self.preparation.fetch_linker(self.linker_xyz_file)
        if not self.preparation_check():
            print("Error: Could not find the required files")
            print("Please redo the preparation steps")
            return

        # check before building
        if not hasattr(self, "supercell"):
            self.supercell = (1, 1, 1)
        self.supercell = list(
            [self.supercell[0], self.supercell[1], self.supercell[2]])

        if self.linker_topic == 2:
            print("ditopic mof builder driver is called")
            start_time = time.time()
            linker_pdb = self.linker_pdb
            template_cif = self.template_cif
            node_pdb = self.node_pdb
            supercell = self.supercell
            self.net = NetOptimizer()
            if hasattr(self, "connection_constant_length"):
                self.net.set_constant_length(self.connection_constant_length)
            if hasattr(self, "rotation_optimizer_maxfun"):
                self.net.set_rotation_optimizer_maxfun(
                    self.rotation_optimizer_maxfun)
            if hasattr(self, "rotation_optimizer_maxiter"):
                self.net.set_rotation_optimizer_maxiter(
                    self.rotation_optimizer_maxiter)
            if hasattr(self, "rotation_optimizer_method"):
                self.net.set_rotation_optimizer_method(
                    self.rotation_optimizer_method)
            if hasattr(self, "rotation_optimizer_eps"):
                self.net.set_rotation_optimizer_eps(
                    self.rotation_optimizer_eps)
            if hasattr(self, "rotation_optimizer_iprint"):
                self.net.set_rotation_optimizer_iprint(
                    self.rotation_optimizer_iprint)
            if hasattr(self, "rotation_optimizer_display"):
                self.net.set_rotation_optimizer_display(
                    self.rotation_optimizer_display)

            if hasattr(self, "saved_optimized_rotations"):
                self.net.load_saved_optimized_rotations(
                    self.saved_optimized_rotations)
            if hasattr(self, "optimized_rotations_filename"):
                self.net.to_save_optimized_rotations(
                    self.optimized_rotations_filename)
            if hasattr(self, "use_saved_rotations_as_initial_guess"):
                self.net.use_saved_rotations_as_initial_guess(
                    self.use_saved_rotations_as_initial_guess)

            self.net.analyze_template_ditopic(template_cif)
            self.net.node_info(node_pdb)
            self.net.linker_info(linker_pdb)
            self.net.optimize()
            print("-" * 80)
            print(
                " " * 15,
                "Building time cost: %.5f seconds " %
                (time.time() - start_time),
            )
            print("-" * 80)
            if hasattr(self, "supercell_cleaved_buffer_plus"):
                cleaved_buffer_plus = self.supercell_cleaved_buffer_plus
            else:
                cleaved_buffer_plus = 0.0
            if hasattr(self, "supercell_cleaved_buffer_minus"):
                cleaved_buffer_minus = self.supercell_cleaved_buffer_minus
            else:
                cleaved_buffer_minus = 0.0

            self.net.set_supercell(supercell)
            self.net.place_edge_in_net()
            self.net.make_supercell_ditopic()
            """bridge node"""
            if self.bridge_node:
                self.net.set_virtual_edge(
                    self.bridge_node,
                    self.bridge_node_search_range,
                    self.bridge_node_max_neighbor,
                )

                self.net.superG = self.net.add_virtual_edge_for_bridge_node(
                    self.net.superG)
            self.net.make_eG_from_supereG_ditopic()
            self.net.main_frag_eG()
            self.archive_eG = self.net.eG.copy()
            self.net.add_xoo_to_edge_ditopic()
            self.net.make_supercell_range_cleaved_eG(
                buffer_plus=cleaved_buffer_plus,
                buffer_minus=cleaved_buffer_minus)
            self.net.find_unsaturated_node_eG()
            # self.net.add_xoo_to_edge_ditopic()

            if hasattr(self, "node_termination"):
                self.net.set_node_terminamtion(self.node_termination)
            # default termination is methyl in data folder
            self.net.set_node_terminamtion(
                str(
                    Path(
                        self.preparation.data_path,
                        "terminations_database",
                        "methyl.pdb",
                    )))
            ##TODO:
            # update  self.net.unsaturated_node  and self.net.matched_vnode_xind
            self.net.add_terminations_to_unsaturated_node()
            self.net.remove_xoo_from_node()
            # self.net = net

        elif self.linker_topic > 2:
            print("multitopic mof builder driver is called")
            start_time = time.time()
            linker_pdb = self.linker_pdb
            linker_center_pdb = self.linker_center_pdb
            template_cif = self.template_cif
            node_pdb = self.node_pdb
            supercell = self.supercell
            self.net = NetOptimizer()
            if hasattr(self, "connection_constant_length"):
                self.net.set_constant_length(self.connection_constant_length)

            if hasattr(self, "rotation_optimizer_maxfun"):
                self.net.set_rotation_optimizer_maxfun(
                    self.rotation_optimizer_maxfun)
            if hasattr(self, "rotation_optimizer_maxiter"):
                self.net.set_rotation_optimizer_maxiter(
                    self.rotation_optimizer_maxiter)
            if hasattr(self, "rotation_optimizer_method"):
                self.net.set_rotation_optimizer_method(
                    self.rotation_optimizer_method)

            if hasattr(self, "rotation_optimizer_eps"):
                self.net.set_rotation_optimizer_eps(
                    self.rotation_optimizer_eps)

            if hasattr(self, "rotation_optimizer_iprint"):
                self.net.set_rotation_optimizer_iprint(
                    self.rotation_optimizer_iprint)
            if hasattr(self, "rotation_optimizer_display"):
                self.net.set_rotation_optimizer_display(
                    self.rotation_optimizer_display)

            if hasattr(self, "saved_optimized_rotations"):
                self.net.load_saved_optimized_rotations(
                    self.saved_optimized_rotations)
            if hasattr(self, "optimized_rotations_filename"):
                self.net.to_save_optimized_rotations(
                    self.optimized_rotations_filename)
            if hasattr(self, "use_saved_rotations_as_initial_guess"):
                self.net.use_saved_rotations_as_initial_guess(
                    self.use_saved_rotations_as_initial_guess)

            self.net.analyze_template_multitopic(template_cif)
            self.net.node_info(node_pdb)
            self.net.linker_info(linker_pdb)
            self.net.linker_center_info(linker_center_pdb)
            self.net.optimize()
            print("-" * 80)
            print(
                " " * 15,
                "Building time cost: %.5f seconds " %
                (time.time() - start_time),
            )
            print("-" * 80)
            if hasattr(self, "supercell_cleaved_buffer_plus"):
                cleaved_buffer_plus = self.supercell_cleaved_buffer_plus
            else:
                cleaved_buffer_plus = 0.0
            if hasattr(self, "supercell_cleaved_buffer_minus"):
                cleaved_buffer_minus = self.supercell_cleaved_buffer_minus
            else:
                cleaved_buffer_minus = 0.0

            self.net.set_supercell(supercell)
            self.net.place_edge_in_net()
            self.net.make_supercell_multitopic()
            self.net.make_eG_from_supereG_multitopic()
            self.net.main_frag_eG()
            self.archive_eG = self.net.eG.copy()
            self.net.add_xoo_to_edge_multitopic()
            self.net.make_supercell_range_cleaved_eG(
                buffer_plus=cleaved_buffer_plus,
                buffer_minus=cleaved_buffer_minus)

            self.net.find_unsaturated_node_eG()
            if hasattr(self, "node_termination"):
                self.net.set_node_terminamtion(self.node_termination)
            # default termination is methyl in data folder
            self.net.set_node_terminamtion(
                str(
                    Path(
                        self.preparation.data_path,
                        "terminations_database",
                        "methyl.pdb",
                    )))
            self.net.add_terminations_to_unsaturated_node()
            self.net.remove_xoo_from_node()
        self.mofG = self.net.eG.copy()
        self.net.extract_node_edge_term()
        self.nodes_eG = self.net.nodes_eG.copy()
        self.edges_eG = self.net.edges_eG.copy()
        self.terms_eG = self.net.terms_eG.copy()

    def get_gro_lines_list(self, graphG):
        merged_node_edge_term = self.net.get_node_edge_term_grolines(
            graphG, self.net.sc_unit_cell)
        head = []
        head.append("MOF\n")
        head.append(str(len(merged_node_edge_term)) + "\n")
        tail = ["20 20 20 \n"]
        return head + merged_node_edge_term + tail

    def write_gromacs_files(self, gro_name=None):
        if gro_name is not None:
            self.gro_name = gro_name
        else:
            self.gro_name = ("mof_" + str(self.mof_family.split(".")[0]) +
                             "_" +
                             self.linker_xyz.strip(".xyz").split("/")[-1] +
                             ".gro")
            print("gro_name is not set, will be saved as: ", self.gro_name)
        grolines = self.get_gro_lines_list(self.mofG)
        Path("output_gros").mkdir(parents=True, exist_ok=True)
        gro_file_path = str(Path("output_gros", self.gro_name))
        print("writing the gromacs file", gro_file_path)
        with open(gro_file_path, "w") as f:
            f.writelines(grolines)
        if hasattr(self, "defective_mofG"):
            defective_gro_path = str(
                Path("output_gros", "defective_" + self.gro_name))
            print("writing the gromacs file", defective_gro_path)
            defective_gro_lines = self.get_gro_lines_list(self.defective_mofG)
            with open(defective_gro_path, "w") as f:
                f.writelines(defective_gro_lines)

    def show(self, width=800, height=600, res_indices=False, res_names=False):
        if hasattr(self, "defective_mofG"):
            grolines = self.get_gro_lines_list(self.defective_mofG)
        else:
            grolines = self.get_gro_lines_list(self.mofG)
        gro_string_show(
            gro_lines_list=grolines,
            w=width,
            h=height,
            res_id=res_indices,
            res_name=res_names,
        )

    # functions for defects are under construction
    def remove(self,
               nodes=[],
               linkers=[],
               update_node_termination=False,
               clean_unsaturated_linkers=False):
        if not hasattr(self, "defective_net"):
            self.defective_net = deepcopy(self.net)

        print(
            "built MOF is saved",
            #    "nodes: ",
            #    len(self.saved_eG.nodes),
            #    "edges: ",
            #    len(self.saved_eG.edges),
        )

        #self.defective_net.eG = self.archive_eG.copy()
        remove_node_list = nodes
        remove_edge_list = linkers

        remove_edge_list = [
            str(int(i) - len(self.nodes_eG)) for i in remove_edge_list
        ]  # TODO: check if it is correct

        self.to_remove_nodes_name = extract_node_name_from_gro_resindex(
            remove_node_list, self.net.nodes_eG)
        self.to_remove_edges_name = extract_node_name_from_gro_resindex(
            remove_edge_list, self.net.edges_eG)

        if hasattr(self, "supercell_cleaved_buffer_plus"):
            cleaved_buffer_plus = self.supercell_cleaved_buffer_plus
        else:
            cleaved_buffer_plus = 0.0
        if hasattr(self, "supercell_cleaved_buffer_minus"):
            cleaved_buffer_minus = self.supercell_cleaved_buffer_minus
        else:
            cleaved_buffer_minus = 0.0

        if self.linker_topic == 2:
            self.defective_net.add_xoo_to_edge_ditopic()
        elif self.linker_topic > 2:
            self.defective_net.add_xoo_to_edge_multitopic()

        if clean_unsaturated_linkers:
            self.to_remove_edges_name.update(self.saved_eG_unsaturated_linker)

        for node_name in self.to_remove_nodes_name:
            if node_name in self.defective_net.eG.nodes():
                self.defective_net.eG.remove_node(node_name)
            else:
                print("node ", node_name, " is not in this MOF")
        for edge_name in self.to_remove_edges_name:
            neighbors = list(self.defective_net.eG.neighbors(edge_name))
            if len(neighbors) == 2:  # ditopic linker case
                self.defective_net.eG.remove_edge(neighbors[0], neighbors[1])
            if edge_name in self.defective_net.eG.nodes():
                self.defective_net.eG.remove_node(edge_name)
            else:
                print("edge ", edge_name, " is not in this MOF")

        #print(self.to_remove_edges_name, " will be removed edge")
        #print(self.to_remove_nodes_name, " will be removed node")
        # update the matched_vnode_xind
        self.defective_net.matched_vnode_xind = update_matched_nodes_xind(
            self.to_remove_nodes_name,
            self.to_remove_edges_name,
            self.defective_net.matched_vnode_xind,
        )
        # sort subgraph by connectivity

        print(
            "defective MOF is updated",
            #    "nodes: ",
            #    len(self.defective_net.eG.nodes),
            #    "edges: ",
            #    len(self.defective_net.eG.edges),
        )

        self.defective_net.make_supercell_range_cleaved_eG(
            buffer_plus=cleaved_buffer_plus, buffer_minus=cleaved_buffer_minus)

        if update_node_termination:
            self.defective_net.find_unsaturated_node_eG()
        else:
            self.defective_net.unsaturated_node = self.saved_eG_unsaturated_node
            self.defective_net.matched_vnode_xind = self.saved_eG_matched_vnode_xind

        self.defective_net.add_terminations_to_unsaturated_node()
        self.defective_net.remove_xoo_from_node()

        self.defective_mofG = self.defective_net.eG.copy()
        self.defective_net.extract_node_edge_term()
        self.defective_mof_nodes_eG = self.defective_net.nodes_eG.copy()
        self.defective_mof_edges_eG = self.defective_net.edges_eG.copy()
        self.defective_mof_terms_eG = self.defective_net.terms_eG.copy()

    def exchange_linkers(self, linkers=[], exchange_linker_pdb=None):
        if exchange_linker_pdb is None:
            raise ValueError(
                "exchange_linker_pdb is not set, please set it before calling exchange_linkers()"
            )
        if hasattr(self, "defective_net"):
            defective_net = deepcopy(self.defective_net)
            nodes_eG = self.defective_net.nodes_eG
            edges_eG = self.defective_net.edges_eG
        else:
            defective_net = deepcopy(self.net)
            nodes_eG = self.nodes_eG
            edges_eG = self.edges_eG

        exchange_edge_list = linkers

        exchange_edge_list = [
            str(int(i) - len(nodes_eG)) for i in exchange_edge_list
        ]  # TODO: this will be updated every time if a defect is made before
        exchange_edges_name = extract_node_name_from_gro_resindex(
            exchange_edge_list, edges_eG)

        defective_net.eG = replace_edges_by_callname(
            exchange_edges_name,
            defective_net.eG,
            defective_net.sc_unit_cell_inv,
            exchange_linker_pdb,
            prefix="R",
        )

        self.defective_net = defective_net
        self.defective_mofG = self.defective_net.eG.copy()

    def write_defective_split_node_gro_again(self, gro_name):
        if not self.preparation.dummy_node:
            print("dummy node is not used, splitting node is not possible")
            return
        print(
            "splitting node and saving gro again, called after write_defect_gro()"
        )
        nodes_eG = self.defective_mof_nodes_eG
        edges_eG = self.defective_mof_edges_eG
        terms_eG = self.defective_mof_terms_eG
        node_split_dict = make_dummy_split_node_dict(self.node_pdb)
        nodes_eGarr = np.vstack(nodes_eG)
        metals_list, hho_list, ho_list, o_list = rename_node_arr(
            node_split_dict, nodes_eGarr)

        merged_split_node_edge_term = []
        line_num = 0
        res_count = 0
        print("writing split node gro")
        for splitted_node in [
                metals_list,
                hho_list,
                ho_list,
                o_list,
                edges_eG,
                terms_eG,
        ]:
            merged_split_node_edge_term, line_num, res_count = (
                merge_metal_list_to_node_array(merged_split_node_edge_term,
                                               splitted_node, line_num,
                                               res_count))

        print("metal_res_num: ", len(metals_list))
        print("hho_res_num: ", len(hho_list))
        print("ho_res_num: ", len(ho_list))
        print("o_res_num: ", len(o_list))
        print("edge_res_num: ", len(edges_eG))
        print("term_res_num: ", len(terms_eG))

        save_node_edge_term_gro(merged_split_node_edge_term, gro_name)

    def write_split_node_gro_again(self, gro_name):
        if not self.preparation.dummy_node:
            print("dummy node is not used, splitting node is not possible")
            return

        print("splitting node and saving gro again, called after write_gro()")

        nodes_eG = self.nodes_eG
        edges_eG = self.edges_eG
        terms_eG = self.terms_eG

        node_split_dict = make_dummy_split_node_dict(self.node_pdb)
        nodes_eGarr = np.vstack(nodes_eG)
        metals_list, hho_list, ho_list, o_list = rename_node_arr(
            node_split_dict, nodes_eGarr)

        merged_split_node_edge_term = []
        line_num = 0
        res_count = 0
        print("writing split node gro")
        for splitted_node in [
                metals_list,
                hho_list,
                ho_list,
                o_list,
                edges_eG,
                terms_eG,
        ]:
            merged_split_node_edge_term, line_num, res_count = (
                merge_metal_list_to_node_array(merged_split_node_edge_term,
                                               splitted_node, line_num,
                                               res_count))

        save_node_edge_term_gro(merged_split_node_edge_term, gro_name)

        print("metal_res_num: ", len(metals_list))
        print("hho_res_num: ", len(hho_list))
        print("ho_res_num: ", len(ho_list))
        print("o_res_num: ", len(o_list))
        print("edge_res_num: ", len(edges_eG))
        print("term_res_num: ", len(terms_eG))

    def md_prepare(self):
        if not hasattr(self, "linker_file_ff"):
            self.linker_file_ff = None

        def make_array(lists):
            arr = np.zeros((len(lists), 7), dtype="object")
            for idx, line in enumerate(lists):
                line = line.strip("/n")
                arr[idx][0] = int(line[0:5])
                arr[idx][1] = re.sub(r" ", "", str(line[5:10]))  # remove space
                arr[idx][2] = re.sub(r" ", "",
                                     str(line[10:15]))  # remove space
                arr[idx][3] = int(line[15:20])
                arr[idx][4] = float(line[20:28]) * 10
                arr[idx][5] = float(line[28:36]) * 10
                arr[idx][6] = float(line[36:44]) * 10
            return arr

        arr, node_split_dict = (
            make_array(self.merged_split_node_edge_term),
            self.node_split_dict,
        )
        new_arr, self.res_info, self.restypes = get_residues_forcefield(
            arr,
            node_split_dict,
            self.dummy_node,
            self.linker_xyz_file,
            self.linker_file_ff,
            self.linker_topic,
            self.preparation.linker_center_frag_nodes_num,
            self.preparation.linker_center_Xs,
            self.preparation.linker_single_frag_nodes_num,
            self.preparation.linker_frag_Xs,
        )

        self.itp_dir = get_itps(
            self.preparation.data_path,
            self.restypes,
            self.node_metal,
            self.net.node_termination,
            sol_list=[],
        )
        self.top_path = genrate_top_file(self.itp_dir,
                                         self.preparation.data_path,
                                         self.res_info)
        self.mdp_dir = copy_mdps(self.preparation.data_path)

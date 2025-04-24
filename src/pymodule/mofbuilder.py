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

from .mofpreparer import MofPreparer
from .mofutils import (
    NetOptimizer,
    extract_node_name_from_gro_resindex,
    replace_edges_by_callname,
    make_dummy_split_node_dict,
    rename_node_arr,
    merge_metal_list_to_node_array,
    save_node_edge_term_gro,
    find_unsaturated_linker,
    update_matched_nodes_xind,
    gro_show,
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

    def remove_nodes(self, remove_node_list):
        self.remove_node_list = remove_node_list

    def remove_linkers(self, remove_edge_list):
        self.remove_edge_list = remove_edge_list

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
                self.net.set_rotation_optimizer_eps(self.rotation_optimizer_eps)
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
                self.net.set_rotation_optimizer_eps(self.rotation_optimizer_eps)

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

    def write_gromacs_files(self, gro_name=None):
        if hasattr(self, "saved_eG"):
            if self.supercell == self.saved_supercell:
                print("saved_eG is found, will write the preserved eG")
                self.net.eG = self.saved_eG
                self.net.write_node_edge_node_gro(self.gro_name)
                return

        if gro_name is not None:
            self.gro_name = gro_name
        else:
            self.gro_name = ("mof_" + str(self.mof_family.split(".")[0]) + "_" +
                             self.linker_xyz.strip(".xyz").split("/")[-1])
            print("gro_name is not set, will be saved as: ",
                  self.gro_name + ".gro")

        print("writing gro file")
        print("nodes:", len(self.net.eG.nodes()), "edges:",
              len(self.net.eG.edges()))
        self.net.write_node_edge_node_gro(self.gro_name)

    def show(self, width=800, height=600, res_indices=False, res_names=False):
        gro_file_path = str(Path("output_gros", self.gro_name))
        gro_show(
            gro_file_path + ".gro",
            w=width,
            h=height,
            res_id=res_indices,
            res_name=res_names,
        )

    # functions are under construction
    def make_defects_missing(self,
                             update_node_term=False,
                             clean_unsaturated_linkers=False):
        self.saved_eG = self.net.eG.copy(
        )  # save the original eG before making defects
        self.saved_eG_unsaturated_node = self.net.unsaturated_node
        self.saved_eG_matched_vnode_xind = self.net.matched_vnode_xind
        self.saved_eG_unsaturated_linker = find_unsaturated_linker(
            self.net.eG, self.linker_topic)
        self.saved_supercell = self.supercell
        # herit the original net to defective net
        self.defective_net = self.net

        print(
            "saved_eG is saved",
            "nodes: ",
            len(self.saved_eG.nodes),
            "edges: ",
            len(self.saved_eG.edges),
        )

        self.defective_net.eG = self.archive_eG.copy()
        remove_node_list = []
        remove_edge_list = []
        if hasattr(self, "remove_node_list"):
            remove_node_list = self.remove_node_list
        if hasattr(self, "remove_edge_list"):
            remove_edge_list = self.remove_edge_list
            remove_edge_list = [
                str(int(i) - len(self.net.nodes_eG)) for i in remove_edge_list
            ]  # TODO: check if it is correct

        to_remove_nodes_name = extract_node_name_from_gro_resindex(
            remove_node_list, self.net.nodes_eG)
        to_remove_edges_name = extract_node_name_from_gro_resindex(
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
            to_remove_edges_name.update(self.saved_eG_unsaturated_linker)

        for node_name in to_remove_nodes_name:
            self.defective_net.eG.remove_node(node_name)
        for edge_name in to_remove_edges_name:
            neighbors = list(self.defective_net.eG.neighbors(edge_name))
            if len(neighbors) == 2:  # ditopic linker case
                self.defective_net.eG.remove_edge(neighbors[0], neighbors[1])
            self.defective_net.eG.remove_node(edge_name)
        self.defective_net.main_frag_eG()
        # update the matched_vnode_xind
        self.defective_net.matched_vnode_xind = update_matched_nodes_xind(
            to_remove_nodes_name,
            to_remove_edges_name,
            self.defective_net.matched_vnode_xind,
        )
        # sort subgraph by connectivity

        print(
            "defective eG is updated",
            "nodes: ",
            len(self.defective_net.eG.nodes),
            "edges: ",
            len(self.defective_net.eG.edges),
        )
        self.defective_net.make_supercell_range_cleaved_eG(
            buffer_plus=cleaved_buffer_plus, buffer_minus=cleaved_buffer_minus)

        if update_node_term:
            self.defective_net.find_unsaturated_node_eG()
        else:
            self.defective_net.unsaturated_node = self.saved_eG_unsaturated_node
            self.defective_net.matched_vnode_xind = self.saved_eG_matched_vnode_xind

        self.defective_net.add_terminations_to_unsaturated_node()
        self.defective_net.remove_xoo_from_node()

    # def set_exchange_node_list(self, exchange_node_list): #avoid node exchange
    #    self.exchange_node_list = exchange_node_list
    def set_exchange_edge_list(self, exchange_edge_list):
        self.exchange_edge_list = exchange_edge_list

    # def set_to_exchange_node_pdb(self, to_exchange_node_pdb): #avoid node exchange
    #    self.to_exchange_node_pdb = to_exchange_node_pdb
    def set_to_exchange_edge_pdb(self, to_exchange_edge_pdb):
        self.to_exchange_edge_pdb = to_exchange_edge_pdb

    def make_defects_exchange(self):
        defective_net = self.net
        # if hasattr(self, 'exchange_node_list'):
        #    exchange_node_list = self.exchange_node_list
        #    exchange_nodes_name = extract_node_name_from_gro_resindex(exchange_node_list, self.net.nodes_eG)
        if hasattr(self, "exchange_edge_list"):
            exchange_edge_list = self.exchange_edge_list
            exchange_edge_list = [
                str(int(i) - len(self.net.nodes_eG)) for i in exchange_edge_list
            ]  # TODO: check if it is correct
            exchange_edges_name = extract_node_name_from_gro_resindex(
                exchange_edge_list, self.net.edges_eG)

        # if hasattr(self, 'to_exchange_node_pdb'):
        #    to_exchange_node_pdb = self.to_exchange_node_pdb
        if hasattr(self, "to_exchange_edge_pdb"):
            to_exchange_edge_pdb = self.to_exchange_edge_pdb
        # TODO:

        if hasattr(self, "exchange_edge_list") and hasattr(
                self, "to_exchange_edge_pdb"):
            print(
                "exchange_edge_list and to_exchange_edge_pdb are set, will exchange the edges"
            )
            defective_net.eG = replace_edges_by_callname(
                exchange_edges_name,
                defective_net.eG,
                defective_net.sc_unit_cell_inv,
                to_exchange_edge_pdb,
                prefix="R",
            )
        # if hasattr(self, 'exchange_node_list') and hasattr(self, 'to_exchange_node_pdb'):
        #    print('exchange_node_list and to_exchange_node_pdb are set, will exchange the nodes')
        #    defective_net.eG = replace_edges_by_callname (exchange_nodes_name,defective_net.eG,defective_net.sc_unit_cell_inv,to_exchange_node_pdb, prefix='R')
        self.defective_net = defective_net

    def set_defect_gro_name(self, defect_gro_name):
        self.defect_gro_name = defect_gro_name

    def write_defective_model_gromacs_file(self):
        # if not hasattr(self, 'defective_net'):
        #    print('defective_net is not set')
        #    print('make_defects_missing() or make_defects_exchange() should be called before write_defect_gro(), or you can write with write_gro()')
        #    return

        if not hasattr(self, "defect_gro_name"):
            self.defect_gro_name = (
                "defective_mof_" + str(self.mof_family.split(".")[0]) + "_" +
                self.linker_xyz.strip(".xyz").split("/")[-1])
            print(
                "defect_gro_name is not set, will be saved as: ",
                self.defect_gro_name + ".gro",
            )
        print("writing defective gro file")

        self.defective_net.write_node_edge_node_gro(self.defect_gro_name)

    def show_defective_model(self,
                             width=800,
                             height=600,
                             res_indices=False,
                             res_names=False):
        gro_file_path = str(Path("output_gros", self.defect_gro_name))
        print("showing the gromacs file", gro_file_path)
        gro_show(
            gro_file_path + ".gro",
            w=width,
            h=height,
            res_id=res_indices,
            res_name=res_names,
        )

    def write_defective_split_node_gro_again(self, gro_name):
        if not self.preparation.dummy_node:
            print("dummy node is not used, splitting node is not possible")
            return
        print(
            "splitting node and saving gro again, called after write_defect_gro()"
        )
        nodes_eG = self.defective_net.nodes_eG
        edges_eG = self.defective_net.edges_eG
        terms_eG = self.defective_net.terms_eG
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

        nodes_eG = self.net.nodes_eG
        edges_eG = self.net.edges_eG
        terms_eG = self.net.terms_eG

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


if __name__ == "__main__":
    ##mof = MofBuilder()
    ##mof.mof_family = "UiO-67"
    ##mof.node_metal = "Zr"
    ##
    ##mof.linker_xyz_file = "database/linker4test/bdc.xyz"
    ##
    ##mof.set_rotation_optimizer_display(False)
    ##
    ### mof.use_saved_optimized_rotations_npy('rota')
    ### save optimized rotations to numpy file for later use
    ##mof.save_optimized_rotations("rota")
    ###mof.set_supercell_cleaved_buffer_minus(0.2)
    ##mof.supercell = (1, 1, 1)
    ##mof.build()
    ##mof.write_gromacs_files()
    ##mof.show(res_indices=True, res_names=False)
    ##
    ##mof.make_defects_missing(clean_unsaturated_linkers=True, update_node_term=True)
    ##mof.write_defective_model_gromacs_file()
    ##mof.show_defective_model(res_indices=True, res_names=False)

    mof = MofBuilder()
    mof.mof_family = "PCN-222"
    mof.node_metal = "Zr"
    mof.linker_xyz_file = "../database/linker4test/tcpp.xyz"
    mof.save_optimized_rotations("rotc")

    mof.supercell = (1, 1, 1)
    mof.build()
    mof.write_gromacs_files()
    mof.show(res_indices=True, res_names=False)
    mof.remove_linkers = [11]
    mof.make_defects_missing(clean_unsaturated_linkers=True,
                             update_node_term=True)
    mof.write_defective_model_gromacs_file()
    mof.show_defective_model(res_indices=True, res_names=False)

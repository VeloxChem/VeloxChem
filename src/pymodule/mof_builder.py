import time
import numpy as np
import os
from mofbuilder_prepare import prepare
from mofbuilder_utils import (
    net_optimizer,
    nn,
    extract_node_name_from_gro_resindex,
    replace_edges_by_callname,
    make_dummy_split_node_dict,
    rename_node_arr,
    merge_metal_list_to_node_array,
    save_node_edge_term_gro,
)


class MofBuilder:
    def __init__(self):
        # call preparation driver
        # find database path which should be decided later
        # load the MOF_topology_dict flie in database folder
        preparation = prepare()
        self.preparation = preparation

    def show_available_mof_families(self):
        # show available mof families and hints for preparation
        self.preparation.list_mof_family()
        print(
            "MOF builder is initialized, please prepare the building material by calling the preparation driver"
        )
        print("1.preparation.select_mof_family()")
        print("2.preparation.select_node_metal()")
        print("3.preparation.use_dummy_node()")
        print("4.preparation.fetch_node()")
        print("5.preparation.fetch_linker()")
        print("***check the preparation status by calling preparation_check()***")

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
            self.node_target_type = preparation.node_metal
            self.linker_pdb = preparation.selected_linker_edge_pdb
            self.linker_center_pdb = (
                preparation.selected_linker_center_pdb
            )  # could be None if ditopic linker
            self.linker_topic = preparation.linker_topic
            self.linker_xyz = preparation.linker_xyz
        else:
            print("Error: Could not find the required files")
            print("Please redo the preparation steps")
            return

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

    def set_remove_node_list(self, remove_node_list):
        self.remove_node_list = remove_node_list

    def set_remove_edge_list(self, remove_edge_list):
        self.remove_edge_list = remove_edge_list

    def save_optimized_rotations(self, filename):
        self.optimized_rotations_filename = filename

    def set_use_saved_optimized_rotations_npy(self, saved_optimized_rotations):
        saved_optimized_rotations = saved_optimized_rotations + ".npy"

        if os.path.exists(saved_optimized_rotations):
            self.saved_optimized_rotations = np.load(
                saved_optimized_rotations, allow_pickle=True
            )
            print("Optimized rotations are loaded from: ", saved_optimized_rotations)
        else:
            print(
                f"Could not find the saved optimized rotations:  {saved_optimized_rotations} will start the optimization from the beginning"
            )
            pass

    def set_use_saved_rotations_as_initial_guess(
        self, use_saved_rotations_as_initial_guess
    ):
        """
        use the saved optimized rotations as initial guess
        """
        if hasattr(self, "saved_optimized_rotations"):
            self.use_saved_rotations_as_initial_guess = (
                use_saved_rotations_as_initial_guess
            )
        else:
            print(
                "saved_optimized_rotations is not found, will start the optimization from the beginning"
            )
            pass

    def build(self):
        # check before building
        if not hasattr(self, "supercell"):
            self.supercell = [1, 1, 1]

        if self.linker_topic == 2:
            print("ditopic mof builder driver is called")
            start_time = time.time()
            linker_pdb = self.linker_pdb
            template_cif = self.template_cif
            node_pdb = self.node_pdb
            supercell = self.supercell
            self.net = net_optimizer()
            if hasattr(self, "rotation_optimizer_maxfun"):
                self.net.set_rotation_optimizer_maxfun(self.rotation_optimizer_maxfun)
            if hasattr(self, "rotation_optimizer_maxiter"):
                self.net.set_rotation_optimizer_maxiter(self.rotation_optimizer_maxiter)
            if hasattr(self, "rotation_optimizer_method"):
                self.net.set_rotation_optimizer_method(self.rotation_optimizer_method)
            if hasattr(self, "rotation_optimizer_eps"):
                self.net.set_rotation_optimizer_eps(self.rotation_optimizer_eps)
            if hasattr(self, "rotation_optimizer_iprint"):
                self.net.set_rotation_optimizer_iprint(self.rotation_optimizer_iprint)
            if hasattr(self, "rotation_optimizer_display"):
                self.net.set_rotation_optimizer_display(self.rotation_optimizer_display)

            if hasattr(self, "saved_optimized_rotations"):
                self.net.load_saved_optimized_rotations(self.saved_optimized_rotations)
            if hasattr(self, "optimized_rotations_filename"):
                self.net.to_save_optimized_rotations(self.optimized_rotations_filename)
            if hasattr(self, "use_saved_rotations_as_initial_guess"):
                self.net.use_saved_rotations_as_initial_guess(
                    self.use_saved_rotations_as_initial_guess
                )

            self.net.analyze_template_ditopic(template_cif)
            self.net.node_info(node_pdb)
            self.net.linker_info(linker_pdb)
            self.net.optimize()
            print("-" * 80)
            print(
                " " * 15,
                "Building time cost: %.5f seconds " % (time.time() - start_time),
            )
            print("-" * 80)

            self.net.set_supercell(supercell)
            self.net.place_edge_in_net()
            self.net.make_supercell_ditopic()
            self.net.make_eG_from_supereG_ditopic()
            self.net.main_frag_eG()
            self.net.make_supercell_range_cleaved_eG()
            self.net.add_xoo_to_edge_ditopic()
            self.net.find_unsaturated_node_eG()
            if hasattr(self, "node_termination"):
                self.net.set_node_terminamtion(self.node_termination)
            # default termination is methyl in data folder
            self.net.set_node_terminamtion(
                os.path.join(
                    self.preparation.data_path, "terminations_database/methyl.pdb"
                )
            )
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
            self.net = net_optimizer()

            if hasattr(self, "rotation_optimizer_maxfun"):
                self.net.set_rotation_optimizer_maxfun(self.rotation_optimizer_maxfun)
            if hasattr(self, "rotation_optimizer_maxiter"):
                self.net.set_rotation_optimizer_maxiter(self.rotation_optimizer_maxiter)
            if hasattr(self, "rotation_optimizer_method"):
                self.net.set_rotation_optimizer_method(self.rotation_optimizer_method)
            if hasattr(self, "rotation_optimizer_eps"):
                self.net.set_rotation_optimizer_eps(self.rotation_optimizer_eps)
            if hasattr(self, "rotation_optimizer_iprint"):
                self.net.set_rotation_optimizer_iprint(self.rotation_optimizer_iprint)
            if hasattr(self, "rotation_optimizer_display"):
                self.net.set_rotation_optimizer_display(self.rotation_optimizer_display)
            if hasattr(self, "saved_optimized_rotations"):
                self.net.load_saved_optimized_rotations(self.saved_optimized_rotations)
            if hasattr(self, "optimized_rotations_filename"):
                self.net.to_save_optimized_rotations(self.optimized_rotations_filename)
            if hasattr(self, "use_saved_rotations_as_initial_guess"):
                self.net.use_saved_rotations_as_initial_guess(
                    self.use_saved_rotations_as_initial_guess
                )

            self.net.analyze_template_multitopic(template_cif)
            self.net.node_info(node_pdb)
            self.net.linker_info(linker_pdb)
            self.net.linker_center_info(linker_center_pdb)
            self.net.optimize()
            print("-" * 80)
            print(
                " " * 15,
                "Building time cost: %.5f seconds " % (time.time() - start_time),
            )
            print("-" * 80)

            self.net.set_supercell(supercell)
            self.net.place_edge_in_net()
            self.net.make_supercell_multitopic()
            self.net.make_eG_from_supereG_multitopic()
            self.net.main_frag_eG()
            self.net.make_supercell_range_cleaved_eG()
            print("_" * 80)
            print("debugging1,before add_xoo_to_edge_multitopic")
            # debugging
            for i in self.net.eG.nodes["EDGE_-9"]["f_points"]:
                if nn(i[0]) == "X":
                    print(i)

            self.net.add_xoo_to_edge_multitopic()
            print("_" * 80)
            print("debugging2,after add_xoo_to_edge_multitopic")
            for i in self.net.eG.nodes["EDGE_-9"]["f_points"]:
                if nn(i[0]) == "X":
                    print(i)

            self.net.find_unsaturated_node_eG()
            if hasattr(self, "node_termination"):
                self.net.set_node_terminamtion(self.node_termination)
            # default termination is methyl in data folder
            self.net.set_node_terminamtion(
                os.path.join(
                    self.preparation.data_path, "terminations_database/methyl.pdb"
                )
            )
            self.net.add_terminations_to_unsaturated_node()
            self.net.remove_xoo_from_node()
            # self.net = net

    def set_gro_name(self, gro_name):
        self.gro_name = gro_name

    def write_gro(self):
        if hasattr(self, "saved_eG"):
            if self.supercell == self.saved_supercell:
                print("saved_eG is found, will write the preserved eG")
                self.net.eG = self.saved_eG
                self.net.write_node_edge_node_gro(self.gro_name)
                return

        if not hasattr(self, "gro_name"):
            self.gro_name = (
                "mof_"
                + str(self.mof_family.split(".")[0])
                + "_"
                + os.path.basename(self.linker_xyz.strip(".xyz"))
            )
            print("gro_name is not set, will be saved as: ", self.gro_name + ".gro")
        print("writing gro file")
        print("nodes:", len(self.net.eG.nodes()), "edges:", len(self.net.eG.edges()))
        self.net.write_node_edge_node_gro(self.gro_name)
        # temp_save_eGterm_gro(net.eG,net.sc_unit_cell) #debugging

    # functions are under construction
    def make_defects_missing(self):
        self.saved_eG = self.net.eG.copy()  # save the original eG before making defects
        self.saved_supercell = self.supercell
        self.defective_net = self.net
        print(
            "saved_eG is saved",
            "nodes: ",
            len(self.saved_eG.nodes),
            "edges: ",
            len(self.saved_eG.edges),
        )

        dG = self.net.eG.copy()
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
            remove_node_list, self.net.nodes_eG
        )
        to_remove_edges_name = extract_node_name_from_gro_resindex(
            remove_edge_list, self.net.edges_eG
        )

        for node_name in to_remove_nodes_name:
            dG.remove_node(node_name)
        for edge_name in to_remove_edges_name:
            neighbors = list(dG.neighbors(edge_name))
            if len(neighbors) == 2:  # ditopic linker case
                dG.remove_edge(neighbors[0], neighbors[1])
            dG.remove_node(edge_name)

        print(
            "defective eG is updated",
            "nodes: ",
            len(dG.nodes),
            "edges: ",
            len(dG.edges),
        )

        self.defective_net.eG = dG
        self.defective_net.main_frag_eG()
        # sort subgraph by connectivity
        self.defective_net.make_supercell_range_cleaved_eG()
        self.defective_net.find_unsaturated_node_eG()
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
                exchange_edge_list, self.net.edges_eG
            )

        # if hasattr(self, 'to_exchange_node_pdb'):
        #    to_exchange_node_pdb = self.to_exchange_node_pdb
        if hasattr(self, "to_exchange_edge_pdb"):
            to_exchange_edge_pdb = self.to_exchange_edge_pdb
        # TODO:

        if hasattr(self, "exchange_edge_list") and hasattr(
            self, "to_exchange_edge_pdb"
        ):
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

    def write_defect_gro(self):
        # if not hasattr(self, 'defective_net'):
        #    print('defective_net is not set')
        #    print('make_defects_missing() or make_defects_exchange() should be called before write_defect_gro(), or you can write with write_gro()')
        #    return

        if not hasattr(self, "defect_gro_name"):
            self.defect_gro_name = (
                "defective_mof_"
                + str(self.mof_family.split(".")[0])
                + "_"
                + os.path.basename(self.linker_xyz.strip(".xyz"))
            )
            print(
                "defect_gro_name is not set, will be saved as: ",
                self.defect_gro_name + ".gro",
            )
        print("writing defective gro file")

        self.defective_net.write_node_edge_node_gro(self.defect_gro_name)

    def write_defective_split_node_gro_again(self, gro_name):
        if not self.preparation.dummy_node:
            print("dummy node is not used, splitting node is not possible")
            return
        print("splitting node and saving gro again, called after write_defect_gro()")
        nodes_eG = self.defective_net.nodes_eG
        edges_eG = self.defective_net.edges_eG
        terms_eG = self.defective_net.terms_eG
        node_split_dict = make_dummy_split_node_dict(self.node_pdb)
        nodes_eGarr = np.vstack(nodes_eG)
        metals_list, hho_list, ho_list, o_list = rename_node_arr(
            node_split_dict, nodes_eGarr
        )

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
                merge_metal_list_to_node_array(
                    merged_split_node_edge_term, splitted_node, line_num, res_count
                )
            )

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
            node_split_dict, nodes_eGarr
        )

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
                merge_metal_list_to_node_array(
                    merged_split_node_edge_term, splitted_node, line_num, res_count
                )
            )

        save_node_edge_term_gro(merged_split_node_edge_term, gro_name)

        print("metal_res_num: ", len(metals_list))
        print("hho_res_num: ", len(hho_list))
        print("ho_res_num: ", len(ho_list))
        print("o_res_num: ", len(o_list))
        print("edge_res_num: ", len(edges_eG))
        print("term_res_num: ", len(terms_eG))


if __name__ == "__main__":
    mof = MofBuilder()
    mof.preparation.select_mof_family("UiO-67")
    mof.preparation.select_node_metal("Zr")
    mof.preparation.use_dummy_node(True)
    mof.preparation.fetch_node()
    mof.preparation.fetch_linker("database/linker4test/ndi.xyz")
    mof.preparation_check()
    mof.set_rotation_optimizer_display(True)

    mof.set_use_saved_optimized_rotations_npy("rota")
    mof.set_use_saved_rotations_as_initial_guess(True)
    # save or update optimized rotations to numpy file for later use
    mof.save_optimized_rotations("rota")
    mof.set_supercell([1, 1, 1])
    mof.build()
    mof.write_gro()
    mof.set_remove_edge_list([30, 31, 32])
    mof.set_remove_node_list([])
    mof.make_defects_missing()
    mof.write_defect_gro()

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

from .molecule import Molecule

from .mofutils import (
    nodepdb2G,
    add_dummy_atoms_nodepdb,
    process_linker_molecule,
    fetch_pdbfile,
    read_mof_top_dict,
    copy_file,
)

from .environment import get_data_path


class MofPreparer:

    def __init__(self):
        # clean up
        if hasattr(self, "mof_top_dict"):
            del self.mof_top_dict
        if hasattr(self, "data_path"):
            del self.data_path

        data_path = get_data_path()
        mof_top_dict = read_mof_top_dict(data_path)
        self.mof_top_dict = mof_top_dict
        self.data_path = data_path

    def list_mof_family(self):
        # print mof_top_dict keys fit to screen
        print("Available MOF Family:")
        print(" ".join(self.mof_top_dict.keys()))

    def select_template_dir(self, template_dir):
        self.set_template_dir = template_dir
        print(f"{template_dir} is selected for template cif files")

    def select_mof_family(self, mof_family):
        self.mof_family = mof_family
        self.node_connectivity = self.mof_top_dict[mof_family][
            "node_connectivity"]
        self.linker_topic = self.mof_top_dict[mof_family]["linker_topic"]
        self.template_cif = self.mof_top_dict[mof_family]["topology"] + ".cif"
        # check if template cif exists
        print(f"{mof_family} selected")
        print(
            f"available metal nodes: {self.mof_top_dict[mof_family]['metal']}")
        print("please select a metal node")
        if not hasattr(self, "set_template_dir"):
            self.set_template_dir = str(
                Path(self.data_path, "template_database"))  # default
            print(f"will search template cif files in {self.set_template_dir}")

        template_cif_file = str(Path(self.set_template_dir, self.template_cif))

        if not Path(template_cif_file).exists():
            print(f"{self.template_cif} not found in template_database")
            print(
                "please download the template cif file and put it in template_database folder, and try again"
            )
            print(
                "or select another MOF family, or upload the template cif file with submit_template method"
            )
            return
        else:
            print(f"{self.template_cif} is found in template_database")
            print(f"{self.template_cif} will be used for MOF building")
            self.selected_template_cif_file = template_cif_file

    def submit_template(
        self,
        template_cif,
        mof_family,
        template_mof_node_connectivity,
        template_node_metal,
        template_linker_topic,
    ):
        # add this item to mof_top_dict in data path
        # check if template cif exists
        if not Path(template_cif).exists():
            print(f"{template_cif} not found")
            print("please upload the template cif file and try again")
            return
        if template_cif.split(".")[1] != "cif":
            print("please upload a cif file")
            return
        if not isinstance(template_mof_node_connectivity, int):
            print("please enter an integer for node connectivity")
            return
        if not isinstance(template_node_metal, str):
            print("please enter a string for node metal")
            return
        if not isinstance(template_linker_topic, int):
            print("please enter an integer for linker topic")
            return

        if mof_family in self.mof_top_dict.keys():
            print(
                f"{mof_family} already exists in the database, the template you submitted will not be used"
            )
            print(
                "if you want to use the template you submitted, please go to the database and replace the existing template"
            )
            return
        else:
            self.mof_top_dict[mof_family] = {
                "node_connectivity": template_mof_node_connectivity,
                "metal": [template_node_metal],
                "linker_topic": template_linker_topic,
                "topology": template_cif.split(".")[0],
            }
            print(f"{mof_family} is added to the database")
            print(f"{mof_family} will be used for MOF building")
            # rewrite mof_top_dict file
            with open(str(Path(self.data_path, "MOF_topology_dict")),
                      "w") as fp:
                head = "MOF            node_connectivity    metal     linker_topic     topology \n"
                fp.write(head)
                for key in self.mof_top_dict.keys():
                    for met in self.mof_top_dict[key]["metal"]:
                        # format is 10s for string and 5d for integer
                        line = "{:15s} {:^16d} {:^12s} {:^12d} {:^18s}".format(
                            key,
                            self.mof_top_dict[key]["node_connectivity"],
                            met,
                            self.mof_top_dict[key]["linker_topic"],
                            self.mof_top_dict[key]["topology"],
                        )
                        fp.write(line + "\n")
            print("mof_top_dict file is updated")
            return str(Path(self.data_path, "template_database", template_cif))

    def select_node_metal(self, node_metal):
        self.node_metal = node_metal
        print(f"{node_metal} node is selected")

    def use_dummy_node(self, dummy_node):
        self.dummy_node = dummy_node

    def fetch_node(self):
        if not hasattr(self, "node_metal"):
            print("please select a metal node with select_node_metal method")
            return
        if not hasattr(self, "node_connectivity"):
            print("please select a MOF family with select_mof_family method")
            return
        if not hasattr(self, "dummy_node"):
            print(
                "if you want to use dummy node, please select dummy node type with select_dummy_node method"
            )
            print("now we are going to fetch the node pdb file")
            self.dummy_node = False

        data_path = self.data_path
        nodes_database_path = str(Path(data_path, "nodes_database"))
        node_pdb = fetch_pdbfile(
            nodes_database_path,
            [str(self.node_connectivity) + "c", self.node_metal],
            ["dummy"],
        )
        # if node_termination == 'methyl' :
        # n_term_file = '../data/terminations_database/methyl.pdb'
        for i in range(len(node_pdb)):
            node_pdb_database = str(
                Path(data_path, "nodes_database", node_pdb[i]))
            target_node_path = str(Path("nodes", node_pdb[i]))
            copy_file(node_pdb_database, target_node_path)

        if self.dummy_node:
            nodeG = nodepdb2G(target_node_path, self.node_metal)
            all_lines, dummy_node_file = add_dummy_atoms_nodepdb(
                target_node_path, self.node_metal, nodeG)
            Path(target_node_path).unlink()
            print(target_node_path, "removed")
            print(f"new dummy node file {dummy_node_file} created")
        else:
            print(f"{target_node_path} is fetched")
            print(
                f"default node without dummy atoms will be used{target_node_path}"
            )

        nodes_dir = str(
            Path(target_node_path).parent
        )  #'nodes' #default  #NOTE: this is the path to the nodes
        print(f"nodes will be saved in {nodes_dir}")

        if self.dummy_node:
            keywords = [
                str(self.node_connectivity) + "c", self.node_metal, "dummy"
            ]
            nokeywords = []
        else:
            keywords = [str(self.node_connectivity) + "c", self.node_metal]
            nokeywords = ["dummy"]

        selected_node_pdb_file = fetch_pdbfile(nodes_dir, keywords,
                                               nokeywords)[0]
        self.selected_node_pdb_file = str(
            Path(nodes_dir, selected_node_pdb_file))
        self.save_nodes_dir = nodes_dir
        return self.selected_node_pdb_file

    def set_split_linker_topic(self, split_linker_topic):
        self.split_linker_topic = split_linker_topic
        print(f"{split_linker_topic} is selected for linker splitting")

    def set_save_edge_dir(self, edges_dir):
        self.save_edges_dir = edges_dir
        print(f"{edges_dir} is selected for saving edges")

    def fetch_linker(self, linker_file):
        self.linker_xyz = linker_file
        if not hasattr(self, "save_edges_dir"):
            self.save_edges_dir = (
                "edges"  # default  #NOTE: this is the path to the edges
            )

        # check if path exists
        if not Path(linker_file).exists():
            raise FileNotFoundError(f"{linker_file} not found")
        if not hasattr(self, "split_linker_topic"):
            split_linker_topic = self.linker_topic
            print(f"{split_linker_topic} is selected for linker splitting")
        else:
            split_linker_topic = self.split_linker_topic
            if split_linker_topic != self.linker_topic:
                # warning
                print(
                    f"{split_linker_topic} topic linker is selected for linker splitting"
                )
                print(
                    f"but {self.linker_topic} topic linker is selected for MOF family"
                )
                print("please make sure they are compatible")

        molecule = Molecule.read_xyz_file(linker_file)

        (
            self.linker_center_frag_nodes_num,
            self.linker_center_Xs,
            self.linker_single_frag_nodes_num,
            self.linker_frag_Xs,
            self.selected_linker_center_pdb,
            self.selected_linker_edge_pdb,
        ) = process_linker_molecule(
            molecule,
            split_linker_topic,
            save_nodes_dir=self.save_nodes_dir,
            save_edges_dir=self.save_edges_dir,
        )

        if self.selected_linker_center_pdb:
            if Path(self.selected_linker_center_pdb).exists():
                print(
                    f"linker center fragment {self.selected_linker_center_pdb} is saved and will be used for MOF building"
                )
        else:
            if self.linker_topic > 2:
                print("Warning: linker center fragment is not saved")

        if Path(self.selected_linker_edge_pdb).exists():
            print(
                f"linker edge fragment {self.selected_linker_edge_pdb} is saved and will be used for MOF building"
            )
        else:
            print("Warning: linker edge fragment is not saved")

    def check_status(self):
        if not hasattr(self, "selected_node_pdb_file"):
            print("node pdb file is not set")
            print("please fetch node with fetch_node method")
            return False
        if not hasattr(self, "selected_linker_center_pdb"):
            print("linker center pdb file is not set")
            print("please fetch linker with fetch_linker method")
            return False
        if not hasattr(self, "selected_linker_edge_pdb"):
            print("linker edge pdb file is not set")
            print("please fetch linker with fetch_linker method")
            return False
        if not hasattr(self, "selected_template_cif_file"):
            print("template cif file is not set")
            print("please select a MOF family with select_mof_family method")
            return False
        if not hasattr(self, "node_metal"):
            print("node metal is not set")
            print("please select a metal node with select_node_metal method")
            return False
        if not hasattr(self, "linker_topic"):
            print("linker topic is not set")
            print("please select a MOF family with select_mof_family method")
            return False

        return True

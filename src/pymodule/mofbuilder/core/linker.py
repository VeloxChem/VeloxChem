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
from typing import Any, List, Optional

import numpy as np
import networkx as nx
import mpi4py.MPI as MPI

from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...molecule import Molecule

from ..io.pdb_reader import PdbReader
from ..io.pdb_writer import PdbWriter


class FrameLinker:
    """Process organic linkers (xyz/molecule): build graph, find center/outer fragments and X atoms for connection.

    Handles ditopic and multitopic linkers; produces linker_center_data, linker_outer_data,
    and X data for placement in the net.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        properties: Dict for optional properties.
        filename: Path to linker file (xyz) or None when using molecule.
        target_directory: Directory for output files.
        new_xyzfilename: Output XYZ path.
        linker_connectivity: Number of connection points (2 = ditopic, etc.).
        pdbreader: PdbReader instance.
        pdbwriter: PdbWriter instance.
        _debug: If True, print extra debug messages.
        linker_data: Raw linker atom data.
        lines: Center fragment lines.
        rows: Outer fragment lines.
        save_files: If True, write output files.
        molecule: VeloxChem (or compatible) molecule when set directly.
        metals: List of metal atom indices in linker.
        isolated_atoms_indices: List of isolated atom indices in linker (not in largest connected component).
        linker_center_data: Center fragment atom data (set by create).
        linker_center_X_data: Center X-atom data.
        linker_outer_data: Outer branch(es) atom data.
        linker_outer_X_data: Outer X-atom data.
        center_class: Center node class (e.g. for multitopic).
        center_nodes: Center node identifiers.
        fake_edge: If True, treat as zero-length edge (e.g. pillar).
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[Any] = None,
        filepath: Optional[str] = None,
    ) -> None:
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)
        self.properties = {}
        self.filename = filepath
        self.target_directory = None
        self.new_xyzfilename = None
        self.linker_connectivity = 2
        self.pdbreader = PdbReader(comm=self.comm, ostream=self.ostream)
        self.pdbwriter = PdbWriter(comm=self.comm, ostream=self.ostream)
        self._debug = False
        self.linker_data = None
        self.lines = None  #center fragment lines
        self.rows = None  #outer fragment lines
        self.save_files = False
        self.molecule = None
        self.metals = []
        self.isolated_atoms_indices = []
        self.linker_center_data = None
        self.linker_center_X_data = None
        self.linker_outer_data = None
        self.linker_outer_X_data = None
        self.center_class = None
        self.center_nodes = None
        self.fake_edge = False

    def check_dirs(self, passfilecheck: bool = True) -> None:
        """Create target_directory; optionally assert linker file exists. Set new_pdbfilename if save_files."""
        if not passfilecheck:
            assert_msg_critical(
                Path(self.filename).exists(),
                f"Linker file {self.filename} not found")
        self.target_directory = Path(self.target_directory or Path.cwd())
        self.target_directory.mkdir(parents=True, exist_ok=True)
        base = Path(self.filename).stem
        if self.save_files:
            self.new_pdbfilename = self.target_directory / f"{base}_processed.pdb"

        if not passfilecheck:
            self.ostream.print_info(
                f"Processing linker file: {self.filename} ...")
        else:
            self.ostream.print_separator()
            self.ostream.print_info(f"Processing linker data ...")
        self.ostream.print_info(f"Linker topic: {self.linker_connectivity}")
        self.ostream.flush()

        if self._debug:
            self.ostream.print_info(
                f"Target directory: {self.target_directory}")

    def _create_lG(self, molecule):
        matrix = molecule.get_connectivity_matrix()
        coords = molecule.get_coordinates_in_angstrom()
        labels = molecule.get_labels()
        dist_matrix = molecule.get_distance_matrix_in_angstrom()
        mass_center_bohr = molecule.center_of_mass_in_bohr()
        bohr_to_angstrom = 0.529177
        mass_center_angstrom = np.asarray(mass_center_bohr) * bohr_to_angstrom
        coords -= mass_center_angstrom

        metal_elements_list = {
            "Ag", "Al", "Au", "Ba", "Be", "Bi", "Ca", "Cd", "Ce", "Co", "Cr",
            "Cs", "Cu", "Fe", "Ga", "Gd", "Hf", "Hg", "In", "Ir", "K", "Li",
            "Mg", "Mn", "Na", "Ni", "Pb", "Pd", "Pt", "Rb", "Rh", "Sc", "Sn",
            "Sr", "Ti", "V", "W", "Y", "Zn", "Zr"
        }

        lG = nx.Graph()
        metals = [
            i for i, label in enumerate(labels) if label in metal_elements_list
        ]
        for i, label in enumerate(labels):
            lG.add_node(i, label=label, coords=coords[i])
        for i in range(len(labels)):
            for j in np.where(matrix[i] == 1)[0]:
                if i not in metals and j not in metals:
                    lG.add_edge(i, j, weight=dist_matrix[i, j])
        if self._debug:
            self.ostream.print_info(f"Number of atoms: {len(labels)}")
            self.ostream.print_info(f"Number of metal atoms: {len(metals)}")
            self.ostream.print_info(
                f"Number of bonds in linker graph: {lG.number_of_edges()}")
            self.ostream.flush()
        self.lG = lG
        self.metals = metals
        self.mass_center_angstrom = mass_center_angstrom

    def _find_center_highly_connected_isolated_cycle(self, lG):
        max_frag_num = 0
        min_frag_size_std = float('inf')
        center_cycle = []
        cycles = list(nx.simple_cycles(lG, length_bound=200))
        for cycle in cycles:
            lG_temp = lG.copy()
            lG_temp.remove_nodes_from(cycle)
            frag_num = nx.number_connected_components(lG_temp)
            frag_sizes = [len(f) for f in nx.connected_components(lG_temp)]
            frag_size_std = np.std(frag_sizes)
            if frag_num > max_frag_num or (frag_num == max_frag_num and
                                           frag_size_std < min_frag_size_std):
                max_frag_num = frag_num
                min_frag_size_std = frag_size_std
                center_cycle = cycle
        return center_cycle

    def _find_center_cycle_nodes(self, lG):
        return self._find_center_highly_connected_isolated_cycle(lG)

    def _find_center_nodes_pair(self, lG, center_nodes):
        if len(center_nodes) > 6:
            centers = nx.barycenter(lG)

        pairs = []
        for i in range(len(centers)):
            for j in range(i, len(centers)):
                l = nx.shortest_path_length(lG, centers[i], centers[j])
                if l == 1:
                    pairs.append([centers[i], centers[j]])

        # loop each pair to find center pair
        for p in pairs:
            a = p[0]
            # b = p[1]
            ds = []
            for n in centers:
                if n not in p:
                    d = nx.shortest_path_length(lG, a, n)
                    ds.append(d)
            if len(set(ds)) < len(ds):
                center_pair = p

        return center_pair

    def _check_two_points_center(self, lG, centers):
        if nx.shortest_path_length(lG, centers[0], centers[1]) != 1:
            return False
        G = lG.copy()
        G.remove_edge(centers[0], centers[1])
        return nx.number_connected_components(G) == 2

    def _in_same_cycle(self, G, nodes):
        for cycle in nx.cycle_basis(G):
            if set(nodes).issubset(cycle):
                return cycle
        return None

    def _find_centers(self, lG):
        barycenter = nx.barycenter(lG)
        normalcenter = nx.center(lG)
        if self._debug:
            self.ostream.print_info(
                f"barycenter: {barycenter}, normalcenter: {normalcenter}")
        if set(barycenter).intersection(normalcenter):
            if self._in_same_cycle(lG, normalcenter) is not None:
                return self._in_same_cycle(lG, normalcenter)
            return barycenter
        return barycenter

    def _distinguish_G_centers(self, lG):
        centers = self._find_centers(lG)
        if len(centers) == 1:
            center_class = "onepoint"
            center_nodes = centers
        elif len(centers) == 2:
            if self._check_two_points_center(lG, centers):
                center_class = "twopoints"
                center_nodes = centers
            else:
                lG.remove_edge(centers[0], centers[1])
                center_class = "cycle"
                center_nodes = self._find_center_cycle_nodes(lG)
        else:
            center_class = "cycle"
            center_nodes = self._find_center_cycle_nodes(lG)
        if self._debug:
            self.ostream.print_info(
                f"center_nodes: {center_nodes}, center_class: {center_class}, n_centers: {len(centers)}"
            )
        self.lG = lG
        self.center_nodes = center_nodes
        self.center_class = center_class

    def _classify_nodes(self):
        lG = self.lG.copy()
        for center_ind, center in enumerate(self.center_nodes):
            lengths = nx.single_source_shortest_path_length(lG, center)
            for k in lengths:
                if center_ind == 0 or lengths[k] < lG.nodes[k].get(
                        "cnodes_l", (None, float('inf')))[1]:
                    lG.nodes[k]["cnodes_l"] = (center, lengths[k])
                elif lengths[k] == lG.nodes[k]["cnodes_l"][1]:
                    lG.nodes[k]["cnodes_l"] = (-1, lengths[k])
        self.lG = lG

    def _get_pairX_outer_frag(self, connected_pairXs, outer_frag_nodes):
        for x in list(connected_pairXs):
            pairXs = [connected_pairXs[x][1], connected_pairXs[x][3]]
            if set(pairXs) < set(outer_frag_nodes):
                break
        return pairXs

    def _lines_of_center_frag(self, subgraph_center_frag, Xs_indices, metals,isolated_atoms_indices):
        labels = self.molecule_labels
        coords = self.molecule_coords
        mass_center_angstrom = self.mass_center_angstrom

        count = 1
        lines = []
        Xs = []
        for cn in list(subgraph_center_frag.nodes):
            label = subgraph_center_frag.nodes[cn]["label"]
            coord = subgraph_center_frag.nodes[cn]["coords"]
            if cn not in Xs_indices:
                name = label + str(count)
            else:
                name = "X" + str(count)
                Xs.append(count - 1)
            count += 1
            lines.append([name, label, coord[0], coord[1], coord[2]])

        for cm in metals:
            label = labels[cm]
            coord = coords[cm] - mass_center_angstrom
            name = label + str(count)
            lines.append([name, label, coord[0], coord[1], coord[2]])
            count += 1
        for ci in isolated_atoms_indices:
            label = labels[ci]
            coord = coords[ci] - mass_center_angstrom
            name = label + str(count)
            lines.append([name, label, coord[0], coord[1], coord[2]])
            count += 1

        return lines, Xs

    def create_pdb(self, filename, lines):
        header = f"REMARK   Generated by MOFbuilder Linker module\n"
        data, x_data = self.pdbreader.expand_arr2data(lines)
        self.pdbwriter.write(filename, header=header, lines=data)
        if self._debug:
            self.ostream.print_info(f"Written PDB file: {filename}")
            self.ostream.flush()

    def _lines_of_single_frag(self, subgraph_single_frag, Xs_indices):
        count = 1
        rows = []
        Xs = []
        for sn in list(subgraph_single_frag.nodes):
            label = subgraph_single_frag.nodes[sn]["label"]
            coord = subgraph_single_frag.nodes[sn]["coords"]
            if sn not in Xs_indices:
                name = label + str(count)
            else:
                name = "X" + str(count)
                Xs.append(count - 1)
            count += 1
            rows.append([name, label, coord[0], coord[1], coord[2]])
        return rows, Xs

    @staticmethod
    def _find_boundary_atom(G, boundary_labels=["O"]):
        boundary_atoms = []
        for n in G.nodes:
            if G.nodes[n]["label"] in boundary_labels:
                #check if it belongs to a -COO group or -COOH group
                #if only 1 edge on node n
                neighbors = list(G.neighbors(n))
                if len(neighbors) == 1:
                    boundary_atoms.append((n, neighbors[0]))
        return boundary_atoms

    @staticmethod
    def _in_frag_labels(f_labels, frag_labels):
        for f in frag_labels:
            if sorted(f) == sorted(f_labels):
                return True
        return False

    @staticmethod
    def _find_boundary_frag(G,
                            center_nodes,
                            frag_labels=[["C", "O", "O"], ["C", "O", "O",
                                                           "H"]]):
        boundary_nodes = FrameLinker._find_boundary_atom(G,
                                                         boundary_labels=["O"])
        center_labels = list(
            set([i for i in sum(frag_labels, []) if i not in ["O", "H"]]))
        boundary_frags = {}
        for n, frag_center in boundary_nodes:
            frag = []
            if G.nodes[frag_center]["label"] in center_labels:
                neighbors = list(G.neighbors(frag_center))
                #n should be added directly
                frag = [n, frag_center]
                #find the closest center node
                center_node, center_path = FrameLinker.find_closest_center_node(
                    G, center_nodes, frag_center)

                for check_n in neighbors:
                    if check_n in frag:
                        continue
                    if nx.shortest_path_length(
                            G, source=check_n,
                            target=center_node) < center_path:
                        continue
                    frag.append(check_n)
                    frag.extend(list(G.neighbors(check_n)))

            frag = list(set(frag))
            f_labels = [G.nodes[i]["label"] for i in frag]
            if FrameLinker._in_frag_labels(f_labels, frag_labels):
                boundary_frags["_".join(str(i) for i in frag)] = {
                    "frag_nodes": frag,
                    "frag_center": frag_center,
                    "center_node": center_node,
                    "center_length": center_path,
                    "labels": f_labels
                }
        return boundary_frags

    @staticmethod
    def find_closest_center_node(G, center_nodes, frag_center):
        center_node = None
        min_path = float('inf')
        for c in center_nodes:
            path_len = nx.shortest_path_length(G, source=frag_center, target=c)
            # to avoid the center pair case, like oxalate
            if path_len == 0:
                continue
            if path_len < min_path:
                min_path = path_len
                center_node = c
        center_path = nx.shortest_path_length(G,
                                              source=frag_center,
                                              target=center_node)

        return center_node, center_path

    def process_linker_molecule(self, molecule, linker_connectivity):
        """
        Processes the linker molecule based on the linker_connectivity and center classification.
        Identifies center nodes, Xs (connection points), fragments, and writes PDB files for each fragment.
        """
        if self.save_files:
            save_nodes_dir = Path(self.target_directory, "nodes")
            save_edges_dir = Path(self.target_directory, "edges")
        self.molecule_coords = molecule.get_coordinates_in_angstrom()
        self.molecule_labels = molecule.get_labels()

        # Remove metal atoms from the graph
        self._create_lG(molecule)
        self.lG.remove_nodes_from(self.metals)
        #only keep the largest connected component for center finding, in case there are isolated subsets in the graph
        if not nx.is_connected(self.lG):
            largest_cc = max(nx.connected_components(self.lG), key=len)
            self.isolated_atoms_indices = self.lG.nodes - largest_cc
            self.lG = self.lG.subgraph(largest_cc).copy()
        self._distinguish_G_centers(self.lG)
        # For large cycles, reduce center nodes to a pair
        if linker_connectivity == 2 and len(self.center_nodes) > 6:
            self.center_nodes = self._find_center_nodes_pair(
                self.lG, self.center_nodes)

        # check center class
        self._classify_nodes()
        if self._debug:
            self.ostream.print_info(f"Linker topic: {linker_connectivity}")
            self.ostream.print_info(f"Center class: {self.center_class}")
            self.ostream.print_info(f"Center nodes: {self.center_nodes}")

        # use center nodes to find the farest, single connected boundary atoms, for carboxylate, only O is-CO(OH) -CO(O)
        frag_labels = [["C", "O", "O"], ["C", "O", "O", "H"]]
        potential_frags = FrameLinker._find_boundary_frag(
            self.lG, self.center_nodes, frag_labels)
        if len(potential_frags) > linker_connectivity:
            #check element count
            lengths = [
                potential_frags[i]['center_length'] for i in potential_frags
            ]
            #check if any length has the count same to linker_connectivity
            length_counts = {len: lengths.count(len) for len in set(lengths)}
            selected_length = None
            for lc in length_counts:
                if length_counts[lc] == linker_connectivity:
                    selected_length = lc
                    break
            if selected_length is not None:
                raise ValueError(
                    "Cannot uniquely determine boundary fragments for multitopic linker."
                )

        if linker_connectivity == 2:
            #remove all 2 fragments
            center_Xs = []
            temp_lG = self.lG.copy()
            for f in potential_frags:
                frag = potential_frags[f]['frag_nodes']
                center_Xs.extend([
                    n for n in self.lG.neighbors(
                        potential_frags[f]['frag_center']) if n not in frag
                ])
                temp_lG.remove_nodes_from(frag)

            if temp_lG.number_of_nodes() == 0:
                #exception for very small linkers, use two X from COO and set edge length to 0 later
                self.fake_edge = True

                temp_lG = self.lG.copy()
                for f in potential_frags:
                    frag = potential_frags[f]['frag_nodes']
                    x_node = potential_frags[f]['frag_center']
                    temp_lG.remove_nodes_from([n for n in frag if n != x_node])
                self.lines, _ = self._lines_of_center_frag(
                    temp_lG, center_Xs, self.metals,self.isolated_atoms_indices)
            else:
                self.lines, _ = self._lines_of_center_frag(
                    temp_lG, center_Xs, self.metals,self.isolated_atoms_indices)
            if self.save_files:
                edge_pdb_name = str(Path(save_edges_dir, "diedge"))
                self.create_pdb(edge_pdb_name, self.lines)

        else:
            #out branch fragments
            #center fragment
            branch_outer_Xs = []
            branch_inner_Xs = []
            center_Xs = []
            temp_lG = self.lG.copy()
            for f in potential_frags:
                frag = potential_frags[f]['frag_nodes']
                outer_X = [
                    n for n in self.lG.neighbors(
                        potential_frags[f]['frag_center']) if n not in frag
                ]
                center_X = potential_frags[f]['center_node']
                inner_X = [
                    n for n in self.lG.neighbors(potential_frags[f]
                                                 ['center_node'])
                    if n not in self.center_nodes + frag
                ]

                if not inner_X:
                    self.fake_edge = True
                    #the -coo group is directly connected to center node
                    #keep the C in -COO as X in outer fragment
                    temp_lG.remove_nodes_from([
                        n for n in frag
                        if n != potential_frags[f]['frag_center']
                    ])
                    temp_lG.remove_edge(potential_frags[f]['frag_center'],
                                        center_X)
                    branch_outer_Xs.append(potential_frags[f]['frag_center'])
                    #no inner X but need creat an edge with zero length, the edge is on same point
                    branch_inner_Xs.append(potential_frags[f]['frag_center'])
                    center_Xs.append(center_X)
                else:
                    branch_outer_Xs.extend(outer_X)
                    branch_inner_Xs.extend(inner_X)
                    center_Xs.append(center_X)
                    temp_lG.remove_edge(center_X, inner_X[0])
                    temp_lG.remove_nodes_from(frag)

            frag_nodes = list(
                sorted(nx.connected_components(temp_lG), key=len,
                       reverse=True))
            center_frag_nodes = []
            outer_frag_nodes = []
            for f in frag_nodes:
                if center_frag_nodes and outer_frag_nodes:
                    break
                if set(self.center_nodes) < set(f):
                    center_frag_nodes = f
                else:
                    outer_frag_nodes = f

            self.lines, _ = self._lines_of_center_frag(
                self.lG.subgraph(center_frag_nodes), center_Xs, self.metals,self.isolated_atoms_indices)
            self.rows, self.frag_Xs = self._lines_of_single_frag(
                self.lG.subgraph(outer_frag_nodes),
                branch_outer_Xs + branch_inner_Xs)

            if self.save_files:
                if linker_connectivity == 3:
                    linker_center_node_pdb_name = str(
                        Path(save_nodes_dir, "tricenter"))
                    self.create_pdb(linker_center_node_pdb_name, self.lines)
                    linker_branch_pdb_name = str(
                        Path(save_edges_dir, "triedge"))
                    self.create_pdb(linker_branch_pdb_name, self.rows)

                elif linker_connectivity == 4:
                    linker_center_node_pdb_name = str(
                        Path(save_nodes_dir, "tetracenter"))
                    self.create_pdb(linker_center_node_pdb_name, self.lines)
                    linker_branch_pdb_name = str(
                        Path(save_edges_dir, "tetraedge"))
                    self.create_pdb(linker_branch_pdb_name, self.rows)

                else:
                    linker_center_node_pdb_name = str(
                        Path(save_nodes_dir, "multicenter"))
                    self.create_pdb(linker_center_node_pdb_name, self.lines)
                    linker_branch_pdb_name = str(
                        Path(save_edges_dir, "multiedge"))
                    self.create_pdb(linker_branch_pdb_name, self.rows)

    def create(self, molecule=None):
        if self.save_files:
            assert_msg_critical(
                self.target_directory is not None,
                "Linker: target_dir is not set. Please set the target directory."
            )
        #assert_msg_critical(self.linker_connectivity in [2, 3, 4] or int(self.linker_connectivity) > 4, "Linker: linker_connectivity should be 2, 3, 4 or >4.")

        if molecule is None:
            assert_msg_critical(
                self.filename is not None,
                "Linker: filename is not set. Please set the filename of the linker molecule."
            )
            self.check_dirs()
            self.molecule = Molecule.read_xyz_file(self.filename)
        else:
            self.molecule = molecule

        self.process_linker_molecule(self.molecule, self.linker_connectivity)
        self.linker_center_data, self.linker_center_X_data = self.pdbreader.expand_arr2data(
            self.lines)
        self.linker_outer_data, self.linker_outer_X_data = self.pdbreader.expand_arr2data(
            self.rows)
        self.ostream.print_info("Linker processing completed.")
        if hasattr(self, "new_pdbfilename"):
            self.ostream.print_info(
                f"Processed linker file is saved as: {self.new_pdbfilename}")
        self.ostream.flush()

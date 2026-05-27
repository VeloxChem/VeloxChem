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

from typing import Any, Dict, List, Optional

import numpy as np
import networkx as nx
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...environment import get_data_path
from mpi4py import MPI
import sys
import re

from ..utils.geometry import (
    unit_cell_to_cartesian_matrix,
    fractional_to_cartesian,
    find_edge_pairings,
)

from pathlib import Path

from ..io.basic import pname, is_list_A_in_B, nn
from ..io.pdb_reader import PdbReader
from ..io.pdb_writer import PdbWriter
from ..io.gro_writer import GroWriter
from ..io.xyz_writer import XyzWriter
from ..io.cif_writer import CifWriter
from ..utils.geometry import cartesian_to_fractional, fractional_to_cartesian
from .superimpose import superimpose_rotation_only


class MofWriter:
    """Convert MOF edge graphs to structure files.

    The writer converts graph node and edge coordinates into the common
    MOFBuilder data layout, merges node/linker/termination residues, and writes
    PDB, GRO, XYZ, or CIF files.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        filename: Base filename without extension for output files.
        G: Edge graph to write.
        frame_cell_info: Cell parameters [a, b, c, alpha, beta, gamma].
        supercell_boundary: Optional [x_min, x_max, y_min, y_max, z_min, z_max] for fractional filter.
        sc_unit_cell: 3x3 supercell unit cell matrix.
        xoo_dict: Dict mapping X index to list of O indices (XOO) per node.
        dummy_atom_node_dict: Optional dict of dummy atom counts per node for renaming.
        residues_info: Dict of residue names and counts (e.g. 'EDGE', ';NODE').
        merged_data: Merged Cartesian data (nodes + edges + terms).
        merged_f_data: Merged fractional data.
        nodes_data: List of node data arrays (set by convert_graph_to_data).
        edges_data: List of edge data arrays.
        terms_data: List of termination data arrays.
        cG: Copy of G with 'data' attached to each node/edge.
        _debug: If True, print additional diagnostics.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[Any] = None,
    ) -> None:
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        self.filename = None
        self.G = None
        self.frame_cell_info = None
        self.supercell_boundary = None
        self.sc_unit_cell = None
        self.xoo_dict = None
        self.dummy_atom_node_dict = None

        self.residues_info = {}
        self.merged_data = None
        self.merged_f_data = None
        self._debug = False

    def _remove_xoo_from_node(
        self, G: nx.Graph, xoo_dict: Dict[int, List[int]]
    ) -> nx.Graph:
        """Remove node-side XOO rows and store ``noxoo_f_points``."""
        eG = G.copy()

        all_xoo_indices = []
        for x_ind, oo_ind in xoo_dict.items():
            all_xoo_indices.append(x_ind)
            all_xoo_indices.extend(oo_ind)

        for n in eG.nodes():
            if pname(n) != "EDGE":
                all_f_points = eG.nodes[n]["f_points"]
                noxoo_f_points = np.delete(all_f_points,
                                           all_xoo_indices,
                                           axis=0)
                eG.nodes[n]["noxoo_f_points"] = noxoo_f_points
        return eG

    def convert_graph_to_data(
        self, G: nx.Graph, sc_unit_cell: np.ndarray
    ) -> None:
        """Convert graph coordinates to Cartesian node, edge, and termination data."""
        rG = self._remove_xoo_from_node(G, self.xoo_dict)

        def arr2data(arr, residue_name=None, residue_number=None, note=None):
            if arr is None or len(arr) == 0:
                return None, None
            if isinstance(arr, list):
                arr = np.vstack(arr)

            data = []
            for i in range(len(arr)):
                atom_type = arr[i, 0]
                atom_label = arr[i, 1]
                if atom_label == "Fr":
                    continue
                value_x = float(arr[i, 2])
                value_y = float(arr[i, 3])
                value_z = float(arr[i, 4])
                atom_number = i + 1
                residue_name = "MOL" if residue_name is None else residue_name
                residue_number = 1 if residue_number is None else residue_number
                charge = 0.0
                spin = 0
                note = nn(atom_type) if note is None else note
                data.append([
                    atom_type, atom_label, atom_number, residue_name,
                    residue_number, value_x, value_y, value_z, spin, charge,
                    note
                ])
            data = np.vstack(data)
            return data

        def get_node_data(n, G, sc_unit_cell):
            node_f_points = G.nodes[n]["noxoo_f_points"]
            res_name = G.nodes[n]["name"]
            res_idx = G.nodes[n]["index"]
            node_f_coords = node_f_points[:, 2:5]
            c_coords = fractional_to_cartesian(node_f_coords, sc_unit_cell)
            node_arr = np.hstack((node_f_points[:, 0:2], c_coords))
            node_data = arr2data(node_arr,
                                 residue_name=res_name,
                                 residue_number=res_idx,
                                 note=n)
            return node_data

        def get_edge_data(n, G, sc_unit_cell):
            edge_f_points = np.vstack(
                (G.nodes[n]["f_points"], G.nodes[n]["xoo_f_points"]))
            res_name = G.nodes[n]["name"]
            res_idx = G.nodes[n]["index"]
            edge_f_coords = edge_f_points[:, 2:5]
            c_coords = fractional_to_cartesian(edge_f_coords, sc_unit_cell)
            edge_arr = np.hstack((edge_f_points[:, 0:2], c_coords))
            edge_data = arr2data(edge_arr,
                                 residue_name=res_name,
                                 residue_number=res_idx,
                                 note=n)
            return edge_data

        cG = self._remove_xoo_from_node(G, self.xoo_dict)
        count = 0
        term_count = 0
        nodes_data = []
        terms_data = []
        edges_data = []
        for n in cG.nodes():
            if pname(n) != "EDGE":
                node_data = get_node_data(n, rG, sc_unit_cell)
                cG.nodes[n]["data"] = node_data
                nodes_data.append(node_data)
                if "term_c_points" in cG.nodes[n]:
                    for term_ind_key, c_positions in cG.nodes[n][
                            "term_c_points"].items():
                        term_name = "T" + cG.nodes[n][
                            "name"]
                        term_count -= 1
                        term_data = arr2data(c_positions,
                                             residue_name=term_name,
                                             residue_number=term_count,
                                             note=term_ind_key)
                        terms_data.append(term_data)
            elif pname(n) == "EDGE":
                edge_data = get_edge_data(n, cG, sc_unit_cell)
                cG.nodes[n]["data"] = edge_data
                edges_data.append(edge_data)

        self.nodes_data = nodes_data
        self.edges_data = edges_data
        self.terms_data = terms_data
        self.cG = cG

    def get_merged_data(
        self, dummy_atom_node_dict: Optional[Dict[str, Any]] = None
    ) -> np.ndarray:
        """Merge Cartesian node, edge, and termination data."""
        nodes_data = self.nodes_data
        edges_data = self.edges_data
        terms_data = self.terms_data
        if dummy_atom_node_dict is not None:
            nodes_data = self._rename_node_name(nodes_data,
                                                dummy_atom_node_dict)
        else:
            nodes_data = np.vstack(nodes_data)
        edges_data = np.vstack(edges_data) if len(
            edges_data) > 0 else np.empty((0, 11))
        terms_data = np.vstack(terms_data) if len(
            terms_data) > 0 else np.empty((0, 11))

        merged_data = np.vstack(
            (nodes_data, edges_data, terms_data
             )) if len(edges_data) > 0 or len(terms_data) > 0 else nodes_data
        self.merged_data = merged_data

        nodes_number = len(self.nodes_data)
        edges_number = len(self.edges_data)
        terms_number = len(self.terms_data)

        if dummy_atom_node_dict is not None:
            self.residues_info[
                'METAL'] = nodes_number * dummy_atom_node_dict.get(
                    'METAL_count', 0)
            self.residues_info[
                'HHO'] = nodes_number * dummy_atom_node_dict.get(
                    'HHO_count', 0)
            self.residues_info['HO'] = nodes_number * dummy_atom_node_dict.get(
                'HO_count', 0)
            self.residues_info['O'] = nodes_number * dummy_atom_node_dict.get(
                'O_count', 0)
        self.residues_info[';NODE'] = nodes_number
        self.residues_info['EDGE'] = edges_number
        self.residues_info['TNODE'] = terms_number

        return merged_data

    def convert_graph_to_fcoords_data(
        self, G: nx.Graph, supercell_boundary: List[float]
    ) -> None:
        """Convert graph coordinates to fractional node and edge data."""
        rG = self._remove_xoo_from_node(G, self.xoo_dict)

        def arr2data(arr, residue_name=None, residue_number=None, note=None):
            if arr is None or len(arr) == 0:
                return None, None
            if isinstance(arr, list):
                arr = np.vstack(arr)

            data = []
            for i in range(len(arr)):
                atom_type = arr[i, 0]
                atom_label = arr[i, 1]
                value_x = float(arr[i, 2])
                value_y = float(arr[i, 3])
                value_z = float(arr[i, 4])
                atom_number = i + 1
                residue_name = "MOL" if residue_name is None else residue_name
                residue_number = 1 if residue_number is None else residue_number
                charge = 0.0
                spin = 0
                note = nn(atom_type) if note is None else note
                data.append([
                    atom_type, atom_label, atom_number, residue_name,
                    residue_number, value_x, value_y, value_z, spin, charge,
                    note
                ])
            data = np.vstack(data)
            return data

        def get_node_fcoords_data(n, G):
            node_f_points = G.nodes[n]["noxoo_f_points"]
            res_name = G.nodes[n]["name"]
            res_idx = G.nodes[n]["index"]
            node_data = arr2data(node_f_points,
                                 residue_name=res_name,
                                 residue_number=res_idx,
                                 note=n)
            return node_data

        def get_edge_fcoords_data(n, G):
            edge_f_points = np.vstack(
                (G.nodes[n]["f_points"], G.nodes[n]["xoo_f_points"]))
            res_name = G.nodes[n]["name"]
            res_idx = G.nodes[n]["index"]
            edge_data = arr2data(edge_f_points,
                                 residue_name=res_name,
                                 residue_number=res_idx,
                                 note=n)
            return edge_data

        def check_supercell_box_range(f_coords, box):
            if f_coords.ndim == 2:
                f_coords = np.mean(f_coords, axis=0)
            x, y, z = map(float, f_coords)
            box_xmax, box_ymax, box_zmax = map(float, box)
            if x < 0 or x >= box[0]:
                return False
            if y < 0 or y >= box[1]:
                return False
            if z < 0 or z >= box[2]:
                return False
            if self._debug:
                self.ostream.print_info(f"f_coords: {f_coords}, box: {box}")
                self.ostream.flush()
            return True

        cG = self._remove_xoo_from_node(G, self.xoo_dict)
        count = 0
        term_count = 0
        nodes_data = []
        terms_data = []
        edges_data = []
        for n in cG.nodes():
            if pname(n) != "EDGE":
                if not check_supercell_box_range(cG.nodes[n]["fcoords"],
                                                 supercell_boundary):
                    continue
                node_f_data = get_node_fcoords_data(n, rG)
                cG.nodes[n]["f_data"] = node_f_data
                nodes_data.append(node_f_data)
            elif pname(n) == "EDGE":
                if not check_supercell_box_range(cG.nodes[n]["fcoords"],
                                                 supercell_boundary):
                    continue
                edge_f_data = get_edge_fcoords_data(n, cG)
                cG.nodes[n]["f_data"] = edge_f_data
                edges_data.append(edge_f_data)

        self.nodes_f_data = np.vstack(nodes_data) if len(
            nodes_data) > 0 else np.empty((0, 11))
        self.edges_f_data = np.vstack(edges_data) if len(
            edges_data) > 0 else np.empty((0, 11))
        self.cG = cG

    def get_merged_fcoords_data(self) -> np.ndarray:
        """Merge fractional node and edge data."""
        nodes_f_data = np.vstack(self.nodes_f_data) if len(
            self.nodes_f_data) > 0 else np.empty((0, 11))
        edges_f_data = np.vstack(self.edges_f_data) if len(
            self.edges_f_data) > 0 else np.empty((0, 11))

        merged_f_data = np.vstack(
            (nodes_f_data,
             edges_f_data)) if len(edges_f_data) > 0 else nodes_f_data
        self.merged_f_data = merged_f_data
        return merged_f_data

    def only_get_merged_data(self) -> tuple:
        """Build and return both Cartesian and fractional merged data."""
        self.convert_graph_to_data(self.G, self.sc_unit_cell)
        self.merged_data = self.get_merged_data(self.dummy_atom_node_dict)

        self.convert_graph_to_fcoords_data(self.G, self.supercell_boundary)
        self.merged_f_data = self.get_merged_fcoords_data()
        return self.merged_data, self.merged_f_data

    def write_pdb(self, skip_merge: bool = False) -> None:
        """Write Cartesian merged data to a PDB file."""
        if self.merged_data is None:
            skip_merge = False

        if skip_merge:
            merged_data = self.merged_data
        else:
            self.convert_graph_to_data(self.G, self.sc_unit_cell)
            merged_data = self.get_merged_data(self.dummy_atom_node_dict)
            self.merged_data = merged_data

        pdb_writer = PdbWriter(comm=self.comm, ostream=self.ostream)
        header = "REMARK   Generated by MOFbuilder\n"
        if self.frame_cell_info is not None:
            a, b, c, alpha, beta, gamma = self.frame_cell_info
            space_group = "P 1"
            z_value = 1
            cryst_line = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} {space_group:<11s}{z_value:4d}\n"
            header += cryst_line
        filename = self.filename
        if self._debug:
            self.ostream.print_info(f"writing pdb file to: {filename}")
            self.ostream.flush()
        pdb_writer.write(filepath=filename, header=header, lines=merged_data)

    def write_xyz(self, skip_merge: bool = False) -> None:
        """Write Cartesian merged data to an XYZ file."""
        if self.merged_data is None:
            skip_merge = False

        if skip_merge:
            merged_data = self.merged_data
        else:
            self.convert_graph_to_data(self.G, self.sc_unit_cell)
            merged_data = self.get_merged_data(self.dummy_atom_node_dict)
            self.merged_data = merged_data

        xyz_writer = XyzWriter(comm=self.comm, ostream=self.ostream)
        header = "Generated by MOFbuilder\n"
        filename = self.filename
        if self._debug:
            self.ostream.print_info(f"writing xyz file to: {filename}")
            self.ostream.flush()
        xyz_writer.write(filepath=filename, header=header, lines=merged_data)

    def write_gro(self, skip_merge: bool = False) -> None:
        """Write Cartesian merged data to a GRO file."""
        if self.merged_data is None:
            skip_merge = False
        if skip_merge:
            merged_data = self.merged_data
        else:
            self.convert_graph_to_data(self.G, self.sc_unit_cell)
            merged_data = self.get_merged_data(self.dummy_atom_node_dict)
            self.merged_data = merged_data

        gro_writer = GroWriter(comm=self.comm, ostream=self.ostream)
        header = "Generated by MOFbuilder\n"
        filename = self.filename
        if self._debug:
            self.ostream.print_info(f"writing gro file to: {filename}")
            self.ostream.flush()
        gro_writer.write(filepath=filename,
                         header=header,
                         lines=merged_data,
                         box=self.frame_cell_info)

    def write_cif(
        self,
        skip_merge: bool = False,
        supercell_boundary: Optional[List[float]] = None,
        frame_cell_info: Optional[List[float]] = None,
    ) -> None:
        """Write fractional merged data to a CIF file."""
        if frame_cell_info is None:
            frame_cell_info = self.frame_cell_info
        if supercell_boundary is None:
            supercell_boundary = self.supercell_boundary
        if self.merged_f_data is None:
            skip_merge = False
        if skip_merge:
            merged_f_data = self.merged_f_data
        else:
            self.convert_graph_to_fcoords_data(self.G, supercell_boundary)
            merged_f_data = self.get_merged_fcoords_data()
            self.merged_f_data = merged_f_data

        assert_msg_critical(frame_cell_info is not None,
                            "frame_cell_info is not provided for cif writing")
        assert_msg_critical(
            supercell_boundary is not None,
            "supercell_boundary is not provided for cif writing")

        cif_writer = CifWriter(comm=self.comm, ostream=self.ostream)
        header = "Generated by MOFbuilder\n"
        filename = self.filename

        if self._debug:
            self.ostream.print_info(f"writing cif file to: {filename}")
            self.ostream.flush()
        header = "Generated by MOFbuilder\n"
        cif_writer.write(filepath=filename,
                         header=header,
                         lines=self.merged_f_data,
                         cell_info=frame_cell_info,
                         supercell_boundary=supercell_boundary)

    def _rename_node_name(
        self,
        nodes_data: List[np.ndarray],
        dummy_atom_node_dict: Optional[Dict[str, Any]],
    ) -> np.ndarray:
        """Rename and group dummy node residues by residue type."""
        if dummy_atom_node_dict is None:
            return np.vstack(nodes_data)
        nodes_num = len(nodes_data)
        metal_count = dummy_atom_node_dict["METAL_count"]
        dummy_res_len = int(dummy_atom_node_dict["dummy_res_len"])
        hho_count = dummy_atom_node_dict["HHO_count"]
        ho_count = dummy_atom_node_dict["HO_count"]
        o_count = dummy_atom_node_dict["O_count"]
        metal_num = metal_count * dummy_res_len
        hho_num = hho_count * 3
        ho_num = ho_count * 2
        o_num = o_count * 1

        residue_specs = ([("METAL", dummy_res_len)] * metal_count +
                         [("HHO", 3)] * hho_count + [("HO", 2)] * ho_count +
                         [("O", 1)] * o_count)

        new_name_list = [
            f"{name}_{i+1}" for i, (name, count) in enumerate(residue_specs)
            for _ in range(count)
        ]
        nodes_data = np.vstack(nodes_data)
        if self._debug:
            self.ostream.print_info(
                f"dummy node split dict: {dummy_atom_node_dict}")
            self.ostream.print_info(f"new_name_list: {new_name_list}")
            self.ostream.flush()

        name_col = np.tile(new_name_list, nodes_num).reshape(-1, 1)
        rename_data = np.hstack((nodes_data[:, 0:3], name_col, nodes_data[:,
                                                                          4:]))
        rename_data = rename_data.reshape(nodes_num, -1, 11)
        metals_data = np.vstack(
            rename_data[:, :metal_num, :]) if metal_num > 0 else np.empty(
                (0, 11))
        hhos_data = np.vstack(
            rename_data[:, metal_num:metal_num +
                        hho_num, :]) if hho_num > 0 else np.empty((0, 11))
        hos_data = np.vstack(
            rename_data[:, metal_num + hho_num:metal_num + hho_num +
                        ho_num, :]) if ho_num > 0 else np.empty((0, 11))
        os_data = np.vstack(
            rename_data[:, metal_num + hho_num +
                        ho_num:, :]) if o_num > 0 else np.empty((0, 11))
        ordered_data = np.vstack((metals_data, hhos_data, hos_data, os_data))

        return ordered_data

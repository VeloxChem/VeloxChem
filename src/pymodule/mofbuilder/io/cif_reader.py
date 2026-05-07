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

import re
import numpy as np
from pathlib import Path
from .basic import convert_fraction_to_decimal, remove_bracket, remove_quotes, remove_tail_number, extract_quote_lines
from .basic import find_keyword, extract_xyz_lines
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from mpi4py import MPI
import sys

class CifReader:
    """
    Class for reading CIF files, extracting symmetry and atomic site information,
    and providing atom coordinates in the primitive cell.

    Attributes:
        comm (MPI.Comm): MPI communicator.
        rank (int): Rank of the process.
        nodes (int): Total number of nodes/processes.
        ostream (OutputStream): Output stream for logging.
        filepath (str or Path): Path to the CIF file.
        data (np.ndarray or None): Array with atom data.
        _debug (bool): Toggle debug logging.
        cell_info (list): Unpacked cell dimensions and angles.
        symmetry_sector (list): Lines containing symmetry operations.
        atom_site_sector (list): Lines defining atomic sites.
        net_name (str): CIF network name.
        spacegroup (str): Space group label from CIF.
        hall_number (str): Hall number label from CIF.
        V_con (str): Connected component value from CIF, if present.
        EC_con (str): Edge connectivity value from CIF, if present.
        ciffile_lines (list[str]): Filtered lines from the CIF file.
        fcoords (np.ndarray): Fractional coordinates for atoms.
        target_fcoords (np.ndarray): Specific fractional coordinates for requested atom type.
    """

    def __init__(self, comm: MPI.Comm = None, ostream: OutputStream = None, filepath: str = None):
        """
        Initialize a CifReader instance.

        Args:
            comm (MPI.Comm, optional): MPI communicator. Defaults to MPI.COMM_WORLD.
            ostream (OutputStream, optional): Output stream for logging. If None, selected by rank.
            filepath (str, optional): Path to the CIF file to be read.
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm: MPI.Comm = comm
        self.rank: int = self.comm.Get_rank()
        self.nodes: int = self.comm.Get_size()
        self.ostream: OutputStream = ostream

        self.filepath: str = filepath
        self.data: np.ndarray | None = None

        self._debug: bool = False

    def _extract_value_str_slice(self, s: list) -> int:
        """
        Extracts a signed integer coefficient from a symmetry operation string segment.

        Args:
            s (list): List of string fragments containing symmetry operation.

        Returns:
            int: The signed integer value (e.g., +1, -1).
        """
        if len(s) == 0:
            return 0
        sign = 1
        mul_value = 1
        s_list = list(s[0])

        if "-" in s_list:
            sign = -1
        if "*" in s_list:
            mul_value = s_list[s_list.index("*") - 1]

        return sign * int(mul_value)

    def _extract_value_from_str(self, s: str) -> tuple[int, int, int, float]:
        """
        Extracts x, y, z coefficients and constants from a symmetry operation string.

        Args:
            s (str): The symmetry operation component string.

        Returns:
            tuple[int, int, int, float]: (x_coeff, y_coeff, z_coeff, constant)
        """
        s = re.sub(r" ", "", s)
        s = re.sub(r"(?<=[+-])", ",", s[::-1])[::-1]
        if s[0] == ",":
            s = s[1:]
        s_list = list(s.split(","))
        x_slice = [s_list[i] for i in range(len(s_list)) if "x" in s_list[i]]
        y_slice = [s_list[i] for i in range(len(s_list)) if "y" in s_list[i]]
        z_slice = [s_list[i] for i in range(len(s_list)) if "z" in s_list[i]]
        const_slice = [
            s_list[i] for i in range(len(s_list))
            if "x" not in s_list[i] and "y" not in s_list[i] and "z" not in s_list[i]
        ]
        x_coeff = self._extract_value_str_slice(x_slice)
        y_coeff = self._extract_value_str_slice(y_slice)
        z_coeff = self._extract_value_str_slice(z_slice)

        if len(const_slice) == 0:
            const = 0
        else:
            const = const_slice[0]
            const = convert_fraction_to_decimal(const)

        return x_coeff, y_coeff, z_coeff, const

    def _extract_transformation_matrix_from_symmetry_operator(self, expr_str: str) -> np.ndarray:
        """
        Converts a symmetry operator string into a 4x4 transformation matrix.

        Args:
            expr_str (str): Symmetry operation as string, e.g., 'x, -y, z+1/2'.

        Returns:
            np.ndarray: 4x4 transformation matrix.
        """
        expr_str = str(expr_str)
        expr_str = expr_str.strip("\n")
        expr_str = expr_str.replace(" ", "")
        split_str = expr_str.split(",")
        transformation_matrix = np.zeros((4, 4))
        transformation_matrix[3, 3] = 1
        for i in range(len(split_str)):
            x_coeff, y_coeff, z_coeff, const = self._extract_value_from_str(split_str[i])
            transformation_matrix[i] = [x_coeff, y_coeff, z_coeff, const]

        return transformation_matrix

    def _extract_symmetry_operation_from_lines(self, symmetry_sector: list[str]) -> list[str]:
        """
        Extracts symmetry operation strings from a sector of the CIF file.

        Args:
            symmetry_sector (list[str]): List of lines in the symmetry sector.

        Returns:
            list[str]: List of symmetry operation strings.
        """
        symmetry_operation = []
        for i in range(len(symmetry_sector)):
            pattern = r"([+-]?\d*\.?\d*)\s*([xyz])"
            match = re.search(pattern, symmetry_sector[i])
            if match:
                string = remove_quotes(symmetry_sector[i].strip("\n"))
                no_space_string = string.replace(" ", "")
                symmetry_operation.append(no_space_string)
        if len(symmetry_operation) < 2:
            if self._debug:
                self.ostream.print_info("P1 cell")
                self.ostream.flush()
            symmetry_operation = ["x,y,z"]
        else:
            if self._debug:
                self.ostream.print_info(
                    f"apply {len(symmetry_operation)}  symmetry operation")
                self.ostream.flush()

        return symmetry_operation

    def _fetch_spacegroup_from_cifinfo(self) -> str:
        """
        Fetches the space group name (Hermann-Mauguin) from the CIF info.

        Returns:
            str: The space group name if found; otherwise, 'P1'.
        """
        pattern = r"_symmetry_space_group_name_H-M\s+'([^']+)'"
        match = re.search(pattern, )
        if match:
            return match.group(1)
        else:
            return "P1"

    def _valid_net_name_line(self, line: str) -> str:
        """
        Checks if a line specifies a network name and returns a valid network name.

        Args:
            line (str): Line from CIF file.

        Returns:
            str: The detected network name.
        """
        if re.search(r"net", line):
            potential_net_name = line.split()[0].split("_")[1]
            if re.sub(r"[0-9]", "", potential_net_name) == "":
                return Path(self.filepath).stem
            else:
                return potential_net_name

    def _valid_spacegroup_line(self, line: str) -> str:
        """
        Checks a line for a space group or potential network name.

        Args:
            line (str): Line from CIF file.

        Returns:
            str: The spacegroup name or network name, defaults to 'P1' if not found.
        """
        if re.search(r"_symmetry_space_group_name_H-M", line):
            space_group = re.search(
                r"_symmetry_space_group_name_H-M\s+'([^']+)'", line)[1]
            return space_group
        elif re.search(r"^data_", line) and line.count("_") >= 3:
            potential_net_name = line.split()[0].split("_")[2]
            return potential_net_name
        return "P1"

    def _valid_hallnumber_line(self, line: str) -> str:
        """
        Checks a line for Hall number information.

        Args:
            line (str): Line from CIF file.

        Returns:
            str: The Hall number as found, or '1' if not found.
        """
        if re.search(r"_symmetry_Int_Tables_number", line):
            hall_number = re.search(r"_symmetry_Int_Tables_number\s+(\d+)", line)[1]
            return hall_number
        elif re.search(r"hall_number:\s*(\d+)", line):
            hall_number = re.search(r"hall_number:\s*(\d+)", line)[1]
            return hall_number
        return "1"

    def read_cif(self, cif_file: str = None) -> None:
        """
        Reads a CIF file, storing relevant crystal and atom site sections in the instance.

        Args:
            cif_file (str, optional): Path to the CIF file. If None, uses self.filepath.

        Raises:
            AssertionError: If the CIF file path does not exist.
        """
        net_flag = False
        spacegroup_flag = False
        hallnumber_flag = False
        vcon_flag = False

        if cif_file is not None:
            self.filepath = cif_file
        assert_msg_critical(
            Path(self.filepath).exists(),
            f"cif file {self.filepath} not found")
        if self._debug:
            self.ostream.print_info(f"Reading cif file {self.filepath}")
            self.ostream.flush()

        def valid_line(line: str) -> bool:
            """Helper to test whether a line is meaningful (not comment or empty)."""
            return line.strip() != "" and not line.strip().startswith("#")

        with open(self.filepath, "r") as f:
            lines = f.readlines()
            nonempty_lines = [line for line in lines if valid_line(line)]

        self.ciffile_lines = nonempty_lines

        for line in nonempty_lines[:100]:
            if net_flag & spacegroup_flag & hallnumber_flag & vcon_flag:
                break
            if not net_flag and re.search(r"net", line):
                self.net_name = self._valid_net_name_line(line)
                net_flag = True
            if not spacegroup_flag:
                self.spacegroup = self._valid_spacegroup_line(line)
                if self.spacegroup != "P1":
                    spacegroup_flag = True
            if not hallnumber_flag:
                self.hall_number = self._valid_hallnumber_line(line)
                if self.hall_number != "1":
                    hallnumber_flag = True
            if not vcon_flag and re.search(r"V_con:\s*(\d+)", line):
                self.V_con = re.search(r"V_con:\s*(\d+)", line)[1]
                vcon_flag = True
                self.EC_con = re.search(r"EC_con:\s*(\d+)",
                                        line)[1] if re.search(
                                            r"EC_con:\s*(\d+)", line) else None

        if hasattr(self, 'net_name'):
            self.ostream.print_info(f"Found net name: {self.net_name}")
        if hasattr(self, 'spacegroup'):
            self.ostream.print_info(f"Spacegroup: {self.spacegroup}")

        self.ostream.flush()

        if self._debug:
            if hasattr(self, 'V_con'):
                self.ostream.print_info(f"Found V_con: {self.V_con}")
            if hasattr(self, 'EC_con'):
                self.ostream.print_info(f"Found EC_con: {self.EC_con}")
            self.ostream.flush()

        keyword1 = r"loop_"
        keyword2 = r"x,\s*y,\s*z"
        keyword3 = r"-x"
        # loop_key marks locations of important 'loop_' and 'x, y, z' tags

        loop_key = [0]
        linenumber = 0
        for i in range(len(nonempty_lines)):
            m = find_keyword(keyword1, nonempty_lines[i]) or (
                find_keyword(keyword2, nonempty_lines[i]) and
                (not find_keyword(keyword3, nonempty_lines[i])))

            a = re.search(r"_cell_length_a", nonempty_lines[i])
            b = re.search(r"_cell_length_b", nonempty_lines[i])
            c = re.search(r"_cell_length_c", nonempty_lines[i])
            alpha = re.search(r"_cell_angle_alpha", nonempty_lines[i])
            beta = re.search(r"_cell_angle_beta", nonempty_lines[i])
            gamma = re.search(r"_cell_angle_gamma", nonempty_lines[i])

            if m:
                loop_key.append(linenumber)
            else:
                if a:
                    cell_length_a = remove_bracket(
                        nonempty_lines[i].split()[1])
                elif b:
                    cell_length_b = remove_bracket(
                        nonempty_lines[i].split()[1])
                elif c:
                    cell_length_c = remove_bracket(
                        nonempty_lines[i].split()[1])
                elif alpha:
                    cell_angle_alpha = remove_bracket(
                        nonempty_lines[i].split()[1])
                elif beta:
                    cell_angle_beta = remove_bracket(
                        nonempty_lines[i].split()[1])
                elif gamma:
                    cell_angle_gamma = remove_bracket(
                        nonempty_lines[i].split()[1])

            linenumber += 1
        loop_key.append(len(nonempty_lines))
        loop_key = list(set(loop_key))
        loop_key.sort()

        cell_info = [
            cell_length_a,
            cell_length_b,
            cell_length_c,
            cell_angle_alpha,
            cell_angle_beta,
            cell_angle_gamma,
        ]

        cif_sectors = []
        for i in range(len(loop_key) - 1):
            cif_sectors.append(nonempty_lines[loop_key[i]:loop_key[i + 1]])
        for i in range(len(cif_sectors)):
            if re.search(keyword2, cif_sectors[i][0]):
                symmetry_sector = cif_sectors[i]

            if len(cif_sectors[i]) > 1:
                if re.search(r"_atom_site_label\s+", cif_sectors[i][1]):
                    atom_site_sector = cif_sectors[i]
        self.cell_info = cell_info
        self.symmetry_sector = symmetry_sector
        self.atom_site_sector = atom_site_sector

    def _limit_value_0_1(self, new_array_metal_xyz: np.ndarray) -> np.ndarray:
        """
        Ensures all entries are within [0,1) using modulo.

        Args:
            new_array_metal_xyz (np.ndarray): Input coordinates.

        Returns:
            np.ndarray: Coordinates wrapped to [0,1).
        """
        new_array_metal_xyz = np.mod(new_array_metal_xyz, 1)
        return new_array_metal_xyz

    def _wrap_fccords_to_0_1(self, fccords: np.ndarray) -> np.ndarray:
        """
        Center/wraps fractional coordinates to [0,1) and removes duplicated sites.

        Args:
            fccords (np.ndarray): Fractional coordinates.

        Returns:
            np.ndarray: Unique, wrapped fractional coordinates.
        """
        fccords = np.unique(np.array(fccords, dtype=float), axis=0)
        fccords = self._limit_value_0_1(fccords)
        fccords += 0.5
        fccords = self._limit_value_0_1(fccords)
        fccords += -0.5
        fccords = np.unique(np.array(fccords, dtype=float), axis=0)
        return fccords

    def _apply_sym_operator(
        self,
        symmetry_operations: list[str],
        array_metal_xyz: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Applies a list of symmetry operations to the coordinates and returns all unique sites.

        Args:
            symmetry_operations (list[str]): The symmetry operations in CIF string syntax.
            array_metal_xyz (np.ndarray): Coordinates to transform.

        Returns:
            tuple[np.ndarray, np.ndarray]: unique coordinates array and corresponding unique indices.
        """
        array_metal_extend_xyz = np.hstack(
            (array_metal_xyz, np.ones((len(array_metal_xyz), 1))))
        cell_array_metal_xyz = np.empty((0, 3))
        for sym_line_idx in range(len(symmetry_operations)):
            transfromation_matrix = self._extract_transformation_matrix_from_symmetry_operator(
                symmetry_operations[sym_line_idx])
            new_extend_xyz = np.matmul(transfromation_matrix,
                                       array_metal_extend_xyz.T).T
            new_xyz = new_extend_xyz[:, 0:3]
            cell_array_metal_xyz = np.vstack((cell_array_metal_xyz, new_xyz))

        round_cell_array_metal_xyz = np.round(
            self._wrap_fccords_to_0_1(cell_array_metal_xyz), 4)
        _, unique_indices = np.unique(round_cell_array_metal_xyz,
                                      axis=0,
                                      return_index=True)
        unique_indices.sort()
        unique_metal_array = round_cell_array_metal_xyz[unique_indices]

        return unique_metal_array, unique_indices

    def _extract_atoms_fcoords_from_lines(
        self,
        atom_site_sector: list[str]
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Extracts atom type, label, and fractional coordinates from the atom_site sector.

        Args:
            atom_site_sector (list[str]): Lines defining atomic sites in the CIF.

        Returns:
            tuple[np.ndarray, np.ndarray]: (array_atom: [atom_type, atom_label], array_xyz: [x, y, z])
        """
        atom_site_lines = []
        keyword = r"_"
        for line in atom_site_sector:
            m = re.search(keyword, line)
            if m is None:
                atom_site_lines.append(line)

        array_atom = np.zeros((len(atom_site_lines), 2), dtype=object)
        array_xyz = np.zeros((len(atom_site_lines), 3))

        for i in range(len(atom_site_lines)):
            for j in [0, 1, 2, 3, 4]:  # only need atom type, atom label, x,y,z
                if j < 2:
                    array_atom[i, j] = remove_tail_number(
                        atom_site_lines[i].split()[j])
                else:
                    array_xyz[i, (j - 2)] = remove_bracket(
                        atom_site_lines[i].split()[j])
        if self._debug:
            self.ostream.print_info(
                f"Found {len(array_atom)} atoms in atom_site_sector")
            self.ostream.print_info(f"array_atom {array_atom}")
            self.ostream.print_info(f"array_xyz {array_xyz}")
            self.ostream.flush()
        return array_atom, array_xyz

    def get_type_atoms_fcoords_in_primitive_cell(
        self,
        target_type: str = None
    ) -> tuple[list, np.ndarray, np.ndarray]:
        """
        Get all atoms of a particular type in the primitive cell, applying symmetry.

        Args:
            target_type (str, optional): Atom type to select (e.g., "V", "E"). If None, uses all.

        Returns:
            tuple[list, np.ndarray, np.ndarray]: (cell_info, atom data, requested atom fractional coordinates)
        """
        array_atom, array_xyz = self._extract_atoms_fcoords_from_lines(
            self.atom_site_sector)
        if target_type is None:
            self.ostream.print_info(
                f"target_type not specified, use {target_type} as default")
            self.ostream.flush()
        if len(self.symmetry_sector) > 1:  # need to apply symmetry operations
            array_metal_xyz = array_xyz[array_atom[:, 0] == target_type]
            array_metal_xyz = np.round(array_metal_xyz, 4)
            symmetry_sector_neat = extract_quote_lines(self.symmetry_sector)
            if len(symmetry_sector_neat) < 2:  # if no quote, then find x,y,z
                symmetry_sector_neat = extract_xyz_lines(self.symmetry_sector)
            symmetry_operations = self._extract_symmetry_operation_from_lines(
                symmetry_sector_neat)
            no_sym_array_metal_xyz, no_sym_indices = self._apply_sym_operator(
                symmetry_operations, array_metal_xyz)
            array_metal_xyz_final = no_sym_array_metal_xyz
            array_atom = np.tile(array_atom,
                                 (len(symmetry_operations), 1))[no_sym_indices]

        else:
            array_metal_xyz = array_xyz[array_atom[:, 0] == target_type]
            array_metal_xyz_final = np.round(array_metal_xyz, 4)

        self.fcoords = self._wrap_fccords_to_0_1(array_metal_xyz_final)
        self.target_fcoords = self._wrap_fccords_to_0_1(array_metal_xyz_final)

        self.data = []
        for i in range(len(self.fcoords)):
            atom_number = i + 1
            atom_type = array_atom[i, 0] + str(atom_number)
            atom_label = array_atom[i, 1]
            residue_name = 'MOL'
            residue_number = 1
            x = self.fcoords[i, 0]
            y = self.fcoords[i, 1]
            z = self.fcoords[i, 2]
            spin = 1.00
            charge = 0.0
            note = ''
            self.data.append([
                atom_type, atom_label, atom_number, residue_name,
                residue_number, x, y, z, spin, charge, note
            ])
        self.data = np.vstack(self.data)
        if self._debug:
            self.ostream.print_info(
                f"Found {len(self.fcoords)} {target_type} atoms in primitive cell"
            )
            self.ostream.print_info(f"fcoords {self.fcoords}")
            self.ostream.flush()

        return self.cell_info, self.data, self.target_fcoords

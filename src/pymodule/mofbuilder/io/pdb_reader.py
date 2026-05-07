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

import numpy as np
from pathlib import Path
from typing import Optional, Any, Tuple
from .basic import nn
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from mpi4py import MPI
import sys




class PdbReader:
    """Class for reading and processing atomic coordinates from PDB files.

    Attributes:
        comm (MPI.Comm): MPI communicator.
        rank (int): MPI rank of current process.
        nodes (int): Total number of MPI processes.
        ostream (OutputStream): Output stream for logging.
        filepath (str | Path | None): Path to the PDB file.
        com_target_type (str): Atom type used for center-of-mass calculation.
        data (np.ndarray | None): Parsed atom data as array;
            see row format below.
        X_data (np.ndarray | None): Rows corresponding to atom type "X".
        node_atoms (np.ndarray | None): Atom types and labels for processed nodes.
        node_ccoords (np.ndarray | None): Coordinates for all node atoms, recentered.
        node_x_ccoords (np.ndarray | None): Coordinates for atoms of type "X", recentered.
        _debug (bool): Toggle debug-level logging.

    Methods:
        read_pdb: Reads the PDB file and extracts atom data.
        expand_arr2data: Converts an array of atom records into the library's standard format.
        process_node_pdb: Specialized recentering and extraction for node atoms.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[OutputStream] = None,
        filepath: Optional[str] = None,
    ):
        """Initialize a PdbReader instance.

        Args:
            comm (MPI.Comm, optional): MPI communicator. Defaults to MPI.COMM_WORLD.
            ostream (OutputStream, optional): Output stream for logging and info.
            filepath (str, optional): Path to the PDB file.
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream
        self.filepath: Optional[str] = filepath
        self.com_target_type: str = "X"
        self.data: Optional[np.ndarray] = None

        self.X_data: Optional[np.ndarray] = None
        self.node_atoms: Optional[np.ndarray] = None
        self.node_ccoords: Optional[np.ndarray] = None
        self.node_x_ccoords: Optional[np.ndarray] = None

        self._debug: bool = False

    def read_pdb(
        self,
        filepath: Optional[str] = None,
        recenter: bool = True,
        com_type: Optional[str] = None,
    ) -> None:
        """Read a PDB file, extract atoms, and optionally recenter coordinates.

        Populates self.data with shape (N, 11): [atom_type, atom_label, atom_number,
        residue_name, residue_number, value_x, value_y, value_z, occupancy, b_factor, note].

        Args:
            filepath (str, optional): Path to the PDB file. Overrides instance filepath if given.
            recenter (bool, optional): If True (default), recenter coordinates to center-of-mass.
            com_type (str, optional): Atom type for COM; if None, uses all atoms.

        Raises:
            AssertionError: If file does not exist.

        Note:
            Only ATOM/HETATM records are parsed; format assumes canonical PDB column mapping.
        """
        if filepath is not None:
            self.filepath = filepath
        assert_msg_critical(
            Path(self.filepath).exists(),
            f"pdb file {self.filepath} not found",
        )
        if self._debug:
            self.ostream.print_info(f"Reading pdb file {self.filepath}")

        inputfile = str(self.filepath)
        with open(inputfile, "r") as fp:
            lines = fp.readlines()

        data = []
        for line in lines:
            line = line.strip()
            if len(line) > 0:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_number = int(line[6:11]) if line[6:11].strip() else 1
                    atom_type = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip() if line[21].strip() else "A"
                    residue_number = int(line[22:26]) if line[22:26].strip() else 1
                    value_x = float(line[30:38])
                    value_y = float(line[38:46])
                    value_z = float(line[46:54])
                    occupancy = float(line[54:60]) if line[54:60].strip() else 0.0
                    b_factor = float(line[60:66]) if line[60:66].strip() else 0.0
                    atom_label = line[76:78].strip()
                    charge = line[78:80].strip()
                    note = nn(atom_type)

                    data.append([
                        atom_type,         # Assigned atom name, e.g., "CA"
                        atom_label,        # Element symbol, e.g., "C"
                        atom_number,       # Serial number
                        residue_name,      # Residue name
                        residue_number,    # Residue index
                        value_x, value_y, value_z,
                        occupancy,         # Occupancy
                        b_factor,          # B-factor
                        note               # Custom mapping
                    ])

        self.data = np.vstack(data)

        def type_data(arr: np.ndarray) -> np.ndarray:
            """Ensure proper dtype for numeric and string columns."""
            arr[:, 2] = arr[:, 2].astype(int)
            arr[:, 4] = arr[:, 4].astype(int)
            arr[:, 5:8] = arr[:, 5:8].astype(float)
            arr[:, 8] = arr[:, 8].astype(float)
            arr[:, 9] = arr[:, 9].astype(float)
            return arr

        self.data = type_data(self.data)

        if recenter:
            # Use all atoms for center calculation if com_type not specified or not present.
            if com_type is None:
                com_type_ccoords = self.data[:, 5:8].astype(float)
            else:
                if com_type not in self.data[:, -1]:
                    self.ostream.print_warning(
                        f"com_type {com_type} not in the pdb file, use all atoms to calculate com"
                    )
                    com_type_ccoords = self.data[:, 5:8].astype(float)
                else:
                    com_type_ccoords = self.data[
                        self.data[:, -1] == com_type][:, 5:8].astype(float)
            com = np.mean(com_type_ccoords, axis=0)
            if self._debug:
                self.ostream.print_info(
                    f"Center of mass type {com_type} at {com}")
            self.data[:, 5:8] = self.data[:, 5:8].astype(float) - com

        self.X_data = self.data[self.data[:, -1] == "X"]

    def expand_arr2data(
        self, arr: Optional[np.ndarray]
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Convert array of minimal records to standard data format.

        Args:
            arr (np.ndarray | list): Array/list of [atom_type, atom_label, x, y, z].

        Returns:
            Tuple[np.ndarray | None, np.ndarray | None]:
                All parsed data (as for self.data), and extracted rows where note is 'X'.

        Example:
            data, x_data = pdb_reader.expand_arr2data([["Fe", "Fe", 0.0, 0.0, 0.0]])
        """
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
            residue_name = "MOL"
            residue_number = 1
            charge = 0.0
            spin = 0
            note = nn(atom_type)
            data.append([
                atom_type, atom_label, atom_number, residue_name,
                residue_number, value_x, value_y, value_z, spin, charge, note
            ])
        data = np.vstack(data)
        X_data = data[data[:, -1] == 'X']
        return data, X_data

    def process_node_pdb(self) -> None:
        """Read the PDB and prepare node-centered atom sets.

        Parses self.data, recenters coordinates for self.com_target_type
        (default 'X'), and stores node atoms, coordinates, and "X" atoms.

        Note:
            Intended for node PDB files (metal nodes) in the MOF framework.

        Populates:
            - self.node_atoms: [atom_type, atom_label] for all atoms
            - self.node_ccoords: centered coordinates for all atoms
            - self.node_x_ccoords: centered coordinates for "X" atoms
        """
        self.read_pdb()
        node_atoms = self.data[:, 0:2]
        node_ccoords = self.data[:, 5:8]
        node_ccoords = node_ccoords.astype(float)
        com_type_indices = [
            i for i in range(len(node_atoms))
            if nn(node_atoms[i, 0]) == self.com_target_type
        ]
        x_indices = [
            j for j in range(len(node_atoms)) if nn(node_atoms[j, 0]) == "X"
        ]
        node_x_ccoords = self.data[x_indices, 5:8]
        node_x_ccoords = node_x_ccoords.astype(float)
        com_type_ccoords = node_ccoords[com_type_indices]
        com_type = np.mean(com_type_ccoords, axis=0)
        node_ccoords = node_ccoords - com_type
        node_x_ccoords = node_x_ccoords - com_type

        if self._debug:
            self.ostream.print_info(
                f"center of mass type {self.com_target_type} at {com_type}")

            self.ostream.print_info(f"number of atoms: {len(node_atoms)}")
            self.ostream.print_info(
                f"number of X atoms: {len(node_x_ccoords)}")

        self.node_atoms = node_atoms
        self.node_ccoords = node_ccoords
        self.node_x_ccoords = node_x_ccoords

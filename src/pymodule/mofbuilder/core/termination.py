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
from typing import Any, Optional

import numpy as np
import networkx as nx
import mpi4py.MPI as MPI

from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...molecule import Molecule

from ..io.basic import nn
from ..io.pdb_reader import PdbReader
from ..io.pdb_writer import PdbWriter


class FrameTermination:
    """Loads and holds termination group geometry (e.g. acetate) for capping unsaturated nodes.

    Reads a PDB file containing X and Y atom types (e.g. X = connection atom,
    Y = O-O center), recenters to Y, and exposes termination_data and X/Y subsets.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        properties: Dictionary of optional properties.
        filename: Path to termination PDB file.
        X_atom_type: Atom type string for connection atom (default "X").
        Y_atom_type: Atom type string for O-O center (default "Y").
        pdbreader: PdbReader instance for reading PDB.
        _debug: If True, print extra debug messages.
        X_data: Subarray of termination_data for X atoms (set by create).
        termination_data: Full atom data from PDB (set by read_termination_file).
        termination_X_data: Rows of termination_data where last column equals X_atom_type (set by create).
        termination_Y_data: Rows of termination_data where last column equals Y_atom_type (set by create).
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
        self.X_atom_type = "X"
        self.Y_atom_type = "Y"
        self.pdbreader = PdbReader(comm=self.comm, ostream=self.ostream)

        self._debug = False
        self.X_data = None
        self.termination_data = None

    def read_termination_file(self) -> None:
        """Read the termination PDB from self.filename and set self.termination_data.

        Optionally prints debug info if _debug is True. Does nothing if filename is None.
        """
        if self.filename is None:
            return None
        assert_msg_critical(
            Path(self.filename).is_file(),
            f"Termination file {self.filename} does not exist.")
        if self._debug:
            self.ostream.print_info(
                f"Reading termination file {self.filename}")
        self.pdbreader.read_pdb(
            self.filename, recenter=True,
            com_type="Y")  # recenter to Y atom which is O-O center for XOO
        self.termination_data = self.pdbreader.data
        if self._debug:
            self.ostream.print_info(
                f"Got {len(self.termination_data)} atoms from termination file."
            )
            self.ostream.flush()

    def create(self) -> None:
        """Load termination file and split data into termination_X_data and termination_Y_data by atom type."""
        self.read_termination_file()
        if self.termination_data is not None:
            self.termination_X_data = self.termination_data[
                self.termination_data[:, -1] == self.X_atom_type]
            self.termination_Y_data = self.termination_data[
                self.termination_data[:, -1] == self.Y_atom_type]

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
import sys
from typing import Optional, Any, List, Sequence
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
import mpi4py.MPI as MPI
from ...errorhandler import assert_msg_critical


class XyzWriter:
    """Writer for XYZ coordinate files.

    Provides functionality to write molecular structure data to XYZ files,
    optionally with MPI parallel awareness and flexible stream/file handling.

    Attributes:
        comm (Any): MPI communicator used for parallel operations.
        rank (int): MPI rank of the current process.
        nodes (int): Total number of MPI nodes.
        ostream (OutputStream): Output stream for logging/info.
        filepath (Optional[str]): Default path to write the XYZ file if not provided per call.
        _debug (bool): If True, debug information is printed.
        file_dir (Optional[Path]): Directory of the current file to write (only set when writing).
    
    Methods:
        write(filepath, header, lines):
            Write the atom coordinates to an XYZ file.
        get_xyzlines(header, lines):
            Return formatted XYZ lines as a list of strings.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[OutputStream] = None,
        filepath: Optional[str] = None,
        debug: bool = False,
    ):
        """
        Initializes the XyzWriter instance.

        Args:
            comm (Optional[Any]): MPI communicator. Defaults to MPI.COMM_WORLD.
            ostream (Optional[OutputStream]): VeloxChem OutputStream for info/debug output.
            filepath (Optional[str]): Default path to the XYZ file.
            debug (bool): Enable debug printing if True.
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
        self.filepath = filepath
        self._debug = debug
        self.file_dir: Optional[Path] = None

    def write(
        self,
        filepath: Optional[str] = None,
        header: str = '',
        lines: Sequence[Sequence[Any]] = [],
    ) -> None:
        """Write atom coordinate lines to an XYZ file.

        Args:
            filepath (Optional[str]): Output XYZ file path. Uses instance default if not specified.
            header (str): Optional header/comment line (is line #2 in XYZ format).
            lines (Sequence[Sequence[Any]]): Each entry must have fields for
                atom_type, atom_label, atom_number, residue_name, residue_number,
                x, y, z, spin, charge, note.

                Only atom_label, x, y, z are directly written to the XYZ file.
        
        Raises:
            AssertionError: If filepath is not specified or invalid.

        Note:
            XYZ file will be written only by the master MPI rank.
            If the file extension is not '.xyz', it is automatically added.
        """
        filepath_final = Path(filepath) if filepath is not None else Path(self.filepath)
        assert_msg_critical(filepath_final is not None, "xyz filepath is not specified")

        # Check if the file directory exists and create it if it doesn't
        self.file_dir = Path(filepath_final).parent
        if self._debug:
            self.ostream.print_info(f"targeting directory: {self.file_dir}")
        self.file_dir.mkdir(parents=True, exist_ok=True)

        if filepath_final.suffix != ".xyz":
            filepath_final = filepath_final.with_suffix(".xyz")

        newxyz: List[str] = []
        newxyz.append(f"{len(lines)}\n")
        newxyz.append(header if header.endswith('\n') else header + '\n')

        for i, values in enumerate(lines):
            atom_label = values[1]
            # atom_number = i + 1  # Not used in XYZ format
            x = float(values[5])
            y = float(values[6])
            z = float(values[7])
            formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                atom_label, x, y, z
            )
            newxyz.append(formatted_line + "\n")

        with open(filepath_final, "w") as fp:
            fp.writelines(newxyz)

    def get_xyzlines(
        self,
        header: str = '',
        lines: Sequence[Sequence[Any]] = [],
    ) -> List[str]:
        """Return formatted XYZ contents as a list of strings.

        Args:
            header (str): Optional header/comment line for XYZ content.
            lines (Sequence[Sequence[Any]]): List of per-atom info, each entry
                having the fields:
                atom_type, atom_label, atom_number, residue_name, residue_number,
                x, y, z, spin, charge, note. Only atom_label, x, y, z are written.

        Returns:
            List[str]: Formatted XYZ contents as a list of strings.

        Example:
            xyz_lines = writer.get_xyzlines(header="Example molecule", lines=my_atomlist)
            # '\n'.join(xyz_lines) is a complete XYZ file or block.
        """
        newxyz: List[str] = []
        newxyz.append(f"{len(lines)}\n")
        newxyz.append(header.strip('\n') + '\n')
        for i, values in enumerate(lines):
            atom_label = values[1]
            # atom_number = i + 1  # Not used in XYZ format
            x = float(values[5])
            y = float(values[6])
            z = float(values[7])
            formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                atom_label, x, y, z
            )
            newxyz.append(formatted_line + "\n")
        return newxyz

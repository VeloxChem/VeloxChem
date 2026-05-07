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
from mpi4py import MPI
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from typing import Optional, Any, Dict
from ..io.xyz_writer import XyzWriter


class Viewer:
    """Visualization interface using py3Dmol for MOFbuilder.

    Attributes:
        comm (Any): MPI communicator (usually `MPI.COMM_WORLD`).
        rank (int): The MPI rank of this process.
        nodes (int): The total number of MPI ranks/nodes.
        ostream (OutputStream): VeloxChem-style output stream for logging/info.
        eG_dict (Optional[Dict[Any, Any]]): Dictionary mapping index to name (graph elements).
        merged_lines (Optional[Any]): Data lines to be visualized; formatted input for XYZWriter.
        eG_name_idx_dict (Optional[Dict[Any, int]]): Reverse-lookup of eG_dict, mapping names to indices.

    Methods:
        __init__: Initialize the Viewer with MPI context and output stream.
        _reverse_eG_dict: Reverse the element graph dictionary for label lookup.
        lines_show: Visualize merged_lines as a 3D molecular scene with labels.
    """

    def __init__(self, comm: Optional[Any] = None, ostream: Optional[Any] = None, filepath: Optional[str] = None):
        """Initialize the Viewer for MOF visualization.

        Args:
            comm (Optional[Any]): MPI communicator. Defaults to `MPI.COMM_WORLD`.
            ostream (Optional[Any]): Output stream for logging. Defaults to VeloxChem's OutputStream.
            filepath (Optional[str]): Path to data file (currently unused).
        """
        self.comm: Any = comm or MPI.COMM_WORLD
        self.rank: int = self.comm.Get_rank()
        self.nodes: int = self.comm.Get_size()

        if ostream is None:
            ostream = OutputStream(sys.stdout if self.rank == mpi_master() else None)
        self.ostream: Any = ostream
        self.eG_dict: Optional[Dict[Any, Any]] = None
        self.merged_lines: Optional[Any] = None
        self.eG_name_idx_dict: Optional[Dict[Any, int]] = None

    def _reverse_eG_dict(self) -> None:
        """Reverse the internally stored eG_dict for fast name-to-index lookup.

        Note:
            eG_dict should be a dictionary where keys are index numbers (usually int as str),
            and values are unique names of the corresponding graph elements.
            After this call, eG_name_idx_dict maps element names to index numbers.
        """
        # The eG dict is a dictionary: key = index, value = name.
        self.eG_name_idx_dict = {v: int(k) for k, v in self.eG_dict.items()}

    def lines_show(
        self,
        w: int = 800,
        h: int = 600,
        res_indices: bool = True,
        res_name: bool = True,
    ) -> None:
        """Display `merged_lines` in a 3D viewer with optional residue names and indices.

        Args:
            w (int): Width of viewer in pixels.
            h (int): Height of viewer in pixels.
            res_indices (bool): If True, show residue indices as labels.
            res_name (bool): If True, show residue names as labels.

        Raises:
            ImportError: If py3Dmol is not installed.

        Note:
            - Requires `py3Dmol` for visualization.
            - Uses `self.merged_lines`, typically a list of atom-info lines output from the builder.
            - Residue labels skip any lines with resname "TNO".
            - If both `res_name` and `res_indices` are True, the label concatenates both.
        """
        try:
            import py3Dmol
            merged_lines = self.merged_lines
            viewer = py3Dmol.view(width=w, height=h)
            xyz_writer = XyzWriter(comm=self.comm, ostream=self.ostream)
            xyz_lines = xyz_writer.get_xyzlines(lines=merged_lines)
            viewer.addModel("".join(xyz_lines), "xyz")
            viewer.setViewStyle({"style": "outline", "width": 0.05})
            viewer.setStyle({"stick": {}, "sphere": {"scale": 0.20}})
            if res_indices or res_name:
                self._reverse_eG_dict()
                old_resnumber = 0
                for line in merged_lines:
                    value_resname = line[3].split('_')[0][:3]
                    if value_resname.strip() == "TNO":
                        continue
                    value_resnumber = self.eG_name_idx_dict[line[10]]
                    if value_resnumber == old_resnumber:
                        continue

                    old_resnumber = value_resnumber

                    value_x = float(line[5])
                    value_y = float(line[6])
                    value_z = float(line[7])

                    text = ""
                    if res_name:
                        text += str(value_resname)
                    if res_indices:
                        text += str(value_resnumber)

                    viewer.addLabel(
                        text,
                        {
                            "position": {
                                "x": value_x,
                                "y": value_y,
                                "z": value_z,
                            },
                            "alignment": "center",
                            "fontColor": "white",
                            "font": "Arial",
                            "fontSize": 12,
                            "backgroundColor": "black",
                            "backgroundOpacity": 0.5,
                        },
                    )
            viewer.render()
            viewer.zoomTo()
            viewer.show()
        except ImportError:
            raise ImportError("Unable to import py3Dmol")

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
from typing import Any, Dict, List, Optional

import numpy as np
import networkx as nx
import mpi4py.MPI as MPI

from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from ...environment import get_data_path


class MofTopLibrary:
    """Lookup and management of MOF families, metals, and topology template CIFs.

    Reads the MOF_topology_dict file to get node connectivity, allowed metals,
    linker topic, and topology name. Supports listing families, selecting a family,
    and submitting new templates.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        data_path: Path to database directory (contains MOF_topology_dict).
        mof_top_dict: Dict mapping MOF family name to node_connectivity, metal list, linker_topic, topology.
        template_directory: Directory containing template CIF files.
        mof_family: Currently selected MOF family name.
        node_connectivity: Node connectivity for selected family.
        node_metal_type: Node metal type (set when used).
        linker_connectivity: Linker topic for selected family.
        net_type: Net type (set when used).
        net_filename: Template CIF filename for selected family.
        selected_template_cif_file: Full path to the selected template CIF (set by select_mof_family).
        _debug: If True, print extra debug messages.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[Any] = None,
        filepath: Optional[str] = None,
    ) -> None:
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # clean up
        if hasattr(self, "mof_top_dict"):
            del self.mof_top_dict
        if hasattr(self, "data_path"):
            del self.data_path

        self.data_path = get_data_path()
        self.mof_top_dict = None
        self.template_directory = None
        self.mof_family = None
        self.node_connectivity = None
        self.node_metal_type = None
        self.linker_connectivity = None
        self.net_type = None

        self._debug = False

    def _read_mof_top_dict(self, data_path: Optional[str] = None) -> None:
        """Load MOF_topology_dict from data_path and populate self.mof_top_dict.

        Args:
            data_path: Directory containing MOF_topology_dict. Uses self.data_path if None.
        """
        if data_path is None:
            data_path = self.data_path
        if Path(data_path, "MOF_topology_dict").exists():
            mof_top_dict_path = str(Path(data_path, "MOF_topology_dict"))
            with open(mof_top_dict_path, "r") as f:
                lines = f.readlines()
            # titles = lines[0].split()
            mofs = lines[1:]
        if self._debug:
            self.ostream.print_info(
                f"MOF_topology_dict path: {mof_top_dict_path}")
            self.ostream.print_info(
                f"Got {len(mofs)} MOF families from MOF_topology_dict")
            self.ostream.print_info(
                f"MOF families: {[mof.split()[0] for mof in mofs]}")
            self.ostream.flush()
        mof_top_dict = {}
        for mof in mofs:
            mof_name = mof.split()[0]
            if mof_name not in mof_top_dict.keys():
                mof_top_dict[mof_name] = {
                    "node_connectivity": int(mof.split()[1]),
                    "metal": [mof.split()[2]],
                    "linker_topic": int(mof.split()[3]),
                    "topology": mof.split()[-1],
                }
            else:
                mof_top_dict[mof_name]["metal"].append(mof.split()[2])
        self.mof_top_dict = mof_top_dict

    def list_mof_families(self) -> None:
        """Print all available MOF family names from the topology dictionary to the output stream."""
        if self.mof_top_dict is None:
            self._read_mof_top_dict(self.data_path)
        print("-" * 80)
        print("\t" * 3, "Available MOF Families:")
        print("-" * 80)
        for mof_family in self.mof_top_dict.keys():
            print(f" - {mof_family}")

    def list_available_metals(self, mof_family: str) -> None:
        """Print available metal types for the given MOF family; prints a warning if family is unknown.

        Args:
            mof_family: MOF family name (e.g. "HKUST-1").
        """
        mof_family = mof_family.upper()
        if self.mof_top_dict is None:
            self._read_mof_top_dict(self.data_path)
        if mof_family not in self.mof_top_dict.keys():
            self.ostream.print_warning(f"{mof_family} not in database")
            self.ostream.print_info("please select a MOF family from below:")
            self.ostream.flush()
            self.list_mof_families()
            return
        self.ostream.print_title(f"Available metals for {mof_family}:")
        for metal in self.mof_top_dict[mof_family.upper()]["metal"]:
            self.ostream.print_info(f" - {metal}")

        self.ostream.flush()

    def select_mof_family(self, mof_family: str) -> None:
        """Set the current MOF family and resolve template CIF path.

        Sets node_connectivity, linker_connectivity, net_filename, and
        selected_template_cif_file. Existence of the template CIF is checked.

        Args:
            mof_family: MOF family name (e.g. "UIO-66").
        """
        self.mof_family = mof_family.upper()
        self.node_connectivity = self.mof_top_dict[mof_family][
            "node_connectivity"]
        self.linker_connectivity = self.mof_top_dict[mof_family][
            "linker_topic"]
        self.net_filename = self.mof_top_dict[mof_family]["topology"] + ".cif"
        # check if template cif exists
        self.ostream.print_info(f"MOF family {mof_family} is selected")
        if self.mof_family not in self.mof_top_dict.keys():
            self.ostream.print_warning(f"{mof_family} not in database")
            self.ostream.print_info("please select a MOF family from below:")
            self.ostream.flush()
            self.list_mof_families()
            return
        self.ostream.print_info(f"node connectivity: {self.node_connectivity}")
        self.ostream.print_info(f"linker topic: {self.linker_connectivity}")
        self.ostream.print_info(
            f"available metal nodes: {self.mof_top_dict[self.mof_family]['metal']}"
        )
        self.ostream.flush()
        if self.template_directory is None:
            self.template_directory = str(
                Path(self.data_path, "template_database"))  # default
            self.ostream.print_info(
                f"Searching template cif files in {self.template_directory}..."
            )
            self.ostream.flush()

        template_cif_file = str(
            Path(self.template_directory, self.net_filename))

        if not Path(template_cif_file).exists():
            self.ostream.print_info(
                f"{self.net_filename} net does not exist in {self.template_directory}"
            )
            self.ostream.print_info(
                "please select another MOF family, or upload the template cif file"
            )

            #TODO: set it as repository for template cif files
            self.ostream.print_info(
                "or download the template cif files from the internet and  set it as the template directory"
            )
            self.ostream.flush()
            return
        else:
            self.ostream.print_info(
                f"{self.net_filename} is found and will be used for MOF building"
            )
            self.ostream.flush()
            self.selected_template_cif_file = template_cif_file

    def submit_template(
        self,
        template_cif: str,
        mof_family: str,
        template_mof_node_connectivity: int,
        template_node_metal: str,
        template_linker_topic: int,
        overwrite: bool = False,
    ) -> Optional[str]:
        """Add or overwrite a MOF family in the topology dict and template_database.

        Validates template_cif path and extension, then updates mof_top_dict and
        rewrites the MOF_topology_dict file. Optionally overwrites existing family.

        Args:
            template_cif: Path to the template CIF file.
            mof_family: MOF family name (e.g. "HKUST-1").
            template_mof_node_connectivity: Node connectivity (int).
            template_node_metal: Metal symbol (str).
            template_linker_topic: Linker topic (int).
            overwrite: If True, replace existing entry for mof_family.

        Returns:
            str: Path to the template CIF in template_database, or None if not written.
        """
        mof_family = mof_family.upper()
        assert_msg_critical(
            Path(self.data_path, "template_database").is_dir(),
            f"template_database directory {Path(self.data_path, 'template_database')} does not exist, please create it first",
        )

        assert_msg_critical(
            Path(template_cif).exists(),
            f"template cif file {template_cif} does not exist, please upload it first",
        )

        assert_msg_critical(
            Path(template_cif).suffix == ".cif",
            f"template cif file {template_cif} is not a cif file, please upload a cif file",
        )

        assert isinstance(template_mof_node_connectivity,
                          int), "please enter an integer for node connectivity"
        assert isinstance(template_node_metal,
                          str), "please enter a string for node metal"
        assert isinstance(template_linker_topic,
                          int), "please enter an integer for linker topic"

        self.ostream.print_info(
            f"Submitting {mof_family} to the database {Path(self.data_path, 'template_database')}"
        )
        self.ostream.print_info(f"template cif file: {template_cif}")
        self.ostream.print_info(
            f"node connectivity: {template_mof_node_connectivity}")
        self.ostream.print_info(f"node metal: {template_node_metal}")
        self.ostream.print_info(f"linker topic: {template_linker_topic}")
        self.ostream.print_info(f"overwrite existing mof family: {overwrite}")
        self.ostream.flush()

        if mof_family in self.mof_top_dict.keys():
            if not overwrite:
                self.ostream.print_warning(
                    f"{mof_family} already exists in the database, the template you submitted will not be used, or you can set overwrite=True to overwrite the existing template"
                )
                return None

        self.mof_top_dict[mof_family] = {
            "node_connectivity": template_mof_node_connectivity,
            "metal": [template_node_metal],
            "linker_topic": template_linker_topic,
            "topology": Path(template_cif).stem,
        }
        self.ostream.print_info(f"{mof_family} is added to the database")
        self.ostream.print_info(f"{mof_family} is ready for MOF building")
        self.ostream.flush()

        # rewrite mof_top_dict file
        with open(str(Path(self.data_path, "MOF_topology_dict")), "w") as fp:
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
        self.ostream.print_info("mof_top_dict file is updated")
        self.ostream.flush()
        return str(Path(self.data_path, "template_database", template_cif))

    def fetch(self, mof_family: Optional[str] = None) -> Optional[str]:
        """Load topology dict, select the given MOF family, and return its template CIF path.

        Args:
            mof_family: MOF family name (e.g. "UIO-66"). If None, lists families and returns None.

        Returns:
            Path to selected_template_cif_file if family is valid; None otherwise.
        """
        mof_family = (mof_family or "").upper()
        self._read_mof_top_dict(self.data_path)
        if not mof_family:
            self.ostream.print_info("please select a MOF family from below:")
            self.ostream.flush()
            self.list_mof_families()
            return None
        else:
            if mof_family not in self.mof_top_dict.keys():
                self.ostream.print_warning(f"{mof_family} not in database")
                self.ostream.print_info(
                    "please select a MOF family from below:")
                self.ostream.flush()
                self.list_mof_families()
                return None
            else:
                self.select_mof_family(mof_family)
                return self.selected_template_cif_file

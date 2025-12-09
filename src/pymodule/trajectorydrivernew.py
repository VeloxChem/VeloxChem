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
from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
import numpy as np
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
import MDAnalysis as mda
import MDAnalysis.transformations as transform

class TrajectoryDriver:
    """
    Driver for trajectory-driven calculations.
    
    This class automates the parsing of trajecoties.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.

    Instance variables:
        - qm_molecule: The QM molecule.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the TrajectoryDriver.
        """
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
        self.size = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def trajectory_parser(self,
                          topology_file: str,
                          trajectory_file: str,
                          num_snapshots: int,
                          qm_region: str):
        """
        Trajectory parser for managing snapshots.
        It returns a list of dictionaries: 
        1 dictionary per snapshot
        key: QM region coordinates
        key: MM region coordinates
        """
        self.universe = mda.Universe(topology_file, trajectory_file)

        total_frames = len(self.universe.trajectory)
        print(f"Total frames in trajectory: {total_frames}")
        if num_snapshots > total_frames:
            raise ValueError(f"Requested number of snapshots ({num_snapshots}) exceeds total frames in trajectory ({total_frames}).")
        
        if num_snapshots == 1:
            start = 0
            stop = 1
            step = 1
        else:
            start = 0
            stop = total_frames
            step = (total_frames -1) // (num_snapshots -1)

        self.start = start
        self.stop = stop
        self.step = step
        #import MDAnalysis.transformations as transform
        qm_atoms = self.universe.select_atoms(qm_region)
        rest = self.universe.select_atoms(f"not ({qm_region})")
        transforms = [
            transform.unwrap(qm_atoms),
            transform.center_in_box(qm_atoms, wrap=True),
            transform.wrap(rest)
        ]
        self.universe.trajectory.add_transformations(*transforms)
        snapshots = []
        for ts in self.universe.trajectory[self.start:self.stop:self.step]:
            qm_coords = qm_atoms.positions / bohr_in_angstrom()
            mm_coords = rest.positions / bohr_in_angstrom()
            snapshot = {
                'frame': ts.frame,
                'qm_coords': qm_coords,
                'mm_coords': mm_coords
            }
            snapshots.append(snapshot)
        return snapshots

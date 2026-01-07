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
    
    This class automates the parsing of molecular dynamics trajectories and
    extracts QM and MM region coordinates for each snapshot.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.

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
                          trajectory_file: str,
                          num_snapshots: int,
                          qm_region: str,
                          topology_file: str = None):
        """
        Parse a molecular dynamics trajectory and extract QM and MM region coordinates.

        Args:
            trajectory_file (str): Path to the trajectory file (e.g., .dcd, .xtc, .trr, etc.)
            num_snapshots (int): Number of snapshots to extract from the trajectory.
                Must be <= total number of frames in trajectory.
            qm_region (str): MDAnalysis selection string defining the QM region atoms.
                Example: 'resname MOL'
            topology_file (str, optional): Path to the topology file (.tpr, etc)
        
        Returns:
            List[Dict]: A list of snapshots, each represented as a dictionary with keys:
                'frame' (int): Frame number from trajectory
                'qm_coords' (np.ndarray): Coordinates of QM region atoms in Bohr
                'mm_coords' (np.ndarray): Coordinates of MM region atoms in Bohr

        Raises:
            ValueError: If num_snapshots exceeds total frames in trajectory.
            FileNotFoundError: If trajectory or topology files are not found.
            
        Example:
            >>> trj_drv = TrajectoryDriver()
            >>> snapshots = trj_drv.trajectory_parser(
                trajectory_file='traj.xtc',
                num_snapshots=4,
                qm_region='resname MOL',
                topology_file='topology.tpr'
            )
            >>> print(snapshots[0].keys())
            dict_keys(['frame', 'qm_coords', 'mm_coords'])
        """
        if trajectory_file.lower().endswith('.pdb'):
            self.universe = mda.Universe(trajectory_file, guess_bonds=True)
        else:
            self.universe = mda.Universe(topology_file, trajectory_file)

        total_frames = len(self.universe.trajectory)
        print(self.universe.trajectory.units)
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
            qm_coords = qm_atoms.positions
            mm_coords = rest.positions
            snapshot = {
                'frame': ts.frame,
                'qm_coords': qm_coords,
                'mm_coords': mm_coords,
                'mm_atom_names': rest.atoms.names,
                'mm_resids': rest.atoms.resids,
                'mm_resnames': rest.atoms.resnames
            }
            snapshots.append(snapshot)
        return snapshots

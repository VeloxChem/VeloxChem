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
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
import numpy as np


class EnvironmentDriver:
    """
    Handles the evironment from the snaphosts of a trajectory.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the EnvironmentDriver.
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

    def write_pot_files(self, snapshots, outdir: str | Path):
        """
        Write environment snapshots to .pot files.

        Generates one .pot file per snapshot.

        Returns:
            None

        Raises:
            KeyError: If required keys are missing from snapshot dictionaries.
            PermissionError: If output directory cannot be created or written to.

        Example:
            >>> env_drv = EnvironmentDriver()
            >>> snapshots = trj_drv.trajectory_parser(...)
            >>> env_drv.write_pot_files(snapshots, outdir='pot_frames')
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for snap in snapshots:
            self._write_single_pot(snap, outdir)

    def _write_single_pot(self, snapshot, outdir: Path):
        """
        Write a single snapshot to a .pot file.

        Generates a .pot file for one snapshot containing the @environment,
        @charges, and @polarizabilities sections formatted for VeloxChem.
        
        Args:
            snapshot (Dict): A single snapshot dictionary containing:
                - 'frame' (int): Frame number
                - 'mm_coords' (np.ndarray): MM region coordinates in Angstrom (N, 3)
                - 'mm_elements' (np.ndarray): Element symbols for each atom (e.g., 'O', 'H')
                - 'mm_resids' (np.ndarray): Residue IDs for each atom
                - 'mm_resnames' (np.ndarray): Residue names (e.g., 'SOL')
            outdir (Path): Output directory path.

        Returns:
            None

        Raises:
            IOError: If .pot file cannot be written.

        """
        frame = snapshot['frame']
        mm_coords = snapshot['mm_coords']
        mm_elements = snapshot['mm_elements']
        mm_resids = snapshot['mm_resids']
        mm_resnames = snapshot['mm_resnames']

        pot_path = outdir / f"frame_{frame:06d}.pot"
        charges = {'O': -0.67444000, 'H': 0.33722000}
        polar = {
            'O': [0.0, 0.0, 5.73935000, 0.0, 0.0, 0.0],
            'H': [0.0, 0.0, 2.30839000, 0.0, 0.0, 0.0],
        }

        with pot_path.open('w') as fh:
            fh.write("@environment\n")
            fh.write("units: angstrom\n")
            fh.write("xyz:\n")
            for (x, y, z), elem, resn, resid in zip(mm_coords, mm_elements, mm_resnames, mm_resids):
                fh.write(f"{elem:<2} {x:12.6f} {y:12.6f} {z:12.6f}  {resn:>3}  {resid}\n")
            fh.write("@end\n\n")

            fh.write("@charges\n")
            for atom in ['O', 'H', 'H']:
                fh.write(f"{atom:<2} {charges[atom]:12.8f}  SOL\n")
            fh.write("@end\n\n")

            fh.write("@polarizabilities\n")
            for atom in ['O', 'H', 'H']:
                vals = polar[atom]
                fh.write(
                    f"{atom:<2} {vals[0]:12.8f} {vals[1]:12.8f} {vals[2]:12.8f} "
                    f"{vals[3]:12.8f} {vals[4]:12.8f} {vals[5]:12.8f}  SOL\n"
                )
            fh.write("@end\n")


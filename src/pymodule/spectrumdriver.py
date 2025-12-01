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
import pyframe

class SpectrumDriver:
    """
    Driver for the sppectra in polarizable embedding (PE) calculations.
    
    This class automates the creation of PE potentials from solvated systems
    and executes PE calculations.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.

    Instance variables:
        - qm_molecule: The QM molecule.
        - solvated system: The complete solvates system from SolvationBuilder.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the SpectrumDriver.
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

    @staticmethod
    def trajectory_parser(topology_file, trajectory_file, 
                          qm_selection, filename_core, 
                          filename_environment, snapshots):
        """"
        Trajectory parser for extracting snapshots.
        """
        
        traj = pyframe.Trajectory(topology_file, 
                                  trajectory_file,
                                  qm_selection=qm_selection,
                                  snapshots=snapshots)
        traj.set_core_region(qm_selection)
        traj.add_region(
            name='test',
            selection='resname SOL or resname Na+',
            use_standard_potentials=True,
            standard_potential_model='SEP'
        )
        traj.create_potential()
        traj.write_core(filename=filename_core)
        traj.write_potential(filename=filename_environment)


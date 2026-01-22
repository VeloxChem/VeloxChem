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

from mpi4py import MPI
import sys
import numpy as np
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
import MDAnalysis as mda
import MDAnalysis.transformations as transform
from MDAnalysis.topology.guessers import guess_atom_element

class EnsembleParser:
    """
    
    This class automates the parsing of molecular dynamics trajectories and
    extracts information from the qm and environment regions for each snapshot.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.
    """
    def __init__(self, comm=None, ostream=None):
        """
        Initialize the EnsembleParser.
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

    def trajectory(self,
                   trajectory_file: str,
                   num_snapshots: int,
                   qm_region: str,
                   env_region: str,
                   topology_file: str | None = None,
                   pe_cutoff: float | None = None,
                   npe_cutoff: float | None = None,):
        """
        Parse a molecular dynamics trajectory and extract QM and MM region data.

        :param trajectory_file:
            Path to the trajectory file (e.g., .dcd, .xtc, .pdb).
        :param topology_file:
            Path to the topology file (e.g., .tpr).
        :param num_snapshots:
            Number of snapshots to extract from the trajectory.
        :param qm_region:
            Selection string defining the QM region (MDAnalysis selection syntax).
        :param env_region:
            Selection string defining the environment region (MDAnalysis selection syntax).
        :param pe_cutoff:
            Cutoff for polarizable embedding (PE) region selection (Angstrom).
        :param npe_cutoff:
            Cutoff for non-polarizable embedding (NPE) region selection (Angstrom).
        :return:
            A list of snapshot dictionaries, each containing:
            - frame (int):
                Frame number.
            - qm_coords (numpy.ndarray):
                QM region Cartesian coordinates, shape (N_qm, 3), in Angstrom.
            - qm_elements (numpy.ndarray):
                Element symbols for each QM atom, shape (N_qm,).
            - pe_coords (numpy.ndarray):
                PE region Cartesian coordinates, shape (N_pe, 3), in Angstrom.
            - pe_elements (numpy.ndarray):
                Element symbols for each PE atom, shape (N_pe,).
            - pe_resids (numpy.ndarray):
                Residue id for each PE atom, shape (N_pe,).
            - pe_resnames (numpy.ndarray):
                Residue name for each PE atom, shape (N_pe,).
            - pe_n_residues (int):
                Number of residues in the PE region.
            - npe_coords (numpy.ndarray):
                NPE region Cartesian coordinates, shape (N_npe, 3), in Angstrom.
            - npe_elements (numpy.ndarray):
                Element symbols for each NPE atom, shape (N_npe,).
            - npe_resids (numpy.ndarray):
                Residue id for each NPE atom, shape (N_npe,).
            - npe_resnames (numpy.ndarray):
                Residue name for each NPE atom, shape (N_npe,).
            - npe_n_residues (int):
                Number of residues in the NPE region.
        """
        
        if num_snapshots <= 0:
            raise ValueError("num_snapshots must be a positive integer")
        
        if pe_cutoff is not None and npe_cutoff is not None:
            if float(npe_cutoff) < float(pe_cutoff):
                raise ValueError("npe_cutoff must be >= pe_cutoff")
        
        if trajectory_file.lower().endswith('.pdb'):
            self.universe = mda.Universe(trajectory_file, guess_bonds=True)
        else:
            if topology_file is None:
                raise ValueError("topology_file is required unless trajectory_file is a .pdb")
            self.universe = mda.Universe(topology_file, trajectory_file)
        
        total_frames = len(self.universe.trajectory)
        self.ostream.print_info(f"Total frames in trajectory: {total_frames}")
        self.ostream.print_blank()
        if num_snapshots > total_frames:
            raise ValueError(
                f"Requested number of snapshots ({num_snapshots}) exceeds total frames ({total_frames})."
            )
    
        if num_snapshots == 1:
            frame_indices = np.array([0], dtype=int)
        else:
            frame_indices = np.linspace(0, total_frames - 1, num_snapshots, dtype=int)
        
        qm_atoms = self.universe.select_atoms(qm_region)
        env_atoms = self.universe.select_atoms(env_region)

        transforms = [
        transform.unwrap(qm_atoms),
        transform.center_in_box(qm_atoms, wrap=True),
        transform.wrap(env_atoms)
        ]
        self.universe.trajectory.add_transformations(*transforms)

        empty_xyz = np.empty((0, 3), dtype=float)
        empty_obj = np.empty((0,), dtype=object)
        empty_int = np.empty((0,), dtype=int)

        snapshots = []
        for iframe in frame_indices:
            self.universe.trajectory[iframe]

            qm_coords = np.asarray(qm_atoms.positions, dtype=float).copy()
            qm_elements = np.asarray(
                [guess_atom_element(n) for n in qm_atoms.names], dtype=object
                )
            
            # Defaults
            pe_coords, pe_elements, pe_resids, pe_resnames, pe_n_residues = (
                empty_xyz, empty_obj, empty_int, empty_obj, 0
            )
            npe_coords, npe_elements, npe_resids, npe_resnames, npe_n_residues = (
                empty_xyz, empty_obj, empty_int, empty_obj, 0
            )

            pe_region = None

            # If neither cutoff is set, interpret as all-NPE environment
            if pe_cutoff is None and npe_cutoff is None:
                npe_region = env_atoms.difference(qm_atoms)
                npe_coords = np.asarray(npe_region.positions, dtype=float).copy()
                npe_elements = np.asarray(
                    [guess_atom_element(n) for n in npe_region.names], dtype=object
                )
                npe_resids = np.asarray(npe_region.resids, dtype=int).copy()
                npe_resnames = np.asarray(npe_region.resnames, dtype=object).copy()
                npe_n_residues = int(npe_region.residues.n_residues)
            
            # PE region
            if pe_cutoff is not None:
                pe_region = self.universe.select_atoms(
                    f"byres ({env_region} and around {float(pe_cutoff)} group qm)",
                    qm=qm_atoms,
                ).difference(qm_atoms)

                pe_coords = np.asarray(pe_region.positions, dtype=float).copy()
                pe_elements = np.asarray(
                    [guess_atom_element(n) for n in pe_region.names], dtype=object
                )
                pe_resids = np.asarray(pe_region.resids, dtype=int).copy()
                pe_resnames = np.asarray(pe_region.resnames, dtype=object).copy()
                pe_n_residues = int(pe_region.residues.n_residues)

            # NPE region
            if npe_cutoff is not None:
                outer_shell = self.universe.select_atoms(
                    f"byres ({env_region} and around {float(npe_cutoff)} group qm)",
                    qm=qm_atoms,
                ).difference(qm_atoms)

                npe_region = outer_shell.difference(pe_region) if pe_region is not None else outer_shell

                npe_coords = np.asarray(npe_region.positions, dtype=float).copy()
                npe_elements = np.asarray(
                    [guess_atom_element(n) for n in npe_region.names], dtype=object
                )
                npe_resids = np.asarray(npe_region.resids, dtype=int).copy()
                npe_resnames = np.asarray(npe_region.resnames, dtype=object).copy()
                npe_n_residues = int(npe_region.residues.n_residues)

            snapshot = {
                    "frame": int(self.universe.trajectory.frame),

                    "qm_coords": qm_coords,
                    "qm_elements": qm_elements,

                    "pe_coords": pe_coords,
                    "pe_elements": pe_elements,
                    "pe_resids": pe_resids,
                    "pe_resnames": pe_resnames,
                    "pe_n_residues": pe_n_residues,

                    "npe_coords": npe_coords,
                    "npe_elements": npe_elements,
                    "npe_resids": npe_resids,
                    "npe_resnames": npe_resnames,
                    "npe_n_residues": npe_n_residues,
            }
            snapshots.append(snapshot)
        return snapshots
    
    def structures(self,
                   structures_file: str):
        """
        Parse a PDB or XYZ that contains different structures

        :param structures_file:
            Path to the structures file (e.g., .pdb, .xyz).
 
        :return:
            A list of snapshot dictionaries, each containing:
            - frame (int):
                Frame number.
            - qm_coords (numpy.ndarray):
                QM region Cartesian coordinates, shape (N_qm, 3), in Angstrom.
            - qm_elements (numpy.ndarray):
                Element symbols for each QM atom, shape (N_qm,).
        """
        self.universe = mda.Universe(structures_file, guess_bonds=True)

        snapshots = []
        for iframe in self.universe.trajectory:
            coords = np.asarray(self.universe.atoms.positions, dtype=float).copy()
            elements = np.asarray(
                [guess_atom_element(n) for n in self.universe.atoms.names], dtype=object
            )

            snapshot = {
                    "frame": int(iframe.frame),
                    "qm_coords": coords,
                    "qm_elements": elements,

            }
            snapshots.append(snapshot)
        return snapshots

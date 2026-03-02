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

    @staticmethod
    def _prefixed_resname(resname: str, prefix: str) -> str:
        """
        Resname in PDBs often go by 3 letters, and terminal residues (e.g., NASN or CASN for
        n-terminal or c-terminal ASN) contain atom parameters slightly different from 
        the non-terminal residues.
        This routine adds a terminal prefix to a residue name, 
        so that the correct parameters can be later assigned.

        Parameters
        ----------
        resname : str
            Residue name as read from the structure/topology (e.g. 'ASN').
        prefix : str
            Terminal prefix to add. Expected values are 'N' or 'C'.

        Returns
        -------
        str
            Prefixed residue name (e.g. 'NASN', 'CASN'), or the original name
            if it already looks prefixed.
        """
        # If the residue name already appears to be terminal-prefixed 
        # (four characters starting with 'N' or 'C'),
        # it is returned unchanged.
        if isinstance(resname, str) and len(resname) == 4 and resname[0] in ("N", "C"):
            return resname
        return f"{prefix}{resname}"
    
    def _terminal_resname_map(self, env_atoms):
        """
        Identifies N- and C-terminal protein residues and returns a renaming map.

        For each protein chain found in the selection, the first residue is 
        considered the N-terminus and the last residue the C-terminus. 
        These residues are assigned terminal-prefixed 
        residue names (e.g. 'ASN' -> 'NASN' or 'CASN').

        The returned mapping is keyed by MDAnalysis `resindex` (unique per
        residue in the Universe) and can be used to update per-atom `resnames`
        arrays (e.g. `AtomGroup.resnames`) based on `AtomGroup.resindices`.

        Parameters
        ----------
        env_atoms : MDAnalysis.core.groups.AtomGroup
            AtomGroup corresponding to the environment selection used in
            trajectory parsing (e.g. `universe.select_atoms(env_region)`).

        Returns
        -------
        dict[int, str]
            Dictionary mapping residue `resindex` to the renamed residue name
            with terminal prefix. Only protein terminal residues are included.
            Returns an empty dictionary if no protein is present in `env_atoms`.
        """
        # Restrict to protein residues
        # prot = env_atoms.select_atoms("protein")
        prot = self.universe.select_atoms("protein")
        if len(prot) == 0:
            return {}

        # Identify chains
        chains = []
        try:
            if len(prot.bonds) > 0:
                frags = [f for f in prot.fragments if len(f.residues) > 1]
                if frags:
                    chains = frags
        except Exception:
            chains = []

        # Fallback: group by chainIDs (PDB) or segids
        if not chains:
            ids = getattr(prot.atoms, "chainIDs", None)
            if ids is None:
                ids = prot.atoms.segids
            ids = np.asarray(ids, dtype=object)

            for uid in np.unique(ids):
                ag = prot[ids == uid]
                if len(ag.residues) > 0:
                    chains.append(ag)

        term_map = {}
        for ch in chains:
            res = ch.residues
            if len(res) == 0:
                continue

            nterm = res[0]
            cterm = res[-1]

            term_map[nterm.resindex] = self._prefixed_resname(nterm.resname, "N")
            term_map[cterm.resindex] = self._prefixed_resname(cterm.resname, "C")

        return term_map

    def structures(self,
                   trajectory_file: str,
                   num_snapshots: int | None = None,
                   qm_region: str = "",
                   env_region: str | None = None,
                   topology_file: str | None = None,
                   pe_cutoff: float | None = None,
                   npe_cutoff: float | None = None,
    ):
        """
        Parse a set of structures and extract QM and MM region data.

        :param trajectory_file:
            Path to the trajectory file:
                - .xtc with a corresponding topology (.tpr) via topology_file
                - .pdb (several configurations, or a single configuration; bonds are guessed        
        :param topology_file:
            Path to the topology file (e.g., .tpr).
        :param num_snapshots:
            Number of snapshots to extract. If None (default), all frames are
            processed. If an integer, frames are sampled evenly.
        :param qm_region:
            Selection string defining the QM region (MDAnalysis selection syntax).
        :param env_region:
            Selection string defining the environment region (MDAnalysis selection syntax).
            If not provided, the environment is taken as all atoms not contained in
            the QM region.
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
            - pe_atom_names (numpy.ndarray):
                Atom names for each PE atom, shape (N_pe,).
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
            - npe_atom_names (numpy.ndarray):
                Atom names for each NPE atom, shape (N_npe,).
            - npe_elements (numpy.ndarray):
                Element symbols for each NPE atom, shape (N_npe,).
            - npe_resids (numpy.ndarray):
                Residue id for each NPE atom, shape (N_npe,).
            - npe_resnames (numpy.ndarray):
                Residue name for each NPE atom, shape (N_npe,).
            - npe_n_residues (int):
                Number of residues in the NPE region.
        """
               
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
        self.ostream.print_blank()

        if num_snapshots is None:
            num_snapshots = total_frames

        if num_snapshots <= 0:
            raise ValueError("num_snapshots must be a positive integer")
        
        if num_snapshots > total_frames:
            raise ValueError(
                f"Requested number of snapshots ({num_snapshots}) exceeds total frames ({total_frames})."
            )
        
        if num_snapshots == total_frames:
            frame_indices = np.arange(total_frames, dtype=int)
        elif num_snapshots == 1:
            frame_indices = np.array([0], dtype=int)
        else:
            frame_indices = np.linspace(0, total_frames - 1, num_snapshots, dtype=int)
           
        qm_atoms = self.universe.select_atoms(qm_region)
        if len(qm_atoms) == 0:
            raise ValueError(f"QM region '{qm_region}' selection is empty")
        
        if env_region is None:
            env_region_sel = f"not ({qm_region})"
        else:
            env_region_sel = str(env_region)

        env_atoms = self.universe.select_atoms(env_region_sel).difference(qm_atoms)

        # Identify terminal protein residues (if any) and assign N*/C* residue names
        # so that terminal variants can be treated as separate residue types downstream.
        term_map = self._terminal_resname_map(env_atoms)

        # MDAnalysis transforms require valid box dimensions (ts.dimensions)
        # (e.g. unwrap/center_in_box/wrap). For single PDBs without box info,
        # skip transformations.
        has_box = False
        try:
            self.universe.trajectory[0]
            dims = getattr(self.universe.trajectory.ts, 'dimensions', None)
            if dims is not None:
                dims = np.asarray(dims, dtype=float)
                if dims.size > 3 and np.all(dims[:3] > 0):
                    has_box = True
        except Exception:
            has_box = False

        if has_box:
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
            
            pe_coords = empty_xyz
            pe_elements = empty_obj
            pe_resids = empty_int
            pe_resnames = empty_obj
            pe_n_residues = 0
            pe_atom_names = empty_obj

            npe_coords = empty_xyz
            npe_elements = empty_obj
            npe_resids = empty_int
            npe_resnames = empty_obj
            npe_n_residues = 0
            npe_atom_names = empty_obj

            pe_region = None

            # PE selection
            if pe_cutoff is not None:
                pe_region = self.universe.select_atoms(
                    f"byres ({env_region_sel} and around {float(pe_cutoff)} group qm)",
                    qm=qm_atoms,
                ).difference(qm_atoms)

                pe_coords = np.asarray(pe_region.positions, dtype=float).copy()
                pe_atom_names = np.asarray(pe_region.names, dtype=object).copy()
                pe_elements = np.asarray(
                    [guess_atom_element(n) for n in pe_region.names], dtype=object
                )
                pe_resids = np.asarray(pe_region.resids, dtype=int).copy()
                pe_resnames = np.asarray(pe_region.resnames, dtype=object).copy()
                
                # Apply terminal residue renaming (NASN/CASN, etc.) if applicable
                if term_map and len(pe_resnames) > 0:
                    pe_residx = np.asarray(pe_region.resindices, dtype=int)
                    for ridx, newname in term_map.items():
                        pe_resnames[pe_residx == ridx] = newname

                pe_n_residues = int(pe_region.residues.n_residues)

            # NPE selection
            if npe_cutoff is not None:
                outer_shell = self.universe.select_atoms(
                    f"byres ({env_region_sel} and around {float(npe_cutoff)} group qm)",
                    qm=qm_atoms,
                ).difference(qm_atoms)

                npe_region = outer_shell.difference(pe_region) if pe_region is not None else outer_shell

                npe_coords = np.asarray(npe_region.positions, dtype=float).copy()
                npe_atom_names = np.asarray(npe_region.names, dtype=object).copy()
                npe_elements = np.asarray(
                    [guess_atom_element(n) for n in npe_region.names], dtype=object
                )
                npe_resids = np.asarray(npe_region.resids, dtype=int).copy()
                npe_resnames = np.asarray(npe_region.resnames, dtype=object).copy()

                # Apply terminal residue renaming (NASN/CASN, etc.) if applicable
                if term_map and len(npe_resnames) > 0:
                    npe_residx = np.asarray(npe_region.resindices, dtype=int)
                    for ridx, newname in term_map.items():
                        npe_resnames[npe_residx == ridx] = newname

                npe_n_residues = int(npe_region.residues.n_residues)                

            # If neither cutoff is set, interpret as all-qm

            snapshot = {
                    "frame": int(self.universe.trajectory.frame),

                    "qm_coords": qm_coords,
                    "qm_elements": qm_elements,

                    "pe_coords": pe_coords,
                    "pe_atom_names": pe_atom_names,
                    "pe_elements": pe_elements,
                    "pe_resids": pe_resids,
                    "pe_resnames": pe_resnames,
                    "pe_n_residues": pe_n_residues,

                    "npe_coords": npe_coords,
                    "npe_atom_names": npe_atom_names,
                    "npe_elements": npe_elements,
                    "npe_resids": npe_resids,
                    "npe_resnames": npe_resnames,
                    "npe_n_residues": npe_n_residues,
            }
            snapshots.append(snapshot)
        return snapshots

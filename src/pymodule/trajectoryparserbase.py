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

from importlib.metadata import version
import re
import sys

from mpi4py import MPI
import numpy as np

from .veloxchemlib import mpi_master, chemical_element_identifier
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical, print_exception_if_debug

try:
    import MDAnalysis as mda
    from MDAnalysis.guesser.default_guesser import DefaultGuesser as mda_DefaultGuesser
    from MDAnalysis.lib.mdamath import make_whole
except ImportError:
    pass


class TrajectoryParserBase:
    """
    Internal shared implementation for parsing PDB structures and
    molecular-dynamics trajectories with MDAnalysis.

    This class is an implementation detail. Use :class:`EmbeddingParser` for
    a single structure or :class:`EnsembleParser` for a trajectory/ensemble.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the trajectory parser base.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
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
    def _pdb_has_conect_records(pdb_file: str) -> bool:
        """
        Checks whether a PDB file contains CONECT records.

        :param pdb_file:
            Path to the PDB file.
        :return:
            True if at least one CONECT record is present.
        """
        with open(pdb_file, "rb") as fh:
            for line in fh:
                if line[:6].strip().upper() == b"CONECT":
                    return True
        return False

    @staticmethod
    def _prefixed_resname(resname: str, prefix: str) -> str:
        """
        Resname in PDBs often go by 3 letters, and terminal residues (e.g., NASN or CASN for
        n-terminal or c-terminal ASN) contain atom parameters slightly different from
        the non-terminal residues.
        This routine adds a terminal prefix to a residue name,
        so that the correct parameters can be later assigned.

        :param resname (str)
            Residue name as read from the structure/topology (e.g. 'ASN').
        :param prefix (str)
            Terminal prefix to add. Expected values are 'N' or 'C'.
        :return (str)
            Prefixed residue name (e.g. 'NASN', 'CASN'), or the original name
            if it already looks prefixed.
        """
        # If the residue name already appears to be terminal-prefixed
        # (four characters starting with 'N' or 'C'),
        # it is returned unchanged.
        if isinstance(resname, str) and len(resname) == 4 and resname[0] in ("N", "C"):
            return resname
        return f"{prefix}{resname}"

    @staticmethod
    def _histidine_resname_from_atoms(resname: str, atom_names: set[str]) -> str | None:
        """
        Infer CHARMM-style histidine protonation-state residue names from atom names.

        Generic HIS is ambiguous in parameter databases. The presence of HD1
        identifies delta-protonated histidine (HSD/HID), HE2 identifies
        epsilon-protonated histidine (HSE/HIE), and both hydrogens identify
        doubly protonated histidine (HSP/HIP).

        :param resname:
            Residue name as read from the structure/topology.
        :param atom_names:
            Set of atom names in the residue.
        :return:
            Refined CHARMM-style histidine residue name, or None if no
            histidine refinement is needed.
        """
        resname = str(resname)
        prefix = ""
        core_resname = resname

        if len(resname) == 4 and resname[0] in ("N", "C") and resname[1:] == "HIS":
            prefix = resname[0]
            core_resname = resname[1:]

        if core_resname != "HIS":
            return None

        has_hd1 = "HD1" in atom_names
        has_he2 = bool({"HE2", "HNE2"} & atom_names)

        if has_hd1 and has_he2:
            return f"{prefix}HSP"
        if has_hd1:
            return f"{prefix}HSD"
        if has_he2:
            return f"{prefix}HSE"

        return None

    @staticmethod
    def _looks_like_n_terminal(residue) -> bool:
        """
        Heuristically detect an N-terminal amino-acid residue from its atom names.

        This is used as a fallback when fragment/chain metadata are not sufficient
        to identify all chains in the system (for example, protein dimers whose
        monomers share the same chain identifier).

        :param residue
            Protein residue to inspect.
        :return (bool)
            True if the residue looks like an N-terminus, False otherwise.
        """
        atom_names = {str(name) for name in residue.atoms.names}

        # Common terminal-ammonium naming conventions
        if {"H1", "H2", "H3"} <= atom_names:
            return True
        if {"HT1", "HT2", "HT3"} <= atom_names:
            return True

        # More permissive fallback used by some force fields/topologies
        if "N" in atom_names and ("H2" in atom_names or "H3" in atom_names):
            return True
        if "N" in atom_names and ("HT2" in atom_names or "HT3" in atom_names):
            return True

        return False

    @staticmethod
    def _looks_like_c_terminal(residue) -> bool:
        """
        Heuristically detect a C-terminal amino-acid residue from its atom names.

        :param residue
            Protein residue to inspect.
        :return (bool)
            True if the residue looks like a C-terminus, False otherwise.
        """
        atom_names = {str(name) for name in residue.atoms.names}

        # Common terminal-carboxylate naming conventions
        if {"OXT", "OT1", "OT2", "OC1", "OC2"} & atom_names:
            return True

        # Additional generic fallback seen in some topologies
        if "O" in atom_names and ({"O1", "O2"} & atom_names):
            return True

        return False

    def _terminal_resname_map(self, mda_universe):
        """
        Identifies N- and C-terminal protein residues and returns a renaming map.

        For each protein chain found in the selection, the first residue is
        considered the N-terminus and the last residue the C-terminus.
        These residues are assigned terminal-prefixed
        residue names (e.g. 'ASN' -> 'NASN' or 'CASN').

        The returned mapping is keyed by MDAnalysis `resindex` (unique per
        residue in the Universe) and can be used to update per-atom `resnames`
        arrays (e.g. `AtomGroup.resnames`) based on `AtomGroup.resindices`.

        :param mda_universe
            The MDAnalysis universe object.
        :return dict[int, str]
            Dictionary mapping residue `resindex` to the renamed residue name
            with terminal prefix. Only protein terminal residues are included.

        Note: Chain detection is first attempted using connectivity information
        (`fragments`) and then chain IDs / segment IDs. As a conservative
        fallback, residues with atom-name patterns characteristic of N- or
        C-termini are also recognized. This helps when multiple protein chains
        are present but the topology does not distinguish them cleanly.
        """
        # Restrict to protein residues
        prot = mda_universe.select_atoms("protein")
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
            print_exception_if_debug()
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

        # Additional fallback based on terminal-specific atom names.
        # This is especially useful when several protein chains are present
        # but fragment/chain metadata collapse them into a single group.
        for res in prot.residues:
            ridx = int(res.resindex)

            if ridx not in term_map and self._looks_like_n_terminal(res):
                term_map[ridx] = self._prefixed_resname(res.resname, "N")

            if ridx not in term_map and self._looks_like_c_terminal(res):
                term_map[ridx] = self._prefixed_resname(res.resname, "C")

        return term_map

    def _protonation_resname_map(self, env_atoms):
        """
        Identify simple protonation-state variants from residue atom names and
        return a renaming map keyed by MDAnalysis ``resindex``.

        The purpose of this routine is analogous to ``_terminal_resname_map``:
        it refines residue names before parameter lookup so that entries such as
        GLU/GLH, ASP/ASH, can be distinguished downstream without
        modifying the writer logic.

        The current heuristics are:

        - GLU + HE1/HE2 -> GLH
        - ASP + HD1/HD2 -> ASH
        - CYS without HG/HG1 -> CYX
        - HIS + HD1/HE2 -> HSD/HSE/HSP

        :param env_atoms (MDAnalysis.core.groups.AtomGroup)
            AtomGroup corresponding to the environment selection used in
            trajectory parsing.

        :return dict[int, str]
            Dictionary mapping residue ``resindex`` to refined residue names.
            Only residues that need renaming are included.
        """
        prot_map = {}

        for res in env_atoms.residues:
            resname = str(res.resname)
            atom_names = {str(name) for name in res.atoms.names}
            prefix = ""
            core_resname = resname

            if len(resname) == 4 and resname[0] in ("N", "C"):
                prefix = resname[0]
                core_resname = resname[1:]

            histidine_resname = self._histidine_resname_from_atoms(
                resname, atom_names)
            if histidine_resname is not None:
                prot_map[res.resindex] = histidine_resname
            elif core_resname == "GLU" and ({"HE1", "HE2"} & atom_names):
                prot_map[res.resindex] = f"{prefix}GLH"
            elif core_resname == "ASP" and ({"HD1", "HD2"} & atom_names):
                prot_map[res.resindex] = f"{prefix}ASH"
            elif core_resname == "CYS" and not ({"HG", "HG1"} & atom_names):
                # Distinguish thiol-less/disulfide cysteine from protonated CYS.
                prot_map[res.resindex] = f"{prefix}CYX"

        return prot_map

    @staticmethod
    def _apply_resname_maps(resnames, resindices, prot_map, term_map,
                            prefixed_resname):
        """
        Applies protonation-state and terminal residue renaming in place.
        """
        if len(resnames) == 0:
            return

        # Apply protonation-state residue renaming (GLH/ASH) if applicable
        if prot_map:
            for ridx, newname in prot_map.items():
                resnames[resindices == ridx] = newname

        # Apply terminal residue renaming (NASN/CASN, etc.) if applicable
        if term_map:
            for ridx, newname in term_map.items():
                term_newname = newname
                if ridx in prot_map:
                    prot_newname = str(prot_map[ridx])
                    term_newname = prefixed_resname(
                        prot_newname, str(newname)[0])
                resnames[resindices == ridx] = term_newname

    def _load_universe(self,
                       trajectory_file=None,
                       topology_file=None,
                       pdb_file=None):
        """
        Creates an MDAnalysis universe for a PDB file or a trajectory.
        """
        if pdb_file is not None and trajectory_file is not None:
            raise ValueError(
                "Provide either pdb_file or trajectory_file, not both")

        if pdb_file is not None:
            pdb_file = str(pdb_file)
            if not pdb_file.lower().endswith(".pdb"):
                raise ValueError("pdb_file must be a .pdb file")
            trajectory_file = pdb_file

        if trajectory_file is None:
            raise ValueError(
                "Either pdb_file or trajectory_file must be provided")

        trajectory_file = str(trajectory_file)

        assert_msg_critical(
            'MDAnalysis' in sys.modules,
            'MDAnalysis is required for EnsembleParser.')

        if trajectory_file.lower().endswith(".pdb"):
            if self._pdb_has_conect_records(trajectory_file):
                # Honor the connectivity records of the file.  Bond guessing
                # would otherwise turn short intermolecular contacts (e.g.
                # hydrogen bonds between a protein and crystallographic
                # waters) into covalent bonds and corrupt the chain and
                # fragment detection used for the terminal residue names.
                mda_universe = mda.Universe(trajectory_file)
            else:
                raw_mda_version = version("MDAnalysis")
                match = re.match(r"^(\d+)\.(\d+)", raw_mda_version)
                if not match:
                    raise RuntimeError(
                        "Cannot parse required major.minor version from "
                        f"MDAnalysis: {raw_mda_version!r}")
                mda_version = int(match.group(1)), int(match.group(2))
                if mda_version >= (2, 8):
                    # https://docs.mdanalysis.org/2.8.0/documentation_pages/core/universe.html
                    # The `guess_bonds` keyword is deprecated and will be removed in MDAnalysis 3.0.
                    # Please pass ("bonds", "angles", "dihedrals") into to_guess or force_guess
                    # instead to guess bonds, angles, and dihedrals respectively.
                    mda_universe = mda.Universe(
                        trajectory_file,
                        to_guess=("bonds", "angles", "dihedrals"),
                    )
                else:
                    mda_universe = mda.Universe(
                        trajectory_file,
                        guess_bonds=True,
                    )
        else:
            if topology_file is None:
                raise ValueError(
                    "topology_file is required unless trajectory_file is a .pdb")
            # Refresh XDR offsets to avid stale .xtc_offsets cache warinings
            # when trajectory metadata changes between runs/environments:
            mda_universe = mda.Universe(
                topology_file, trajectory_file, refresh_offsets=True)

        return mda_universe, trajectory_file

    @staticmethod
    def _seek_trajectory_frame(mda_universe, iframe):
        """
        Moves the trajectory cursor to a frame and returns its timestep.
        """
        return mda_universe.trajectory[iframe]

    def _frame_indices(self,
                       mda_universe,
                       trajectory_file,
                       num_snapshots=None,
                       start=None,
                       end=None,
                       last_snapshot_only=False,
                       sampling_times=None):
        """
        Determines the trajectory frames to be processed.
        """
        total_frames = len(mda_universe.trajectory)
        is_pdb = trajectory_file.lower().endswith(".pdb")

        if sampling_times is not None:
            sampling_times = np.asarray(sampling_times, dtype=float).reshape(-1)
            frame_indices = []
            for iframe in range(total_frames):
                t = float(self._seek_trajectory_frame(
                    mda_universe, iframe).time)
                if np.min(np.abs(sampling_times - t), initial=np.inf) <= 1.0e-6:
                    frame_indices.append(iframe)
            if not frame_indices:
                raise ValueError(
                    "No trajectory frames found at the requested sampling "
                    "times.")
            return np.array(frame_indices, dtype=int)

        if last_snapshot_only:
            return np.array([total_frames - 1], dtype=int)

        use_time_window = (start is not None or end is not None)
        if use_time_window and is_pdb:
            raise ValueError(
                "Start/end time window selection is only supported for "
                "time-resolved trajectories (e.g. .xtc), not .pdb input.")

        if use_time_window:
            traj_start = float(
                self._seek_trajectory_frame(mda_universe, 0).time)
            traj_end = float(
                self._seek_trajectory_frame(
                    mda_universe, total_frames - 1).time)

            if start is None:
                start = traj_start
            if end is None:
                end = traj_end

            if float(end) < float(start):
                raise ValueError(
                    "End time must be greater than or equal to start time")

            window_frames = []
            for iframe in range(total_frames):
                t = float(self._seek_trajectory_frame(
                    mda_universe, iframe).time)
                if float(start) <= t <= float(end):
                    window_frames.append(iframe)

            if len(window_frames) == 0:
                raise ValueError(
                    "No trajectory frames found in the requested time window "
                    f"[{start}, {end}] ps.")

            available_frames = np.array(window_frames, dtype=int)
        else:
            available_frames = np.arange(total_frames, dtype=int)

        if num_snapshots is None:
            num_snapshots = len(available_frames)

        if num_snapshots <= 0:
            raise ValueError("num_snapshots must be a positive integer")

        if num_snapshots > len(available_frames):
            raise ValueError(
                f"Requested number of snapshots ({num_snapshots}) exceeds "
                f"available frames ({len(available_frames)}).")

        if num_snapshots == len(available_frames):
            return available_frames
        if num_snapshots == 1:
            return np.array([available_frames[0]], dtype=int)

        sample_idx = np.linspace(
            0, len(available_frames) - 1, num_snapshots, dtype=int)
        return available_frames[sample_idx]

    @staticmethod
    def _match_molecule_to_atoms(mda_universe, molecule, tolerance=0.1):
        """
        Matches the atoms of a molecule to atoms of the universe by element and
        Cartesian coordinates at the current frame.

        :return:
            An MDAnalysis AtomGroup in molecule atom order, or None if no
            complete match is found.
        """
        if molecule is None or molecule.number_of_atoms() == 0:
            return None

        labels = [str(x).upper() for x in molecule.get_labels()]
        target_elements = np.array(
            [chemical_element_identifier(x) for x in labels], dtype=int)
        target_coords = np.asarray(
            molecule.get_coordinates_in_angstrom(), dtype=float)

        guesser = mda_DefaultGuesser(mda_universe)
        universe_elements = np.array([
            chemical_element_identifier(
                str(guesser.guess_atom_element(name)).upper())
            for name in mda_universe.atoms.names
        ], dtype=int)
        universe_coords = np.asarray(mda_universe.atoms.positions, dtype=float)

        used = np.zeros(len(universe_coords), dtype=bool)
        matched = []
        for elem, xyz in zip(target_elements, target_coords):
            candidates = np.where((universe_elements == elem) & (~used))[0]
            if candidates.size == 0:
                return None
            distances = np.linalg.norm(
                universe_coords[candidates] - xyz, axis=1)
            nearest = int(candidates[np.argmin(distances)])
            if distances.min() > tolerance:
                return None
            used[nearest] = True
            matched.append(nearest)

        if not matched:
            return None
        return mda_universe.atoms[matched]

    def _select_cutoff_residues(self,
                                mda_universe,
                                env_atoms,
                                qm_atoms,
                                qm_coords,
                                cutoff,
                                env_region_sel):
        """
        Selects environment residues within a cutoff of the QM region.
        """
        if qm_atoms is not None and len(qm_atoms) > 0:
            touch = mda_universe.select_atoms(
                f"({env_region_sel} and around {float(cutoff)} group qm)",
                qm=qm_atoms,
            ).difference(qm_atoms)
        else:
            qm_coords = np.asarray(qm_coords, dtype=float)
            env_coords = np.asarray(env_atoms.positions, dtype=float)
            if env_coords.size == 0:
                return np.empty((0,), dtype=int)
            distances = np.linalg.norm(
                env_coords[:, None, :] - qm_coords[None, :, :], axis=2)
            touch = env_atoms[np.any(distances <= float(cutoff), axis=1)]

        return np.unique(np.asarray(touch.resindices, dtype=int))

    def _collect_structures(self,
                            trajectory_file: str | None = None,
                            num_snapshots: int | None = None,
                            qm_region: str = "",
                            qm_charge: int = 0,
                            qm_multiplicity: int = 1,
                            env_region: str | None = None,
                            topology_file: str | None = None,
                            pe_cutoff: float | None = None,
                            npe_cutoff: float | None = None,
                            start: float | None = None,
                            end: float | None = None,
                            last_snapshot_only: bool = False,
                            pdb_file: str | None = None,
                            sampling_times=None,
                            qm_molecule=None):
        """
        Parses structures and returns a list of snapshot dictionaries.
        """
        mda_universe, trajectory_file = self._load_universe(
            trajectory_file=trajectory_file,
            topology_file=topology_file,
            pdb_file=pdb_file)

        if pe_cutoff is not None and npe_cutoff is not None:
            if float(npe_cutoff) < float(pe_cutoff):
                raise ValueError("npe_cutoff must be >= pe_cutoff")

        if qm_molecule is not None:
            if qm_charge == 0:
                qm_charge = int(qm_molecule.get_charge())
            if qm_multiplicity == 1:
                qm_multiplicity = int(qm_molecule.get_multiplicity())

        qm_charge = int(qm_charge)
        qm_multiplicity = int(qm_multiplicity)

        if qm_multiplicity <= 0:
            raise ValueError("qm_multiplicity must be a positive integer")

        self.ostream.print_blank()

        frame_indices = self._frame_indices(
            mda_universe,
            trajectory_file,
            num_snapshots=num_snapshots,
            start=start,
            end=end,
            last_snapshot_only=last_snapshot_only,
            sampling_times=sampling_times)

        qm_region = "" if qm_region is None else str(qm_region).strip()
        external_qm = None
        qm_atoms = None

        if qm_region:
            qm_atoms = mda_universe.select_atoms(qm_region)
            if len(qm_atoms) == 0:
                raise ValueError(
                    f"QM region '{qm_region}' selection is empty")
        elif qm_molecule is not None and qm_molecule.number_of_atoms() > 0:
            self._seek_trajectory_frame(mda_universe, int(frame_indices[0]))
            qm_atoms = self._match_molecule_to_atoms(mda_universe, qm_molecule)
            if qm_atoms is None:
                # The PDB/trajectory describes only the environment. Use the
                # supplied molecule as the fixed QM subsystem.
                external_qm = qm_molecule
        else:
            qm_atoms = mda_universe.select_atoms("resname MOL")
            if len(qm_atoms) == 0:
                raise ValueError(
                    "QM region selection is empty. Provide qm_region or a QM "
                    "molecule when the structure contains no residue named "
                    "'MOL'.")

        if env_region is None or not str(env_region).strip():
            if qm_region:
                env_region_sel = f"not ({qm_region})"
            else:
                env_region_sel = "all"
        else:
            env_region_sel = str(env_region)

        self._seek_trajectory_frame(mda_universe, int(frame_indices[0]))
        env_atoms = mda_universe.select_atoms(env_region_sel)
        if qm_atoms is not None:
            env_atoms = env_atoms.difference(qm_atoms)

        # Identify terminal protein residues (if any) and assign N*/C* residue names
        # so that terminal variants can be treated as separate residue types downstream.
        prot_map = self._protonation_resname_map(env_atoms)
        term_map = self._terminal_resname_map(mda_universe)

        # The unwrap, center_in_box and wrap operations require valid box
        # dimensions (ts.dimensions). For single PDBs without box info, skip
        # them.
        has_box = False
        try:
            dims = getattr(
                self._seek_trajectory_frame(
                    mda_universe, int(frame_indices[0])), "dimensions", None)
            if dims is not None:
                dims = np.asarray(dims, dtype=float)
                if dims.size > 3 and np.all(dims[:3] > 0):
                    has_box = True
        except Exception:
            print_exception_if_debug()
            has_box = False

        atom_guesser = mda_DefaultGuesser(mda_universe)

        snapshots = []
        for iframe in frame_indices:
            ts = self._seek_trajectory_frame(mda_universe, iframe)

            if has_box and external_qm is None:

                # apply the unwrap, center_in_box and wrap operations directly
                # on the coordinates, instead of through the MDAnalysis
                # transformations, which adjust the thread limits of the math
                # libraries of the process as a side effect

                for fragment in qm_atoms.fragments:
                    make_whole(fragment)

                box_center = np.sum(ts.triclinic_dimensions, axis=0) / 2

                shift = box_center - qm_atoms.center_of_geometry(wrap=True)

                ts.positions += shift

                env_atoms.wrap()

            if external_qm is not None:
                qm_coords = np.asarray(
                    external_qm.get_coordinates_in_angstrom(), dtype=float).copy()
                qm_atom_names = np.asarray(
                    [str(x) for x in external_qm.get_labels()], dtype=object)
                qm_elements = qm_atom_names.copy()
                qm_selection = None
            else:
                qm_coords = np.asarray(qm_atoms.positions, dtype=float).copy()
                qm_atom_names = np.asarray(qm_atoms.names, dtype=object).copy()
                qm_elements = np.asarray(
                    [atom_guesser.guess_atom_element(n) for n in qm_atoms.names],
                    dtype=object)
                qm_selection = qm_atoms

            pe_coords = None
            pe_elements = None
            pe_resids = None
            pe_resindices = None
            pe_resnames = None
            number_residues_pe = 0
            pe_atom_names = None

            npe_coords = None
            npe_elements = None
            npe_resids = None
            npe_resindices = None
            npe_resnames = None
            number_residues_npe = 0
            npe_atom_names = None

            pe_residx_unique = np.empty((0,), dtype=int)

            # PE selection
            if pe_cutoff is not None:
                pe_residx_unique = self._select_cutoff_residues(
                    mda_universe,
                    env_atoms,
                    qm_selection,
                    qm_coords,
                    pe_cutoff,
                    env_region_sel,
                )

                if pe_residx_unique.size > 0:
                    pe_region = env_atoms[
                        np.isin(env_atoms.resindices, pe_residx_unique)]
                else:
                    pe_region = env_atoms[:0]

                pe_coords = np.asarray(
                    pe_region.positions, dtype=float).copy()
                pe_atom_names = np.asarray(
                    pe_region.names, dtype=object).copy()
                pe_elements = np.asarray([
                    atom_guesser.guess_atom_element(n)
                    for n in pe_region.names
                ], dtype=object)
                pe_resids = np.asarray(pe_region.resids, dtype=int).copy()
                pe_resindices = np.asarray(
                    pe_region.resindices, dtype=int).copy()
                pe_resnames = np.asarray(
                    pe_region.resnames, dtype=object).copy()

                self._apply_resname_maps(pe_resnames, pe_resindices, prot_map,
                                         term_map, self._prefixed_resname)

                number_residues_pe = int(pe_region.residues.n_residues)

            # NPE selection
            if npe_cutoff is not None:
                outer_residx_unique = self._select_cutoff_residues(
                    mda_universe,
                    env_atoms,
                    qm_selection,
                    qm_coords,
                    npe_cutoff,
                    env_region_sel,
                )
                npe_residx_unique = np.setdiff1d(
                    outer_residx_unique,
                    pe_residx_unique,
                    assume_unique=False,
                )

                if npe_residx_unique.size > 0:
                    npe_region = env_atoms[
                        np.isin(env_atoms.resindices, npe_residx_unique)]
                else:
                    npe_region = env_atoms[:0]

                npe_coords = np.asarray(
                    npe_region.positions, dtype=float).copy()
                npe_atom_names = np.asarray(
                    npe_region.names, dtype=object).copy()
                npe_elements = np.asarray([
                    atom_guesser.guess_atom_element(n)
                    for n in npe_region.names
                ], dtype=object)
                npe_resids = np.asarray(npe_region.resids, dtype=int).copy()
                npe_resindices = np.asarray(
                    npe_region.resindices, dtype=int).copy()
                npe_resnames = np.asarray(
                    npe_region.resnames, dtype=object).copy()

                self._apply_resname_maps(npe_resnames, npe_resindices,
                                         prot_map, term_map,
                                         self._prefixed_resname)

                number_residues_npe = int(npe_region.residues.n_residues)

            # Normalize missing PE/NPE selections to fresh empty arrays
            # per snapshot
            if pe_coords is None:
                pe_coords = np.empty((0, 3), dtype=float)
                pe_atom_names = np.empty((0,), dtype=object)
                pe_elements = np.empty((0,), dtype=object)
                pe_resids = np.empty((0,), dtype=int)
                pe_resindices = np.empty((0,), dtype=int)
                pe_resnames = np.empty((0,), dtype=object)

            if npe_coords is None:
                npe_coords = np.empty((0, 3), dtype=float)
                npe_atom_names = np.empty((0,), dtype=object)
                npe_elements = np.empty((0,), dtype=object)
                npe_resids = np.empty((0,), dtype=int)
                npe_resindices = np.empty((0,), dtype=int)
                npe_resnames = np.empty((0,), dtype=object)

            # If neither cutoff is set, interpret as all-qm

            snapshot = {
                "frame": int(mda_universe.trajectory.frame),
                "trajectory_time": float(ts.time),

                "qm_coords": qm_coords,
                "qm_atom_names": qm_atom_names,
                "qm_elements": qm_elements,
                "qm_charge": qm_charge,
                "qm_multiplicity": qm_multiplicity,

                "pe_coords": pe_coords,
                "pe_atom_names": pe_atom_names,
                "pe_elements": pe_elements,
                "pe_resids": pe_resids,
                "pe_resindices": pe_resindices,
                "pe_resnames": pe_resnames,
                "number_residues_pe": number_residues_pe,

                "npe_coords": npe_coords,
                "npe_atom_names": npe_atom_names,
                "npe_elements": npe_elements,
                "npe_resids": npe_resids,
                "npe_resindices": npe_resindices,
                "npe_resnames": npe_resnames,
                "number_residues_npe": number_residues_npe,
            }
            snapshots.append(snapshot)

        return snapshots

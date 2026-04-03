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
from pathlib import Path
import numpy as np

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .molecule import Molecule


def _resolve_element(token):
    """Return an element symbol from a token that is either a symbol or an
    atomic number (integer string such as '6' for carbon)."""
    try:
        z = int(token)
        return Molecule([z], [[0.0, 0.0, 0.0]], units='bohr').get_label(0)
    except (ValueError, Exception):
        pass  # not an integer — treat as symbol
    return token.capitalize()


class QMTrajectoryParser:
    """
    Parser for QM-only molecular dynamics trajectories.

    Converts trajectory files into sequentially ordered snapshot dictionaries
    suitable for :class:`QMTrajectoryDriver`. Unlike :class:`EnsembleParser`,
    this class handles no environment (PE/NPE) region — every selected atom is
    treated as part of the QM system.

    Supported formats:

    - **Multi-frame XYZ** (``.xyz``): parsed natively, no external dependency
      required.  This is the most common output format from AIMD codes such as
      CP2K, ORCA, and NWChem.
    - **GROMACS trajectory** (``.xtc`` + ``.tpr`` topology): parsed via
      MDAnalysis.
    - **PDB multi-model** (``.pdb``): parsed via MDAnalysis.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialises the QMTrajectoryParser.
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
        self.size = self.comm.Get_size()
        self.ostream = ostream

    @staticmethod
    def _parse_xyz_frames(trajectory_file):
        """
        Parse all frames from a multi-frame XYZ file.

        The expected format is the standard sequential multi-frame XYZ layout::

            N
            comment line (ignored)
            Elem  x  y  z
            ...
            N
            comment line
            Elem  x  y  z
            ...

        Blank lines between frames are tolerated.  Coordinates are assumed to
        be in Ångström, which is the standard XYZ convention.

        :param trajectory_file:
            Path to the XYZ file.
        :return:
            List of ``(elements, coords_angstrom)`` tuples, one per frame.
            ``elements`` is a ``numpy.ndarray`` of shape ``(N,)`` containing
            element symbol strings; ``coords_angstrom`` is a
            ``numpy.ndarray`` of shape ``(N, 3)``.
        """
        frames = []
        path = Path(trajectory_file)

        with path.open("r") as fh:
            lines = fh.readlines()

        i = 0
        while i < len(lines):
            # Skip blank lines between frames.
            line = lines[i].strip()
            if not line:
                i += 1
                continue

            try:
                n_atoms = int(line)
            except ValueError:
                raise ValueError(
                    f"Expected atom count at line {i + 1}, got: {line!r}"
                )

            i += 1  # skip the comment/title line
            if i >= len(lines):
                break
            i += 1  # advance to the first atom line

            if i + n_atoms > len(lines):
                raise ValueError(
                    f"XYZ file truncated at frame {len(frames) + 1}: expected "
                    f"{n_atoms} atom lines but only {len(lines) - i} remain."
                )

            elements = []
            coords = []
            for j in range(n_atoms):
                parts = lines[i + j].split()
                if len(parts) < 4:
                    raise ValueError(
                        f"Malformed XYZ atom line {i + j + 1}: "
                        f"{lines[i + j].rstrip()!r}"
                    )
                elements.append(_resolve_element(parts[0]))
                coords.append(
                    [float(parts[1]), float(parts[2]), float(parts[3])]
                )

            frames.append((
                np.asarray(elements, dtype=object),
                np.asarray(coords, dtype=float),
            ))
            i += n_atoms

        return frames

    def structures(
        self,
        trajectory_file: str,
        topology_file: str | None = None,
        num_frames: int | None = None,
        qm_region: str | None = None,
        qm_charge: int = 0,
        qm_multiplicity: int = 1,
        start_frame: int | None = None,
        end_frame: int | None = None,
        start: float | None = None,
        end: float | None = None,
    ) -> list[dict]:
        """
        Parse a trajectory and return sequentially ordered snapshot
        dictionaries for :class:`QMTrajectoryDriver`.

        Dispatches to a built-in XYZ reader for ``.xyz`` files, and to
        MDAnalysis for ``.xtc``/``.pdb`` files.

        :param trajectory_file:
            Path to the trajectory file (``*.xyz``, ``*.xtc``, or ``*.pdb``).
        :param topology_file:
            Path to the topology file.  Required for ``.xtc``; ignored for
            ``.xyz`` and ``.pdb``.
        :param num_frames:
            Maximum number of frames to extract.  Frames are taken
            **sequentially from the start** of the selected window — no
            sub-sampling is applied.  If ``None``, all frames in the window
            are returned.
        :param qm_region:
            MDAnalysis atom-selection string defining the QM atoms.  Only used
            for MDAnalysis-based trajectories (``.xtc``/``.pdb``).  If
            ``None``, all atoms are treated as the QM region.  Ignored for
            ``.xyz`` files (all atoms are always QM).
        :param qm_charge:
            Total charge of the QM system.  Default is 0.
        :param qm_multiplicity:
            Spin multiplicity of the QM system.  Default is 1 (singlet).
        :param start_frame:
            First frame to include (0-based index, inclusive).  Only used for
            ``.xyz`` files.  Default is 0.
        :param end_frame:
            One-past-the-last frame index (exclusive).  Only used for ``.xyz``
            files.  Default is the total number of frames in the file.
        :param start:
            Start time in ps (inclusive).  Only used for MDAnalysis-based
            trajectories.  Ignored for ``.xyz`` files.
        :param end:
            End time in ps (inclusive).  Only used for MDAnalysis-based
            trajectories.  Ignored for ``.xyz`` files.
        :return:
            List of snapshot dictionaries, each containing:

            - ``frame`` (int): original frame number from the trajectory file.
            - ``frame_index`` (int): 0-based sequential index in the returned list.
            - ``trajectory_name`` (str): base name of the trajectory file.
            - ``total_frames`` (int): total number of snapshots returned.
            - ``qm_coords`` (numpy.ndarray): shape ``(N, 3)``, in Å.
            - ``qm_elements`` (numpy.ndarray): shape ``(N,)``, element symbols.
            - ``qm_charge`` (int): total charge of the QM system.
            - ``qm_multiplicity`` (int): spin multiplicity.
            - ``qm_atom_indices`` (numpy.ndarray): shape ``(N,)``; sequential
              0-based indices for XYZ files, or MDAnalysis universe indices
              for ``.xtc``/``.pdb`` files.
        """
        qm_charge = int(qm_charge)
        if int(qm_multiplicity) <= 0:
            raise ValueError("qm_multiplicity must be a positive integer.")
        qm_multiplicity = int(qm_multiplicity)

        trajectory_file = str(trajectory_file)
        trajectory_name = Path(trajectory_file).name

        if trajectory_file.lower().endswith(".xyz"):
            return self._structures_xyz(
                trajectory_file, trajectory_name, num_frames,
                qm_charge, qm_multiplicity, start_frame, end_frame,
            )
        else:
            return self._structures_mda(
                trajectory_file, topology_file, trajectory_name, num_frames,
                qm_region, qm_charge, qm_multiplicity, start_frame, end_frame,
            )

    def _structures_xyz(
        self,
        trajectory_file,
        trajectory_name,
        num_frames,
        qm_charge,
        qm_multiplicity,
        start_frame,
        end_frame,
    ):
        """Parse a multi-frame XYZ trajectory and return snapshots."""
        raw_frames = self._parse_xyz_frames(trajectory_file)
        total_available = len(raw_frames)

        if total_available == 0:
            raise ValueError(f"No frames found in {trajectory_file!r}.")

        s = 0 if start_frame is None else max(0, int(start_frame))
        e = (
            total_available
            if end_frame is None
            else min(total_available, int(end_frame))
        )

        if s >= e:
            raise ValueError(
                f"start_frame ({s}) must be less than end_frame ({e})."
            )

        window = raw_frames[s:e]

        if num_frames is not None:
            num_frames = int(num_frames)
            if num_frames <= 0:
                raise ValueError("num_frames must be a positive integer.")
            if num_frames > len(window):
                raise ValueError(
                    f"Requested num_frames ({num_frames}) exceeds the "
                    f"{len(window)} frames available in the window "
                    f"[{s}, {e})."
                )
            window = window[:num_frames]

        total_frames = len(window)
        snapshots = []

        for seq_idx, (elements, coords) in enumerate(window):
            orig_frame = s + seq_idx
            snapshots.append({
                "frame": orig_frame,
                "frame_index": seq_idx,
                "trajectory_name": trajectory_name,
                "total_frames": total_frames,
                "qm_coords": coords.copy(),
                "qm_elements": elements.copy(),
                "qm_charge": qm_charge,
                "qm_multiplicity": qm_multiplicity,
                "qm_atom_indices": np.arange(elements.size, dtype=int),
            })

        if self.rank == mpi_master():
            self.ostream.print_info(
                f"QMTrajectoryParser: {total_frames} frames parsed from "
                f"{trajectory_name} (XYZ)."
            )
            self.ostream.print_blank()

        return snapshots

    def _structures_mda(
        self,
        trajectory_file,
        topology_file,
        trajectory_name,
        num_frames,
        qm_region,
        qm_charge,
        qm_multiplicity,
        start_frame,
        end_frame,
    ):
        """Parse a trajectory via MDAnalysis and return snapshots."""
        try:
            import MDAnalysis as mda
            from MDAnalysis.topology.guessers import guess_atom_element
        except ImportError:
            raise ImportError(
                "MDAnalysis is required for .xtc/.pdb trajectory parsing. "
                "Install it with: conda install -c conda-forge mdanalysis"
            )

        if trajectory_file.lower().endswith(".pdb"):
            universe = mda.Universe(trajectory_file, guess_bonds=True)
        else:
            if topology_file is None:
                raise ValueError(
                    "topology_file is required for non-XYZ trajectories "
                    "(e.g. .xtc)."
                )
            universe = mda.Universe(topology_file, trajectory_file)

        total_traj_frames = len(universe.trajectory)

        if qm_region is not None and str(qm_region).strip():
            qm_atoms = universe.select_atoms(qm_region)
            if len(qm_atoms) == 0:
                raise ValueError(
                    f"QM region selection '{qm_region}' matched no atoms."
                )
        else:
            qm_atoms = universe.atoms

        s = 0 if start_frame is None else max(0, int(start_frame))
        e = (
            total_traj_frames
            if end_frame is None
            else min(total_traj_frames, int(end_frame))
        )
        if s >= e:
            raise ValueError(
                f"start_frame ({s}) must be less than end_frame ({e})."
            )
        frame_indices = np.arange(s, e, dtype=int)

        if num_frames is not None:
            num_frames = int(num_frames)
            if num_frames <= 0:
                raise ValueError("num_frames must be a positive integer.")
            if num_frames > len(frame_indices):
                raise ValueError(
                    f"Requested num_frames ({num_frames}) exceeds the "
                    f"{len(frame_indices)} frames available in the selected "
                    f"window."
                )
            frame_indices = frame_indices[:num_frames]

        total_frames = len(frame_indices)
        snapshots = []

        for seq_idx, iframe in enumerate(frame_indices):
            universe.trajectory[iframe]

            qm_coords = np.asarray(qm_atoms.positions, dtype=float).copy()
            _elems = []
            for atom in qm_atoms:
                el = getattr(atom, 'element', '').strip()
                _elems.append(el.capitalize() if el else guess_atom_element(atom.name))
            qm_elements = np.asarray(_elems, dtype=object)
            qm_atom_indices = np.asarray(qm_atoms.indices, dtype=int).copy()

            snapshots.append({
                "frame": int(universe.trajectory.frame),
                "frame_index": seq_idx,
                "trajectory_name": trajectory_name,
                "total_frames": total_frames,
                "qm_coords": qm_coords,
                "qm_elements": qm_elements,
                "qm_charge": qm_charge,
                "qm_multiplicity": qm_multiplicity,
                "qm_atom_indices": qm_atom_indices,
            })

        if self.rank == mpi_master():
            self.ostream.print_info(
                f"QMTrajectoryParser: {total_frames} frames parsed from "
                f"{trajectory_name}."
            )
            self.ostream.print_blank()

        return snapshots

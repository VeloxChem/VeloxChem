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

from .trajectoryparserbase import TrajectoryParserBase


class EnsembleParser(TrajectoryParserBase):
    """
    This class automates the parsing of molecular dynamics trajectories and
    extracts information from the qm and environment regions for each snapshot.

    It shares the per-frame parsing implementation with :class:`EmbeddingParser`
    through :class:`TrajectoryParserBase`, but is intentionally not a subclass of
    :class:`EmbeddingParser`: an ensemble is not a kind of single embedding,
    and the public single-structure API should not be used once per frame
    (that would reload the universe and repeat the topology setup).

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the EnsembleParser.
        """
        super().__init__(comm=comm, ostream=ostream)

    def structures(self,
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
        Parse a set of structures and extract QM and MM region data.

        :param trajectory_file:
            Path to the trajectory file:
                - .xtc with a corresponding topology (.tpr) via topology_file
                - .pdb (several configurations, or a single configuration;
                  CONECT records are honored when present, otherwise bonds are
                  guessed)
        :param pdb_file:
            Path to a PDB file. This is a convenience alternative to
            trajectory_file and cannot be used together with trajectory_file.
        :param topology_file:
            Path to the topology file (e.g., .tpr).
        :param num_snapshots:
            Number of snapshots to extract. If None (default), all frames are
            processed. If an integer, frames are sampled evenly.
        :param qm_region:
            Selection string defining the QM region (MDAnalysis selection syntax).
        :param qm_charge:
            Total charge of the QM region. Default is 0.
        :param qm_multiplicity:
            Spin multiplicity of the QM region. Default is 1 (singlet).
        :param env_region:
            Selection string defining the environment region (MDAnalysis selection syntax).
            If not provided, the environment is taken as all atoms not contained in
            the QM region.
        :param pe_cutoff:
            Cutoff for polarizable embedding (PE) region selection (Angstrom).
        :param npe_cutoff:
            Cutoff for non-polarizable embedding (NPE) region selection (Angstrom).
        :param start:
            Start time of the trajectory window in ps. If None, the first
            trajectory time is used. Intended for time-resolved trajectories
            such as .xtc.
        :param end:
            End time of the trajectory window in ps. If None, the last
            trajectory time is used. Intended for time-resolved trajectories
            such as .xtc.
        :param last_snapshot_only:
            If True, process only the last snapshot in the trajectory. When
            enabled, start, end, and num_snapshots are ignored.
        :param sampling_times:
            Optional sequence of trajectory times (ps) to extract exactly.
            When provided, it takes precedence over last_snapshot_only,
            start/end and num_snapshots.
        :param qm_molecule:
            Optional QM molecule. If qm_region is not given, the molecule is
            matched to the structure; when the structure contains only the
            environment, the molecule is used as the fixed QM subsystem.

        :return:
            A list of snapshot dictionaries, each containing:
            - frame (int):
                Frame number.
            - trajectory_time (float):
                Time of the frame in ps.
            - qm_coords (numpy.ndarray):
                QM region Cartesian coordinates, shape (N_qm, 3), in Angstrom.
            - qm_atom_names (numpy.ndarray):
                Atom names for each QM atom, shape (N_qm,).
            - qm_elements (numpy.ndarray):
                Element symbols for each QM atom, shape (N_qm,).
            - qm_charge (int):
                Total charge of the QM region.
            - qm_multiplicity (int):
                Spin multiplicity of the QM region.
            - pe_coords (numpy.ndarray):
                PE region Cartesian coordinates, shape (N_pe, 3), in Angstrom.
            - pe_atom_names (numpy.ndarray):
                Atom names for each PE atom, shape (N_pe,).
            - pe_elements (numpy.ndarray):
                Element symbols for each PE atom, shape (N_pe,).
            - pe_resids (numpy.ndarray):
                Residue id for each PE atom, shape (N_pe,).
            - pe_resindices (numpy.ndarray):
                MDAnalysis residue index for each PE atom, shape (N_pe,).
                This is unique per residue in the Universe and is used
                internally to keep PE/NPE residue partitioning robust even
                when `resid` values are reused across chains.
            - pe_resnames (numpy.ndarray):
                Residue name for each PE atom, shape (N_pe,).
            - number_residues_pe (int):
                Number of residues in the PE region.
            - npe_coords (numpy.ndarray):
                NPE region Cartesian coordinates, shape (N_npe, 3), in Angstrom.
            - npe_atom_names (numpy.ndarray):
                Atom names for each NPE atom, shape (N_npe,).
            - npe_elements (numpy.ndarray):
                Element symbols for each NPE atom, shape (N_npe,).
            - npe_resids (numpy.ndarray):
                Residue id for each NPE atom, shape (N_npe,).
            - npe_resindices (numpy.ndarray):
                MDAnalysis residue index for each NPE atom, shape (N_npe,).
                This is unique per residue in the Universe and is used
                internally to keep PE/NPE residue partitioning robust even
                when `resid` values are reused across chains.
            - npe_resnames (numpy.ndarray):
                Residue name for each NPE atom, shape (N_npe,).
            - number_residues_npe (int):
                Number of residues in the NPE region.
        """
        return self._collect_structures(
            trajectory_file=trajectory_file,
            num_snapshots=num_snapshots,
            qm_region=qm_region,
            qm_charge=qm_charge,
            qm_multiplicity=qm_multiplicity,
            env_region=env_region,
            topology_file=topology_file,
            pe_cutoff=pe_cutoff,
            npe_cutoff=npe_cutoff,
            start=start,
            end=end,
            last_snapshot_only=last_snapshot_only,
            pdb_file=pdb_file,
            sampling_times=sampling_times,
            qm_molecule=qm_molecule)

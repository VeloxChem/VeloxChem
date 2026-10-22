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


class EmbeddingParser(TrajectoryParserBase):
    """
    Parses a single structure for embedding calculations.

    This class extracts the QM subsystem and the polarizable (PE) and
    non-polarizable (NPE) environment residues within the requested cutoff
    radii from a PDB file or a single trajectory frame.  It can be used for
    single-point embedding calculations without preparing a potfile.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

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
                   last_snapshot_only: bool = True,
                   pdb_file: str | None = None,
                   sampling_times=None,
                   qm_molecule=None):
        """
        Parses a single structure and returns a list containing one snapshot.

        The arguments are the same as for :meth:`EnsembleParser.structures`.
        By default only the last frame is returned.  Pass a ``sampling_times``
        sequence (or set ``last_snapshot_only=False`` together with
        ``num_snapshots=1``) to select a particular frame.
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

    def parse(self, *args, **kwargs):
        """
        Parses one structure and returns the snapshot dictionary.

        This is a convenience wrapper around :meth:`structures`.
        """
        kwargs.setdefault('last_snapshot_only', True)
        return self.structures(*args, **kwargs)[0]

    def parse_snapshot(self, *args, **kwargs):
        """
        Alias for :meth:`parse`.
        """
        return self.parse(*args, **kwargs)

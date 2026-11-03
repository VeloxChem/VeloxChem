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

import numpy as np

from .trajectoryparserbase import TrajectoryParserBase
from .embeddingpotfilewriter import EmbeddingPotfileWriter
from .fikaembedding import fika_embedding_sanity_check
from .inputparser import InputParser
from .molecule import Molecule
from .veloxchemlib import bohr_in_angstrom


class EmbeddingParser(TrajectoryParserBase):
    """
    Parses a single structure for embedding calculations.

    This class extracts the QM subsystem and the polarizable (PE) and
    non-polarizable (NPE) environment residues within the requested cutoff
    radii from a PDB file or a single trajectory frame.  It can be used for
    single-point embedding calculations and potfile preparation.

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

    def write_pot_file(self, snapshot, filename, pe_model=None, npe_model=None):
        """Write a parsed embedding snapshot to a potential file.

        :param snapshot:
            Snapshot dictionary returned by :meth:`parse`.
        :param filename:
            Output potential-file path.
        :param pe_model:
            PE model name or sequence of names, such as ``"SEP"`` or
            ``["CP3", "SEP"]``.
        :param npe_model:
            NPE model name or sequence of names, such as ``"tip3p"`` or
            ``["ff19sb", "tip3p"]``.
        :return:
            The output path, or ``None`` when the snapshot has no environment.
        """
        writer = EmbeddingPotfileWriter(self.comm, self.ostream)
        writer.set_env_models(pe_model=pe_model, npe_model=npe_model)
        return writer.write_pot_file(snapshot, filename)

    def read_pdb(self,
                 pdb_file,
                 *,
                 potfile,
                 qm_region="",
                 qm_charge=0,
                 qm_multiplicity=1,
                 env_region=None,
                 pe_cutoff=None,
                 npe_cutoff=None,
                 pe_model=("CP3", "SEP"),
                 npe_model=("ff19sb", "tip3p"),
                 embedding_settings=None):
        """Read a single-frame PDB and prepare a potfile-backed embedding.

        :param pdb_file:
            Input PDB file containing the QM subsystem and environment.
        :param potfile:
            Path of the potential file to create.
        :param qm_region:
            MDAnalysis selection for the QM subsystem. When omitted, the
            parser selects atoms with residue name ``MOL``.
        :param qm_charge:
            Total charge of the QM subsystem.
        :param qm_multiplicity:
            Spin multiplicity of the QM subsystem.
        :param env_region:
            Optional MDAnalysis selection restricting the environment.
        :param pe_cutoff:
            PE cutoff in angstrom.
        :param npe_cutoff:
            NPE cutoff in angstrom.
        :param pe_model:
            PE model name or sequence of names. Defaults to CP3 and SEP.
        :param npe_model:
            NPE model name or sequence of names. Defaults to ff19SB and TIP3P.
        :param embedding_settings:
            Optional Fika settings to add to the embedding dictionary. The
            embedding method is always ``"fika"``.
        :return:
            A tuple containing the QM :class:`Molecule` and the embedding
            dictionary suitable for ``driver.embedding``.
        """
        if Path(pdb_file).suffix.lower() != ".pdb":
            raise ValueError(
                "EmbeddingParser.read_pdb requires a .pdb file")

        if pe_cutoff is None and npe_cutoff is None:
            raise ValueError(
                "At least one of pe_cutoff or npe_cutoff must be provided.")

        snapshots = self.structures(
            pdb_file=pdb_file,
            qm_region=qm_region,
            qm_charge=qm_charge,
            qm_multiplicity=qm_multiplicity,
            env_region=env_region,
            pe_cutoff=pe_cutoff,
            npe_cutoff=npe_cutoff,
            last_snapshot_only=False,
        )
        if len(snapshots) != 1:
            raise ValueError(
                "EmbeddingParser.read_pdb requires a single-frame PDB file; "
                f"found {len(snapshots)} frames.")

        snapshot = snapshots[0]
        molecule = Molecule(snapshot["qm_elements"], snapshot["qm_coords"],
                            "angstrom")
        molecule.set_charge(int(snapshot["qm_charge"]))
        molecule.set_multiplicity(int(snapshot["qm_multiplicity"]))
        if not molecule.check_multiplicity():
            raise ValueError(
                "Incompatible QM charge and multiplicity for the selected "
                "QM subsystem.")

        settings = dict(embedding_settings or {})
        method = str(settings.setdefault("embedding_method", "fika")).lower()
        if method != "fika":
            raise ValueError("embedding_method must be 'fika'.")
        settings["embedding_method"] = "fika"
        embedding = {
            "settings": settings,
            "inputs": {
                "potfile": str(potfile),
            },
        }
        fika_embedding_sanity_check(embedding)

        potfile_path = self.write_pot_file(
            snapshot,
            potfile,
            pe_model=pe_model,
            npe_model=npe_model,
        )
        if potfile_path is None or not potfile_path.is_file():
            raise ValueError(
                "The selected PDB environment is empty; no potfile was "
                "created.")

        embedding["inputs"]["potfile"] = str(potfile_path)
        return molecule, embedding

    @staticmethod
    def _environment_xyz_from_potfile(potfile):
        """Build an XYZ string from a potfile ``@environment`` section."""
        data = InputParser(potfile).input_dict
        environment = data.get("environment")
        if not isinstance(environment, dict) or "xyz" not in environment:
            raise ValueError(
                f"{potfile}: missing @environment section with xyz")

        units = str(environment.get("units", "angstrom")).lower()
        if len(units) >= 3 and units == "angstrom"[:len(units)]:
            scale = 1.0
        elif units in ("au", "bohr"):
            scale = bohr_in_angstrom()
        else:
            raise ValueError(f"{potfile}: unknown units {units}")

        xyz_lines = environment["xyz"]
        if isinstance(xyz_lines, str):
            xyz_lines = [xyz_lines]

        atoms = []
        for line_number, line in enumerate(xyz_lines, start=1):
            fields = line.split()
            if len(fields) not in (6, 7):
                raise ValueError(
                    f"{potfile}: xyz line {line_number}: needs element, "
                    "three coordinates, residue, resid and atom name")
            try:
                coordinates = [float(value) * scale for value in fields[1:4]]
            except ValueError:
                raise ValueError(
                    f"{potfile}: xyz line {line_number}: invalid coordinates")
            if not all(np.isfinite(coordinates)):
                raise ValueError(
                    f"{potfile}: xyz line {line_number}: invalid coordinates")
            atoms.append((fields[0], coordinates))

        output = [str(len(atoms)), "embedding environment"]
        for element, (x, y, z) in atoms:
            output.append(f"{element} {x:.10f} {y:.10f} {z:.10f}")
        return "\n".join(output) + "\n"

    def show_embedding(self,
                       molecule,
                       potfile=None,
                       embedding=None,
                       width=600,
                       height=500,
                       mm_opacity=0.8):
        """Display a molecule and its potfile-defined environment.

        Exactly one of ``potfile`` and ``embedding`` must be supplied. When an
        embedding dictionary is given, its ``inputs.potfile`` entry is used.
        The potfile is the sole source of environment atoms and coordinates.

        :param molecule:
            QM molecule to display.
        :param potfile:
            Potential file containing an ``@environment`` section.
        :param embedding:
            Fika embedding dictionary containing ``inputs.potfile``.
        :param width:
            Viewer width in pixels.
        :param height:
            Viewer height in pixels.
        :param mm_opacity:
            Opacity of the environment atoms and bonds.
        :return:
            None.
        """
        if (potfile is None) == (embedding is None):
            raise ValueError(
                "show_embedding expects exactly one of potfile or embedding")

        if embedding is not None:
            if not isinstance(embedding, dict):
                raise TypeError("embedding must be a dictionary")
            inputs = embedding.get("inputs")
            if not isinstance(inputs, dict) or not inputs.get("potfile"):
                raise ValueError(
                    "show_embedding: embedding does not contain inputs.potfile"
                )
            potfile = inputs["potfile"]

        environment_xyz = self._environment_xyz_from_potfile(potfile)

        try:
            import py3Dmol
        except ImportError:
            raise ImportError("Unable to import py3Dmol")

        viewer = py3Dmol.view(width=width, height=height)
        viewer.setViewStyle({"style": "outline", "width": 0.02})

        viewer.addModel(molecule.get_xyz_string())
        viewer.setStyle({"model": 0}, {
            "stick": {},
            "sphere": {
                "scale": 0.25
            }
        })

        viewer.addModel(environment_xyz)
        viewer.setStyle(
            {"model": 1},
            {
                "stick": {
                    "radius": 0.1,
                    "opacity": mm_opacity
                },
                "sphere": {
                    "scale": 0.08,
                    "opacity": mm_opacity
                },
            },
        )
        viewer.zoomTo()
        viewer.show()

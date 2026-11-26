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

import csv
from io import StringIO
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

from mpi4py import MPI
import numpy as np

from .environment import get_data_path
from .errorhandler import assert_msg_critical
from .outputstream import OutputStream
from .veloxchemlib import mpi_master


class EmbeddingPotfileWriter:
    """Write embedding snapshots to Fika ``.pot`` files.

    The writer owns the PE/NPE parameter databases and the model selection
    needed to turn snapshots produced by :class:`EmbeddingParser` or
    :class:`EnsembleParser` into potential files.
    """

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD
        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.ostream = ostream

        db_dir = get_data_path() / "environment_parameters"
        if self.rank == mpi_master():
            sep_db = self._load_pe_db(db_dir / "pe_sep.csv")
            cp3_db = self._load_pe_db(db_dir / "pe_cp3.csv")
            tip3p_db = self._build_tip3p_db()
            ff19sb_db = self._build_ff19sb_db(
                db_dir / "protein.ff19SB.xml",
                db_dir / "ff19SB_aliases.csv",
            )
        else:
            sep_db = None
            cp3_db = None
            tip3p_db = None
            ff19sb_db = None

        self.pe_models = {
            "SEP": {
                "db": self.comm.bcast(sep_db, root=mpi_master())
            },
            "CP3": {
                "db": self.comm.bcast(cp3_db, root=mpi_master())
            },
        }
        self.npe_models = {
            "tip3p": {
                "db": self.comm.bcast(tip3p_db, root=mpi_master())
            },
            "ff19sb": {
                "db": self.comm.bcast(ff19sb_db, root=mpi_master())
            },
        }
        self.pe_model = None
        self.npe_model = None

    def _load_pe_db(self, csv_path: Path) -> dict:
        """Load a SEP or CP3 parameter database."""
        txt_path = csv_path.with_suffix(".txt")
        if not csv_path.is_file():
            raise FileNotFoundError(f"Parameter file not found: {csv_path}")
        if not txt_path.is_file():
            raise FileNotFoundError(f"Parameter file not found: {txt_path}")

        db = {}
        with (csv_path.open("r", newline="") as
              csv_file, txt_path.open("r", newline="") as txt_file):
            while True:
                line = csv_file.readline()
                if line.startswith("----") and line.rstrip().endswith("----"):
                    break

            for row, txt_row in zip(csv.DictReader(csv_file),
                                    csv.DictReader(txt_file)):
                resname = str(row["RESNAME"]).strip()
                atom_name = str(row["ATOMNAME"]).strip()
                if not resname or not atom_name:
                    continue

                assert_msg_critical(
                    resname == str(txt_row["res_name"]).strip() and
                    atom_name == str(txt_row["atom_name"]).strip(),
                    f"{type(self).__name__}: Inconsistent residue name or "
                    "atom name in parameter files",
                )
                db.setdefault(resname, {})[atom_name] = {
                    "element": str(txt_row["element"]).strip(),
                    "charge": float(row["q"]),
                    "polar": [
                        float(row[key])
                        for key in ("axx", "axy", "axz", "ayy", "ayz", "azz")
                    ],
                }
        return db

    @staticmethod
    def _build_ff19sb_db(ff19sb_params_path, ff19sb_aliases_path):
        """Build the ff19SB charge database."""
        name_aliases = {}
        with ff19sb_aliases_path.open("r", newline="") as aliases_file:
            for row in csv.DictReader(aliases_file):
                resname = str(row["res_name"]).strip()
                name_aliases.setdefault(resname, {})[str(
                    row["alternative_name"]).strip()] = str(
                        row["original_name"]).strip()

        root = ET.parse(str(ff19sb_params_path)).getroot()
        atom_types = {}
        residues = {}
        for node in root:
            if node.tag == "AtomTypes":
                for atom_type in node:
                    if atom_type.tag == "Type":
                        atom_types[atom_type.
                                   attrib["name"]] = atom_type.attrib["element"]
            elif node.tag == "Residues":
                for residue in node:
                    if residue.tag == "Residue":
                        residues[residue.attrib["name"]] = [
                            atom.attrib
                            for atom in residue
                            if atom.tag == "Atom"
                        ]

        db = {}
        for resname, atoms in residues.items():
            residue_db = db.setdefault(resname, {})
            for atom in atoms:
                residue_db[atom["name"]] = {
                    "element": atom_types[atom["type"]],
                    "charge": float(atom["charge"]),
                }

        for resname, aliases in name_aliases.items():
            if resname in db:
                for alternative, original in aliases.items():
                    db[resname][alternative] = db[resname][original]
        return db

    @staticmethod
    def _build_tip3p_db():
        """Build TIP3P water and ion charge parameters."""
        water = {
            "OW": {
                "element": "O",
                "charge": -0.834
            },
            "HW": {
                "element": "H",
                "charge": 0.417
            },
            "HW1": {
                "element": "H",
                "charge": 0.417
            },
            "HW2": {
                "element": "H",
                "charge": 0.417
            },
        }
        db = {name: water for name in ("WAT", "HOH", "SOL")}
        ions = {
            "Na": {
                "element": "Na",
                "charge": 1.0
            },
            "K": {
                "element": "K",
                "charge": 1.0
            },
            "Ca": {
                "element": "Ca",
                "charge": 2.0
            },
            "Mg": {
                "element": "Mg",
                "charge": 2.0
            },
            "Zn": {
                "element": "Zn",
                "charge": 2.0
            },
            "Cl": {
                "element": "Cl",
                "charge": -1.0
            },
            "Br": {
                "element": "Cl",
                "charge": -1.0
            },
        }
        for name, params in ions.items():
            db[name] = {name: params}
            db[name.upper()] = {name.upper(): params}
            suffix = "+" if params["charge"] > 0 else "-"
            if abs(params["charge"]) == 2.0:
                suffix = "2" + suffix
            ion_name = name + suffix
            db[ion_name] = {ion_name: params}
        return db

    @staticmethod
    def _normalize_model_names(model_names):
        if model_names is None:
            return []
        if isinstance(model_names, str):
            return [model_names]
        return list(model_names)

    @staticmethod
    def _merge_model_dbs(model_names, available_models,
                         model_kind: str) -> dict:
        merged_db = {}
        for model_name in model_names:
            if model_name not in available_models:
                raise KeyError(f"Unknown {model_kind} model '{model_name}'. \n"
                               f"Available: {sorted(available_models)}")
            for resname, atom_db in available_models[model_name]["db"].items():
                merged_residue = merged_db.setdefault(resname, {})
                for atom_name, params in atom_db.items():
                    if (atom_name in merged_residue and
                            merged_residue[atom_name] != params):
                        raise ValueError(
                            f"Conflicting {model_kind} parameters for "
                            f"{resname}/{atom_name} when combining models "
                            f"{model_names}")
                    merged_residue.setdefault(atom_name, dict(params))
        return merged_db

    def set_env_models(
            self,
            pe_model: str | list[str] | tuple[str, ...] | None = None,
            npe_model: str | list[str] | tuple[str, ...] | None = None):
        """Select the PE and/or NPE parameter models used for writing."""
        if pe_model is None and npe_model is None:
            raise ValueError(
                "At least one of pe_model or npe_model must be provided.")

        pe_names = self._normalize_model_names(pe_model)
        npe_names = self._normalize_model_names(npe_model)
        self.pe_model = ({
            "db": self._merge_model_dbs(pe_names, self.pe_models, "PE"),
            "model_names": pe_names,
        } if pe_names else None)
        self.npe_model = ({
            "db": self._merge_model_dbs(npe_names, self.npe_models, "NPE"),
            "model_names": npe_names,
        } if npe_names else None)

    @staticmethod
    def _first_residue_atom_pattern(atom_names, residue_ids, resnames,
                                    target_resname: str) -> list[str]:
        atom_names = np.asarray(atom_names, dtype=object)
        residue_ids = np.asarray(residue_ids, dtype=int)
        resnames = np.asarray(resnames, dtype=object)
        mask = (resnames == target_resname)
        if not np.any(mask):
            return []
        first_residue_id = int(residue_ids[mask][0])
        indices = np.where(mask & (residue_ids == first_residue_id))[0]
        return [str(atom_names[index]) for index in indices]

    @staticmethod
    def _normalize_resname_for_npe_db(db: dict, resname: str) -> str:
        resname = str(resname)
        if resname in db:
            return resname
        aliases = {
            "TIP3": "WAT",
            "TIP3P": "WAT",
            "HSD": "HID",
            "HSE": "HIE",
            "HSP": "HIP",
            "NHSD": "NHID",
            "NHSE": "NHIE",
            "NHSP": "NHIP",
            "CHSD": "CHID",
            "CHSE": "CHIE",
            "CHSP": "CHIP",
            "HYP": "HIP",
            "NHYP": "NHIP",
            "CHYP": "CHIP",
            "CYSH": "CYS",
        }
        mapped = aliases.get(resname, resname)
        if mapped in db:
            return mapped
        if len(resname) > 3 and resname[0] in {"N", "C"}:
            candidate = resname[0] + aliases.get(resname[1:], resname[1:])
            if candidate in db:
                return candidate
        return resname

    @staticmethod
    def _resolve_atom_name_for_npe_db(atom_name: str, available_atoms,
                                      resname: str) -> str | None:
        atom_name = str(atom_name)
        available = {str(atom) for atom in available_atoms}
        if atom_name in available:
            return atom_name

        candidates = []
        if resname in {"TIP3", "TIP3P", "WAT", "HOH", "SOL"}:
            if atom_name in ("OH2", "O"):
                candidates.append("OW")
            elif atom_name == "H1":
                candidates.extend(("HW1", "HW"))
            elif atom_name == "H2":
                candidates.extend(("HW2", "HW"))
        if atom_name == "HN":
            candidates.append("H")
        elif atom_name == "H":
            candidates.append("HN")
        if atom_name.startswith("HT") and len(
                atom_name) == 3 and atom_name[2] in "123":
            candidates.append(f"H{atom_name[2]}")
        if atom_name.startswith("H") and len(
                atom_name) == 2 and atom_name[1] in "123":
            candidates.append(f"HT{atom_name[1]}")
        if atom_name in {"OT1", "OC1"}:
            candidates.append("O")
        if atom_name in {"OT2", "OC2"}:
            candidates.append("OXT")
        if atom_name == "OXT":
            candidates.append("OT2")
        if atom_name == "HG1":
            candidates.append("HG")
        return next((name for name in candidates if name in available), None)

    @staticmethod
    def _normalize_resname_for_pe_db(db: dict, resname: str) -> str:
        resname = str(resname)
        if resname in db:
            return resname
        mapped = {"TIP3": "WAT", "TIP3P": "WAT"}.get(resname, resname)
        return mapped if mapped in db else resname

    @staticmethod
    def _resolve_atom_name_for_pe_db(atom_name: str, available_atoms,
                                     resname: str) -> str | None:
        atom_name = str(atom_name)
        available = {str(atom) for atom in available_atoms}
        if atom_name in available:
            return atom_name

        candidates = []
        if resname in {"TIP3", "TIP3P", "WAT", "HOH", "SOL"}:
            if atom_name in ("OH2", "O"):
                candidates.append("OW")
            elif atom_name == "H1":
                candidates.extend(("HW1", "HW"))
            elif atom_name == "H2":
                candidates.extend(("HW2", "HW"))
        if atom_name == "HN":
            candidates.append("H")
        elif atom_name == "H":
            candidates.append("HN")
        if atom_name.startswith("HT") and len(
                atom_name) == 3 and atom_name[2] in "123":
            candidates.append(f"H{atom_name[2]}")
        if atom_name.startswith("H") and len(
                atom_name) == 2 and atom_name[1] in "123":
            candidates.append(f"HT{atom_name[1]}")

        terminal_aliases = {
            "O1": ("OT1", "OC1", "O"),
            "O2": ("OT2", "OC2", "OXT"),
            "OXT": ("OT2", "OC2", "O2"),
            "OT1": ("O1", "O"),
            "OC1": ("O1", "O"),
            "OT2": ("O2", "OXT"),
            "OC2": ("O2", "OXT"),
            "O": ("OT1",),
        }
        candidates.extend(terminal_aliases.get(atom_name, ()))
        if atom_name == "HG1":
            candidates.append("HG")
        elif atom_name == "HG":
            candidates.append("HG1")

        match = re.match(r"^H([A-Z]+)(\d)(\d)$", atom_name)
        if match:
            candidates.append(
                f"{match.group(3)}H{match.group(1)}{match.group(2)}")
        match = re.match(r"^([123])H([A-Z]+)(\d)$", atom_name)
        if match:
            candidates.append(
                f"H{match.group(2)}{match.group(3)}{match.group(1)}")
        match = re.match(r"^H([A-Z]+)([123])$", atom_name)
        if match:
            candidates.append(f"{match.group(2)}H{match.group(1)}")
        match = re.match(r"^([123])H([A-Z]+)$", atom_name)
        if match:
            candidates.append(f"H{match.group(2)}{match.group(1)}")
        match = re.match(r"^HA([12])$", atom_name)
        if match:
            candidates.append(f"{match.group(1)}HA")
        match = re.match(r"^([12])HA$", atom_name)
        if match:
            candidates.append(f"HA{match.group(1)}")

        seen = set()
        for candidate in candidates:
            if candidate in seen:
                continue
            seen.add(candidate)
            if candidate in available:
                return candidate
        return None

    def _get_npe_params(self, atom_name: str, resname: str) -> dict:
        if self.npe_model is None:
            raise RuntimeError(
                "NPE model is not set. Call set_env_models(npe_model=...).")
        db = self.npe_model["db"]
        db_resname = self._normalize_resname_for_npe_db(db, resname)
        if db_resname not in db:
            raise KeyError(f"No NPE parameters for residue name '{resname}' "
                           f"(normalized to '{db_resname}').")
        resolved = self._resolve_atom_name_for_npe_db(atom_name,
                                                      db[db_resname].keys(),
                                                      db_resname)
        if resolved is None:
            raise KeyError(f"No NPE charge for {resname}/{atom_name}. "
                           f"Available: {sorted(db[db_resname].keys())}")
        return db[db_resname][resolved]

    @staticmethod
    def _ensure_no_split_residues_between_pe_and_npe(snapshot: dict):
        pe_ids = np.asarray(snapshot.get("pe_resindices",
                                         snapshot.get("pe_resids", [])),
                            dtype=int)
        npe_ids = np.asarray(snapshot.get("npe_resindices",
                                          snapshot.get("npe_resids", [])),
                             dtype=int)
        if pe_ids.size == 0 or npe_ids.size == 0:
            return
        shared = np.intersect1d(np.unique(pe_ids), np.unique(npe_ids))
        if shared.size:
            frame = int(snapshot.get("frame", -1))
            raise ValueError(
                "Snapshot contains residues split between PE and NPE regions "
                f"(frame {frame}, shared residue ids: {shared.tolist()}). "
                "PE/NPE partitioning must be residue-exclusive.")

    @classmethod
    def _validate_residue_atom_patterns(cls, atom_names, residue_ids, resnames,
                                        region: str):
        atom_names = np.asarray(atom_names, dtype=object)
        residue_ids = np.asarray(residue_ids, dtype=int)
        resnames = np.asarray(resnames, dtype=object)
        for resname in np.unique(resnames):
            mask = (resnames == resname)
            unique_ids = np.unique(residue_ids[mask])
            if unique_ids.size <= 1:
                continue
            reference = cls._first_residue_atom_pattern(atom_names, residue_ids,
                                                        resnames, str(resname))
            for residue_id in unique_ids[1:]:
                indices = np.where(mask & (residue_ids == residue_id))[0]
                current = [str(atom_names[index]) for index in indices]
                if current != reference:
                    raise ValueError(
                        f"Inconsistent {region} atom pattern detected for "
                        f"residue name '{resname}'. Residue id "
                        f"{int(unique_ids[0])} has pattern {reference}, while "
                        f"residue id {int(residue_id)} has pattern {current}. "
                        f"{region} residues with the same residue name must "
                        "contain the same ordered atom list.")

    @staticmethod
    def _snapshot_arrays(snapshot):
        arrays = {}
        for region in ("pe", "npe"):
            arrays[region] = {
                "coords": np.asarray(snapshot.get(f"{region}_coords", []),
                                     dtype=float),
                "elements": np.asarray(snapshot.get(f"{region}_elements", []),
                                       dtype=object),
                "resids": np.asarray(snapshot.get(f"{region}_resids", []),
                                     dtype=int),
                "resindices": np.asarray(snapshot.get(
                    f"{region}_resindices",
                    snapshot.get(f"{region}_resids", [])),
                                         dtype=int),
                "resnames": np.asarray(snapshot.get(f"{region}_resnames", []),
                                       dtype=object),
                "atom_names": np.asarray(snapshot.get(f"{region}_atom_names",
                                                      []),
                                         dtype=object),
            }
        return arrays

    def validate_snapshot(self, snapshot: dict):
        """Validate one parsed snapshot before writing."""
        arrays = self._snapshot_arrays(snapshot)
        has_pe = arrays["pe"]["coords"].size > 0
        has_npe = arrays["npe"]["coords"].size > 0
        if has_pe and self.pe_model is None:
            raise RuntimeError(
                "Snapshots contain PE atoms but pe_model is not set. "
                "Call set_env_models(pe_model=...).")
        if has_npe and self.npe_model is None:
            raise RuntimeError(
                "Snapshots contain NPE atoms but npe_model is not set. "
                "Call set_env_models(npe_model=...).")

        for region in ("pe", "npe"):
            data = arrays[region]
            if data["coords"].size == 0:
                continue
            atom_count = data["coords"].shape[0]
            for key in ("atom_names", "elements", "resids", "resnames",
                        "resindices"):
                if data[key].size != atom_count:
                    raise ValueError(
                        f"{region}_{key} is missing or wrong length in snapshots."
                    )
            self._validate_residue_atom_patterns(data["atom_names"],
                                                 data["resindices"],
                                                 data["resnames"],
                                                 region.upper())
        self._ensure_no_split_residues_between_pe_and_npe(snapshot)
        return arrays

    def validate_partition(self, snapshot: dict):
        """Ensure that PE and NPE regions do not split a residue."""
        self._ensure_no_split_residues_between_pe_and_npe(snapshot)

    @staticmethod
    def _source_value(value):
        if isinstance(value, (list, tuple)):
            value = ", ".join(str(item) for item in value)
        value = str(value)
        if "\n" in value or "\r" in value:
            raise ValueError("potfile source values must be single-line")
        return value

    def _write_source(self, output, snapshot, source):
        provenance = dict(snapshot.get("provenance", {}))
        if source is not None:
            if not isinstance(source, dict):
                raise TypeError(
                    "potfile source provenance must be a dictionary")
            provenance.update(source)

        provenance.setdefault("schema_version", 1)
        provenance.setdefault("generator", "VeloxChem EmbeddingPotfileWriter")
        provenance.setdefault("frame", snapshot.get("frame"))
        provenance.setdefault("trajectory_time_ps",
                              snapshot.get("trajectory_time"))
        provenance.setdefault("qm_charge", snapshot.get("qm_charge"))
        provenance.setdefault("qm_multiplicity",
                              snapshot.get("qm_multiplicity"))
        if self.pe_model is not None:
            provenance["pe_models"] = self.pe_model["model_names"]
        if self.npe_model is not None:
            provenance["npe_models"] = self.npe_model["model_names"]

        preferred_order = (
            "schema_version", "generator", "pdb_file", "topology_file",
            "trajectory_file", "frame", "trajectory_time_ps", "qm_region",
            "qm_charge", "qm_multiplicity", "environment_region",
            "pe_cutoff_angstrom", "npe_cutoff_angstrom", "pe_models",
            "npe_models",
        )
        keys = [key for key in preferred_order if key in provenance]
        keys.extend(sorted(set(provenance) - set(keys)))

        output.write("@source\n")
        for key in keys:
            value = provenance[key]
            if value is None or (isinstance(value, str) and not value):
                continue
            if isinstance(value, (list, tuple)) and not value:
                continue
            output.write(f"{key}: {self._source_value(value)}\n")
        output.write("@end\n\n")

    def render_potfile(self, snapshot: dict, source=None) -> str | None:
        """Render a snapshot as potfile text, including provenance."""

        arrays = self.validate_snapshot(snapshot)
        pe = arrays["pe"]
        npe = arrays["npe"]
        if pe["coords"].size == 0 and npe["coords"].size == 0:
            return None

        unique_resnames = {}
        for region, data in (("pe", pe), ("npe", npe)):
            seen = set()
            unique_resnames[region] = []
            for name in map(str, data["resnames"].tolist()):
                if name not in seen:
                    seen.add(name)
                    unique_resnames[region].append(name)

        output = StringIO()
        self._write_source(output, snapshot, source)
        output.write("@environment\nunits: angstrom\nxyz:\n")
        for region, data in (("pe", pe), ("npe", npe)):
            for coords, element, resname, resid, atom_name in zip(
                    data["coords"], data["elements"], data["resnames"],
                    data["resids"], data["atom_names"]):
                x, y, z = coords
                output.write(
                    f"{str(element):<2} {x:12.6f} {y:12.6f} {z:12.6f}  "
                    f"{str(resname)}_{region}  {int(resid)}  "
                    f"{str(atom_name)}\n")
        output.write("@end\n\n@charges\n")

        pe_db = None if self.pe_model is None else self.pe_model["db"]
        for resname in unique_resnames["pe"]:
            db_resname = self._normalize_resname_for_pe_db(pe_db, resname)
            if db_resname not in pe_db:
                raise KeyError(f"No PE parameters for residue name '{resname}'")
            pattern = self._first_residue_atom_pattern(pe["atom_names"],
                                                       pe["resindices"],
                                                       pe["resnames"], resname)
            for atom_name in pattern:
                resolved = self._resolve_atom_name_for_pe_db(
                    atom_name, pe_db[db_resname].keys(), db_resname)
                if resolved is None:
                    raise KeyError(
                        f"No PE params for {resname}/{atom_name}. "
                        f"Available: {sorted(pe_db[db_resname].keys())}")
                params = pe_db[db_resname][resolved]
                output.write(
                    f"{params['element']:<2} {params['charge']:12.8f}  "
                    f"{resname}_pe\n")

        for resname in unique_resnames["npe"]:
            pattern = self._first_residue_atom_pattern(npe["atom_names"],
                                                       npe["resindices"],
                                                       npe["resnames"], resname)
            for atom_name in pattern:
                params = self._get_npe_params(atom_name, resname)
                output.write(
                    f"{params['element']:<2} {params['charge']:12.8f}  "
                    f"{resname}_npe\n")
        output.write("@end\n\n")

        if pe["coords"].size:
            output.write("@polarizabilities\n")
            for resname in unique_resnames["pe"]:
                db_resname = self._normalize_resname_for_pe_db(pe_db, resname)
                pattern = self._first_residue_atom_pattern(
                    pe["atom_names"], pe["resindices"], pe["resnames"], resname)
                for atom_name in pattern:
                    resolved = self._resolve_atom_name_for_pe_db(
                        atom_name, pe_db[db_resname].keys(), db_resname)
                    if resolved is None:
                        raise KeyError(
                            f"No PE params for {resname}/{atom_name}. "
                            f"Available: {sorted(pe_db[db_resname].keys())}")
                    params = pe_db[db_resname][resolved]
                    polar = params["polar"]
                    output.write(f"{params['element']:<2} {polar[0]:12.8f} "
                                 f"{polar[1]:12.8f} {polar[2]:12.8f} "
                                 f"{polar[3]:12.8f} {polar[4]:12.8f} "
                                 f"{polar[5]:12.8f}  {resname}_pe\n")
            output.write("@end\n")
        return output.getvalue()

    def _potfile_text(self, snapshot: dict) -> str | None:
        """Backward-compatible alias for :meth:`render_potfile`."""

        return self.render_potfile(snapshot)

    def write_pot_file(self,
                       snapshot: dict,
                       filename: str | Path,
                       source=None):
        """Write one snapshot to an explicitly named potential file."""
        text = self.render_potfile(snapshot, source=source)
        path = Path(filename)
        if text is not None and self.rank == mpi_master():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        self.comm.barrier()
        return path if text is not None else None

    def write_pot_files(self, snapshots, outdir: str | Path = "pot_frames"):
        """Write one frame-numbered potential file per environment snapshot."""
        if isinstance(snapshots, dict):
            snapshots = [snapshots]
        paths = []
        outdir = Path(outdir)
        for snapshot in snapshots:
            frame = int(snapshot["frame"])
            path = self.write_pot_file(snapshot,
                                       outdir / f"pe_frame_{frame:06d}.pot")
            if path is not None:
                paths.append(path)
        return paths

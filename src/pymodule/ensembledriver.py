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
import numpy as np
import xml.etree.ElementTree as ET
import sys
import csv
import re

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .lreigensolver import LinearResponseEigenSolver
from .cppsolver import ComplexResponseSolver
from .spectrumaverager import SpectrumAverager
from .environment import get_data_path
from .errorhandler import assert_msg_critical


class EnsembleDriver:
    """
    Implements ensemble calculations over a list of snapshots produced by
    :class EnsembleParser.

    Each snapshot contains teh QM subsysem and (optionally) environment atoms.
    The driver can run, for each snapshot:

    - A SCF calculation (restricted or unrestricted)
    - Polarizable embedding (PE) via per-snapshot .pot files.
    - Non-polarizable embedding (NPE) via point charges.
    - Optional linear response properties using
      :class LinearResponseEigenSolver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - pe_models: Dictionary of available PE models (e.g., SEP, CP3).
        - npe_models: Dictionary of available NPE models (e.g., tip3p, ff19sb).
        - pe_model: Selected PE model dictctionary.
        - npe_model: Selected NPE model dictionary.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - excited_states: If True, run linear response after SCF.
          response is run as well when any response option is set (e.g., :attr: `nstates`).
        - nstates: Number of excited states to compute in linear response.
        - nto: Whether to run NOT analysis in linear response.
        - core_excitation: Whether to compute core excitations in linear response.
        - num_core_orbitals: Number of core orbitals to consider for core excitations in linear response.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the ensemble driver.

        PE (SEP/CP3) and NPE (TIP3P/ff19sb) parameter tables are read from the
        ``database/environment_parameters`` directory.

        See this reference, Figure 4, for a summary of 
        an overview of SEP/CP3 parametrizations:
        https://doi.org/10.1021/acs.jctc.5c01719
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

        db_dir = get_data_path() / "environment_parameters"

        sep_parameters_file = db_dir / "pe_sep.csv"
        cp3_parameters_file = db_dir / "pe_cp3.csv"
        ff19sb_parameters_file = db_dir / "protein.ff19SB.xml"
        ff19sb_aliases_file = db_dir / "ff19SB_aliases.csv"

        if self.rank == mpi_master():
            self._sep_db = self._load_pe_db(sep_parameters_file)
            self._cp3_db = self._load_pe_db(cp3_parameters_file)
            self._tip3p_db = self._build_tip3p_db()
            self._ff19sb_db = self._build_ff19sb_db(ff19sb_parameters_file,
                                                    ff19sb_aliases_file)
        else:
            self._sep_db = None
            self._cp3_db = None
            self._tip3p_db = None
            self._ff19sb_db = None

        self._sep_db = self.comm.bcast(self._sep_db, root=mpi_master())
        self._cp3_db = self.comm.bcast(self._cp3_db, root=mpi_master())
        self._tip3p_db = self.comm.bcast(self._tip3p_db, root=mpi_master())
        self._ff19sb_db = self.comm.bcast(self._ff19sb_db, root=mpi_master())

        self.pe_models = {
            "SEP": {"db": self._sep_db},
            "CP3": {"db": self._cp3_db},
        }
        self.npe_models = {
            "tip3p": {"db": self._tip3p_db},
            "ff19sb": {"db": self._ff19sb_db},
        }

        self.pe_model = None
        self.npe_model = None

    def _load_pe_db(self, csv_path: Path) -> dict:
        """
        Loads the raw PE parameter table:
            RESNAME,ATOMNAME,q,axx,axy,axz,ayy,ayz,azz

        :param csv_path:
            Path to the raw PE parameter file.

        :return:
            db[res_name][atom_name] = {"element": str, "charge": float, "polar": [6 floats]}
        """

        txt_path = csv_path.with_suffix(".txt")

        if not csv_path.is_file():
            raise FileNotFoundError(f"Parameter file not found: {csv_path}")
        if not txt_path.is_file():
            raise FileNotFoundError(f"Parameter file not found: {txt_path}")

        db = {}

        with csv_path.open("r", newline="") as fh, txt_path.open("r", newline="") as fh_txt:
            # Skip first few lines of comments, up to "----"
            while True:
                line = fh.readline()
                if line.startswith("----") and line.rstrip().endswith("----"):
                    break

            reader = csv.DictReader(fh)
            txt_reader = csv.DictReader(fh_txt)

            for row, txt_row in zip(reader, txt_reader):
                resn = str(row["RESNAME"]).strip()
                atom = str(row["ATOMNAME"]).strip()
                if not resn or not atom:
                    continue

                resn_from_txt = str(txt_row["res_name"]).strip()
                atom_from_txt = str(txt_row["atom_name"]).strip()
                assert_msg_critical(
                    resn == resn_from_txt and atom == atom_from_txt,
                    f'{type(self).__name__}: Inconsistent residue name or ' +
                    'atom name in parameter files')

                elem = str(txt_row["element"]).strip()

                q = float(row["q"])

                pol6 = [
                    float(row["axx"]),
                    float(row["axy"]),
                    float(row["axz"]),
                    float(row["ayy"]),
                    float(row["ayz"]),
                    float(row["azz"]),
                ]

                if resn not in db:
                    db[resn] = {}

                db[resn][atom] = {"element": elem, "charge": q, "polar": pol6}

        return db

    @staticmethod
    def _build_ff19sb_db(ff19sb_params_path, ff19sb_aliases_path):
        """
        Builds ff19SB parameters.

        :param ff19sb_params_path:
            The path to ff19SB parameters file.
        :param ff19sb_aliases_path:
            The path to ff19SB name aliases file.
        :return:
            The dictionary of ff19SB parameters.
        """

        name_aliases = {}
        with ff19sb_aliases_path.open("r", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                resn = str(row["res_name"]).strip()
                alt_name = str(row["alternative_name"]).strip()
                name = str(row["original_name"]).strip()
                if resn not in name_aliases:
                    name_aliases[resn] = {}
                name_aliases[resn][alt_name] = name

        tree = ET.parse(str(ff19sb_params_path))
        root = tree.getroot()

        data = {
            'version': None,
            'refs': [],
            'atom_types': {},
            'residues': {},
        }

        for node in root:

            if node.tag == 'Info':
                for subnode in node:
                    if subnode.tag == 'Source':
                        info = subnode.attrib
                        data['version'] = info['Source']
                    elif subnode.tag == 'Reference':
                        data['refs'].append(subnode.text)

            elif node.tag == 'AtomTypes':
                for subnode in node:
                    if subnode.tag == 'Type':
                        info = subnode.attrib
                        atomtype = info['name']
                        data['atom_types'][atomtype] = {'element': info['element']}

            elif node.tag == 'Residues':
                for subnode in node:
                    if subnode.tag == 'Residue':
                        info = subnode.attrib
                        resname = info['name']
                        res_atoms = []
                        for subsubnode in subnode:
                            if subsubnode.tag == 'Atom':
                                atom_info = subsubnode.attrib
                                res_atoms.append(atom_info)
                        data['residues'][resname] = res_atoms


        db = {}

        for resname in data['residues']:
            db[resname] = {}

            for atom in data['residues'][resname]:
                atomname = atom['name']
                db[resname][atomname] = {}

                atomtype = atom['type']
                element = data['atom_types'][atomtype]['element']

                db[resname][atomname]['element'] = element
                db[resname][atomname]['charge'] = float(atom['charge'])

        for resname in db:
            if resname in name_aliases:
                for alt_name, name in name_aliases[resname].items():
                    db[resname][alt_name] = db[resname][name]

        return db

    @staticmethod
    def _build_tip3p_db():
        """
        Builds tip3p parameters.

        :return:
            The tip3p parameters.
        """

        db = {}

        # W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein,
        # J. Chem. Phys. 79, 926-935 (1983)
        tip3p_params = {
            'OW': {'element': 'O', 'charge': -0.834},
            'HW': {'element': 'H', 'charge': 0.417},
            'HW1': {'element': 'H', 'charge': 0.417},
            'HW2': {'element': 'H', 'charge': 0.417},
       }

        db['WAT'] = tip3p_params
        db['HOH'] = tip3p_params
        db['SOL'] = tip3p_params

        ions_params = {
            'Na': {'element': 'Na', 'charge': 1.0},
            'K': {'element': 'K', 'charge': 1.0},
            'Ca': {'element': 'Ca', 'charge': 2.0},
            'Mg': {'element': 'Mg', 'charge': 2.0},
            'Zn': {'element': 'Zn', 'charge': 2.0},
            'Cl': {'element': 'Cl', 'charge': -1.0},
            'Br': {'element': 'Cl', 'charge': -1.0},
        }

        for key, val in ions_params.items():
            db[key] = {key: val}
            db[key.upper()] = {key.upper(): val}
            if val['charge'] == 2.0:
                db[f'{key}2+'] = {f'{key}2+': val}
            elif val['charge'] == 1.0:
                db[f'{key}+'] = {f'{key}+': val}
            elif val['charge'] == -1.0:
                db[f'{key}-'] = {f'{key}-': val}

        return db

    def _get_npe_params(self, atom_name: str, resname: str) -> dict:
        """
        Resolve an NPE residue/atom name and return the corresponding parameters.
        """
        if self.npe_model is None:
            raise RuntimeError("NPE model is not set. Call set_env_models(npe_model=...).")

        db = self.npe_model["db"]
        raw_resn = str(resname)
        raw_atom = str(atom_name)
        db_resn = self._normalize_resname_for_npe_db(db, raw_resn)

        if db_resn not in db:
            raise KeyError(
                f"No NPE parameters for residue name '{raw_resn}' "
                f"(normalized to '{db_resn}')."
            )

        resolved_atom = self._resolve_atom_name_for_npe_db(raw_atom, db[db_resn].keys())
        if resolved_atom is None:
            raise KeyError(
                f"No NPE charge for {raw_resn}/{raw_atom}. "
                f"Available: {sorted(db[db_resn].keys())}"
            )

        return db[db_resn][resolved_atom]

    @staticmethod
    def _normalize_model_names(model_names):
        """
        Normalize a model selection to a list of model names.
        Accepts a single string, a list/tuple of strings, or None.
        """
        if model_names is None:
            return []
        if isinstance(model_names, str):
            return [model_names]
        return list(model_names)

    @staticmethod
    def _merge_model_dbs(model_names, available_models, model_kind: str) -> dict:
        """
        Merge parameter databases from one or more environment models.

        :param model_names:
            List of model names to merge.
        :param available_models:
            Dictionary of available models.
        :param model_kind:
            Model kind used in error messages (e.g., "PE" or "NPE").

        :return:
            Merged parameter database.
        """
        merged_db = {}
        for model_name in model_names:
            if model_name not in available_models:
                raise KeyError(
                    f"Unknown {model_kind} model '{model_name}'. \nAvailable: {sorted(available_models)}"
                )
            model_db = available_models[model_name]["db"]
            for resn, atom_db in model_db.items():
                merged_res_db = merged_db.setdefault(resn, {})
                for atom, params in atom_db.items():
                    if atom in merged_res_db:
                        old = merged_res_db[atom]
                        if old != params:
                            raise ValueError(
                                f"Conflicting {model_kind} parameters for {resn}/{atom} "
                                f"when combining models {model_names}"
                            )
                    else:
                        merged_res_db[atom] = dict(params)
        return merged_db
        
    def set_env_models(self,
                       pe_model: str | list[str] | tuple[str, ...] | None = None,
                       npe_model: str | list[str] | tuple[str, ...] | None = None):
        """
        Set PE and/or NPE environment models

        Valid scenarios:
        PE models = CP3 (for protein), SEP (for water and ions).
        NPE models = ff19sb (for protein), tip3p (for water).
        User can select any combination:

        - PE only:  e.g., set_env_models(pe_model="CP3")
        - NPE only: e.g., set_env_models(npe_model="ff19sb")
        - Both:     e.g., set_env_models(pe_model="CP3", npe_model="ff19sb")
        - Multiple: e.g., set_env_models(pe_model=["SEP", "CP3"], 
                                         npe_model=["tip3p", "ff19sb"]),
                                         e.g. a system that contains a protein,
                                         water, and described with both pe and npe.

        :param pe_model:
            Name of the PE parameter model (e.g., "SEP", "CP3") or a list
            of PE model names.
        :param npe_model:
            Name of the NPE parameter model (e.g., "tip3p", "ff19sb") or a list
            of NPE model names.
        :return:
            None.
        :raises ValueError:
            If neither is provided or if selected models contain conflicting
            parameters for the same residue/atom entry.
        :raises KeyError:
            If an unknown model name is requested.
        """
        if pe_model is None and npe_model is None:
            raise ValueError("At least one of pe_model or npe_model must be provided.")
        
        pe_model_names = self._normalize_model_names(pe_model)
        npe_model_names = self._normalize_model_names(npe_model)

        if pe_model_names:
            self.pe_model = {
                "db": self._merge_model_dbs(pe_model_names, self.pe_models, "PE"),
                "model_names": pe_model_names,
            }
        else:
            self.pe_model = None

        if npe_model_names:
            self.npe_model = {
                "db": self._merge_model_dbs(npe_model_names, self.npe_models, "NPE"),
                "model_names": npe_model_names,
            }
        else:
            self.npe_model = None
    
    @staticmethod
    def _first_residue_atom_pattern(atom_names, residue_ids, resnames, target_resname: str) -> list[str]:
        """
        Return atom_name pattern for one residue instance of target_resname,
        preserving the order in the arrays.
        :param atom_names:
            The array of atom names.
        :param residue_ids:
            The MDAnalysis 'resindex' array, which is unique per residue in the Universe.
            'resid' values may be reused across chains and are therefore not robust enough for internal
+           residue identity checks.
        :param resnames:
            The array of residue names.
        :param target_resname:
            The name of the target residue.
        :return:
            A list of atom names for the first instance of the target residue.
        """
        atom_names = np.asarray(atom_names, dtype=object)
        residue_ids = np.asarray(residue_ids, dtype=int)
        resnames = np.asarray(resnames, dtype=object)

        mask = (resnames == target_resname)
        if not np.any(mask):
            return []

        first_residue_id = int(residue_ids[mask][0])
        idx = np.where(mask & (residue_ids == first_residue_id))[0]
        return [str(atom_names[i]) for i in idx]

    @staticmethod
    def _normalize_resname_for_npe_db(db: dict, resname: str) -> str:
        """
        Normalize residue name to match the selected NPE database.

        This is primarily used to bridge CHARMM-style residue naming (e.g. HSD/HSE/HSP)
        to AMBER-style naming (HID/HIE/HIP) for ff19sb.
        :param db:
            The NPE database to check against.
        :param resname:
            The residue name to normalize.
        :return:
            The normalized residue name.
        """
        resname = str(resname)

        # Fast path
        if resname in db:
            return resname

        # Common CHARMM -> AMBER mappings
        res_alias = {
            # Histidine tautomers / charge states
            "HSD": "HID",
            "HSE": "HIE",
            "HSP": "HIP",
            "NHSD": "NHID",
            "NHSE": "NHIE",
            "NHSP": "NHIP",
            "CHSD": "CHID",
            "CHSE": "CHIE",
            "CHSP": "CHIP",
            # Some workflows call protonated histidine HYP; map if present
            "HYP": "HIP",
            "NHYP": "NHIP",
            "CHYP": "CHIP",
            # Cysteine neutral thiol naming variants
            "CYSH": "CYS",
        }

        mapped = res_alias.get(resname, resname)
        if mapped in db:
            return mapped

        # If terminal prefixes exist, try mapping the core residue name
        # (e.g. N + HSD -> NHID) or (C + HSE -> CHIE)
        if len(resname) > 3 and (resname[0] in {"N", "C"}):
            prefix = resname[0]
            core = resname[1:]
            core_mapped = res_alias.get(core, core)
            cand = prefix + core_mapped
            if cand in db:
                return cand
        return resname

    @staticmethod
    def _resolve_atom_name_for_npe_db(atom_name: str, available_atoms) -> str | None:
        """
        Resolve CHARMM-style atom names to names present in the NPE database.

        :param atom_name:
            The atom name to resolve.
        :param available_atoms:
            The set of atom names available in the NPE database for the given residue.
        :return:
            The resolved atom name if found, or None if no match is found.
        """
        atom_name = str(atom_name)
        avail = set(str(a) for a in available_atoms)

        # Exact match first
        if atom_name in avail:
            return atom_name

        candidates: list[str] = []

        # Backbone amide proton: CHARMM often uses HN, AMBER often uses H (varies by residue in db)
        if atom_name == "HN":
            candidates.append("H")
        elif atom_name == "H":
            candidates.append("HN")

        # N-terminus ammonium hydrogens: HT1/HT2/HT3 (CHARMM) vs H1/H2/H3 (AMBER)
        if atom_name.startswith("HT") and len(atom_name) == 3 and atom_name[2] in "123":
            candidates.append(f"H{atom_name[2]}")
        if atom_name.startswith("H") and len(atom_name) == 2 and atom_name[1] in "123":
            candidates.append(f"HT{atom_name[1]}")

        # C-terminus oxygens: OT1/OT2 (CHARMM) vs O / OXT (AMBER)
        if atom_name in {"OT1", "OC1"}:
            candidates.append("O")
        if atom_name in {"OT2", "OC2"}:
            candidates.append("OXT")
        if atom_name == "OXT":
            candidates.append("OT2")

        # Cysteine thiol hydrogen: HG1 (CHARMM) vs HG (AMBER) for ff19sb cysteine
        if atom_name == "HG1":
            candidates.append("HG")

        # Try candidates in order
        for cand in candidates:
            if cand in avail:
                return cand

        return None


    @staticmethod
    def _resolve_atom_name_for_pe_db(atom_name: str, available_atoms) -> str | None:
        """
        Resolve atom names to names present in the PE database.

        This bridges naming differences coming from different topology conventions,
        e.g. terminal oxygens O1/O2 (CHARMM topologies) vs OT1/OT2 or OC1/OC2
        in CP3/SEP tables.

        The resolver tries:
        1) exact match
        2) common aliases (HN <-> H, HT1-3 <-> H1-3)
        3) terminal oxygen aliases (O1/O2 <-> OT1/OT2, OC1/OC2, O/OXT)
        4) hydrogen digit conventions (HB1 <-> 1HB, HD11 <-> 1HD1, etc.)
        5) cysteine thiol hydrogen naming alias (HG1 <-> HG)

        :param atom_name:
            Atom name from the trajectory/topology.
        :param available_atoms:
            Atom names available in the selected PE db for the residue.
        :return:
            A matching atom name in the db, or None.
        """

        atom_name = str(atom_name)
        avail = set(str(a) for a in available_atoms)

        # Exact match first
        if atom_name in avail:
            return atom_name

        candidates: list[str] = []

        # Backbone amide proton name
        if atom_name == "HN":
            candidates.append("H")
        elif atom_name == "H":
            candidates.append("HN")

        # N-terminus ammonium hydrogens: HT1/HT2/HT3 (CHARMM) vs H1/H2/H3 (AMBER/GROMACS)
        if atom_name.startswith("HT") and len(atom_name) == 3 and atom_name[2] in "123":
            candidates.append(f"H{atom_name[2]}")
        if atom_name.startswith("H") and len(atom_name) == 2 and atom_name[1] in "123":
            candidates.append(f"HT{atom_name[1]}")

        # Terminal carboxylate oxygens:
        # Many GROMACS topologies use O1/O2. CP3/SEP often use OT1/OT2 (and sometimes OC1/OC2).
        if atom_name == "O1":
            candidates.extend(["OT1", "OC1", "O"])
        elif atom_name == "O2":
            candidates.extend(["OT2", "OC2", "OXT"])
        elif atom_name == "OXT":
            candidates.extend(["OT2", "OC2", "O2"])
        elif atom_name in {"OT1", "OC1"}:
            candidates.append("O1")
        elif atom_name in {"OT2", "OC2"}:
            candidates.append("O2")

        # Some CHARMM conventions: OT1/OT2 <-> O/OXT (useful if db uses O/OXT instead)
        if atom_name in {"OT1", "OC1"}:
            candidates.append("O")
        if atom_name in {"OT2", "OC2"}:
            candidates.append("OXT")
        if atom_name == "O":
            candidates.append("OT1")
        if atom_name == "OXT":
            candidates.append("OT2")

        # Cysteine thiol hydrogen naming:
        # CHARMM/GROMAS often uses HG1, while CP3 entires use HG.
        if atom_name == "HG1":
            candidates.append("HG")
        elif atom_name == "HG":
            candidates.append("HG1")

        # Hydrogen digit conventions:
        #   - CHARMM: HD11, HD12, HD13  <->  AMBER: 1HD1, 2HD1, 3HD1
        #   - Also HB1 <-> 1HB, HG2 <-> 2HG, etc.
        m = re.match(r"^H([A-Z]+)(\d)(\d)$", atom_name)
        if m:
            base = m.group(1)
            carbon_idx = m.group(2)
            hyd_idx = m.group(3)
            candidates.append(f"{hyd_idx}H{base}{carbon_idx}")  # HD12 -> 2HD1

        m = re.match(r"^([123])H([A-Z]+)(\d)$", atom_name)
        if m:
            hyd_idx = m.group(1)
            base = m.group(2)
            carbon_idx = m.group(3)
            candidates.append(f"H{base}{carbon_idx}{hyd_idx}")  # 2HD1 -> HD12

        # Single-digit variants like HB1 <-> 1HB
        m = re.match(r"^H([A-Z]+)([123])$", atom_name)
        if m:
            base = m.group(1)
            idx = m.group(2)
            candidates.append(f"{idx}H{base}")  # HB1 -> 1HB

        m = re.match(r"^([123])H([A-Z]+)$", atom_name)
        if m:
            idx = m.group(1)
            base = m.group(2)
            candidates.append(f"H{base}{idx}")  # 1HB -> HB1

        # Glycine alpha hydrogens: HA1/HA2 <-> 1HA/2HA
        m = re.match(r"^HA([12])$", atom_name)
        if m:
            candidates.append(f"{m.group(1)}HA")
        m = re.match(r"^([12])HA$", atom_name)
        if m:
            candidates.append(f"HA{m.group(1)}")

        # Remove duplicates while preserving order
        seen = set()
        uniq_candidates = []
        for c in candidates:
            if c not in seen:
                seen.add(c)
                uniq_candidates.append(c)

        for cand in uniq_candidates:
            if cand in avail:
                return cand

        return None

    @staticmethod
    def _ensure_no_split_residues_between_pe_and_npe(snap: dict):
        """
        Ensure that no residue is split across PE and NPE regions.
        """
        pe_ids = np.asarray(
            snap.get("pe_resindices", snap.get("pe_resids", [])),
            dtype=int,
        )
        npe_ids = np.asarray(
            snap.get("npe_resindices", snap.get("npe_resids", [])),
            dtype=int,
        )

        if pe_ids.size == 0 or npe_ids.size == 0:
            return

        shared = np.intersect1d(np.unique(pe_ids), np.unique(npe_ids))
        if shared.size > 0:
            frame = int(snap.get("frame", -1))
            raise ValueError(
                "Snapshot contains residues split between PE and NPE regions "
                f"(frame {frame}, shared residue ids: {shared.tolist()}). "
                "PE/NPE partitioning must be residue-exclusive."
            )

    @classmethod
    def _validate_pe_residue_atom_patterns(cls, atom_names, residue_ids, resnames):
        """
        Validate that all PE residues sharing the same residue name also share the
        same ordered atom-name pattern.

        The current PE writer stores one @charges/@polarizabilities template per
        residue name. If two residue instances with the same resname contain
        different atom subsets, the generated .pot file becomes ambiguous and will
        fail later.

        :param cls:
            The class object (used for calling the _first_residue_atom_pattern method).
        :param atom_names:
            The array of atom names for the PE atoms.
        :param residue_ids:
            The array of residue ids for the PE atoms (e.g., MDAnalysis resindex).
        :param resnames:
            The array of residue names for the PE atoms.
        """
        atom_names = np.asarray(atom_names, dtype=object)
        residue_ids = np.asarray(residue_ids, dtype=int)
        resnames = np.asarray(resnames, dtype=object)

        if atom_names.size == 0:
            return

        for resn in np.unique(resnames):
            mask = (resnames == resn)
            unique_residue_ids = np.unique(residue_ids[mask])

            if unique_residue_ids.size <= 1:
                continue

            ref_pattern = cls._first_residue_atom_pattern(
                atom_names, residue_ids, resnames, str(resn)
            )

            for residue_id in unique_residue_ids[1:]:
                idx = np.where(mask & (residue_ids == residue_id))[0]
                cur_pattern = [str(atom_names[i]) for i in idx]
                if cur_pattern != ref_pattern:
                    raise ValueError(
                        "Inconsistent PE atom pattern detected for residue name "
                        f"'{resn}'. Residue id {int(unique_residue_ids[0])} has "
                        f"pattern {ref_pattern}, while residue id {int(residue_id)} "
                        f"has pattern {cur_pattern}. PE residues with the same "
                        "residue name must contain the same ordered atom list."
                    )


    def _build_point_charges(self, coords_ang, atom_names, resnames) -> np.ndarray | None:
        """
        Build point charges array expected by SCF driver: shape (6, N), coords in bohr.

        Charges are taken from the selected NPE database by (resname, atom_name).
        This routine performs lightweight CHARMM->AMBER normalization so that
        CHARMM-style names in the trajectory (e.g. HN, OT1/OT2, HSD/HSE/HSP)
        can be resolved against AMBER-style databases such as ff19sb.

        :param coords_ang:
            Coordinates of the NPE atoms in angstrom, shape (N, 3).
        :param atom_names:
            Atom names of the NPE atoms, shape (N,).
        :param resnames:
            Residue names of the NPE atoms, shape (N,).
        :return:
            Point charges array of shape (6, N) with coordinates in bohr and charges from the NPE database,
            or None if there are no NPE atoms.
        """
        coords_ang = np.asarray(coords_ang, dtype=float)
        if coords_ang.size == 0:
            return None

        if self.npe_model is None:
            raise RuntimeError("Snapshot contains NPE atoms but npe_model is not set.")

        atom_names = np.asarray(atom_names, dtype=object)
        resnames = np.asarray(resnames, dtype=object)

        q = np.empty(atom_names.size, dtype=float)

        for i, (atom, resn) in enumerate(zip(atom_names, resnames)):
            q[i] = self._get_npe_params(atom, resn)["charge"]

        pc = np.zeros((6, coords_ang.shape[0]), dtype=float)
        pc[0:3, :] = coords_ang.T / bohr_in_angstrom()
        pc[3, :] = q
        return pc

    def write_pot_files(
        self,
        snapshots,
        outdir: str | Path = "pot_frames",
    ):
        """
        Write PE environment snapshots to .pot files.
        Generates one .pot file per snapshot.

        The file contains the @environment,
        @charges, and @polarizabilities sections.

        :param snapshots:
            A list of snapshot dictionaries.
        :param outdir:
            Output directory where the .pot files will be written.
        :param pe_model:
            PE model name
        :raises KeyError:
           If a residue/atom type in the snapshots is missing from the selected PE database.
        """
        if self.pe_model is None:
            raise RuntimeError("PE model is not set. Call set_env_models(pe_model=..., npe_model=...) first.")

        if isinstance(snapshots, dict):
            snapshots = [snapshots]
        
        # Write PE files only on master rank for now.
        # Other ranks wait at the barrier.
        if self.rank != mpi_master():
            self.comm.barrier()
            return
        
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        pe_db = self.pe_model["db"]

        for snap in snapshots:
            frame = int(snap["frame"])

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            pe_elements = np.asarray(snap.get("pe_elements", []), dtype=object)
            pe_resids = np.asarray(snap.get("pe_resids", []), dtype=int)
            pe_resindices = np.asarray(
                snap.get("pe_resindices", snap.get("pe_resids", [])),
                dtype=int,
            )
            pe_resnames = np.asarray(snap.get("pe_resnames", []), dtype=object)
            pe_atom_names = np.asarray(snap.get("pe_atom_names", []), dtype=object)

            npe_coords = np.asarray(snap.get("npe_coords", []), dtype=float)
            npe_elements = np.asarray(snap.get("npe_elements", []), dtype=object)
            npe_resids = np.asarray(snap.get("npe_resids", []), dtype=int)
            npe_resindices = np.asarray(
                snap.get("npe_resindices", snap.get("npe_resids", [])),
                dtype=int,
            )
            npe_resnames = np.asarray(snap.get("npe_resnames", []), dtype=object)
            npe_atom_names = np.asarray(snap.get("npe_atom_names", []), dtype=object)
 
            if pe_coords.size == 0 and npe_coords.size == 0:                
                continue

            if pe_coords.size > 0 and pe_atom_names.size != pe_coords.shape[0]:
                raise ValueError(
                    "pe_atom_names is missing or wrong length in snapshots. "
                    "Required for CP3/SEP when same element has different parameters."
                )

            if pe_coords.size > 0 and pe_elements.size != pe_coords.shape[0]:
                raise ValueError("pe_elements is missing or wrong length in snapshots.")

            if pe_coords.size > 0 and pe_resids.size != pe_coords.shape[0]:
                raise ValueError("pe_resids is missing or wrong length in snapshots.")
            
            if pe_coords.size > 0 and pe_resnames.size != pe_coords.shape[0]:
                raise ValueError("pe_resnames is missing or wrong length in snapshots.")            

            if pe_coords.size > 0 and pe_resindices.size != pe_coords.shape[0]:
                raise ValueError(
                    "pe_resindices is missing or wrong length in snapshots. "
                    "Required for robust residue-based PE validation and writing."
                )
            
            if npe_coords.size > 0 and npe_atom_names.size != npe_coords.shape[0]:
                raise ValueError(
                    "npe_atom_names is missing or wrong length in snapshots. "
                    "Required for NPE charges in mixed PE/NPE pot files."
                )

            if npe_coords.size > 0 and npe_elements.size != npe_coords.shape[0]:
                raise ValueError("npe_elements is missing or wrong length in snapshots.")

            if npe_coords.size > 0 and npe_resids.size != npe_coords.shape[0]:
                raise ValueError("npe_resids is missing or wrong length in snapshots.")

            if npe_coords.size > 0 and npe_resnames.size != npe_coords.shape[0]:
                raise ValueError("npe_resnames is missing or wrong length in snapshots.")

            if npe_coords.size > 0 and npe_resindices.size != npe_coords.shape[0]:
                raise ValueError(
                    "npe_resindices is missing or wrong length in snapshots. "
                    "Required for robust residue-based NPE validation and writing."
                )
            
            self._ensure_no_split_residues_between_pe_and_npe(snap)
            if pe_coords.size > 0:
                self._validate_pe_residue_atom_patterns(
                    pe_atom_names,
                    pe_resindices,
                    pe_resnames,
                )

            pe_resname_set = []
            seen = set()
            for r in pe_resnames.tolist():
                r = str(r)
                if r not in seen:
                    seen.add(r)
                    pe_resname_set.append(r)

            npe_resname_set = []
            seen = set()
            for r in npe_resnames.tolist():
                r = str(r)
                if r not in seen:
                    seen.add(r)
                    npe_resname_set.append(r)

            pot_path = outdir / f"pe_frame_{frame:06d}.pot"
            

            with pot_path.open("w") as fh:
                fh.write("@environment\n")
                fh.write("units: angstrom\n")
                fh.write("xyz:\n")
                for (x, y, z), elem, resn, resid, atom_name in zip(
                        pe_coords, pe_elements, pe_resnames, pe_resids,
                        pe_atom_names):
                    fh.write(
                        f"{str(elem):<2} {x:12.6f} {y:12.6f} {z:12.6f}  "
                        f"{str(resn)}_pe  {int(resid)}  {str(atom_name)}\n")
                for (x, y, z), elem, resn, resid, atom_name in zip(
                        npe_coords, npe_elements, npe_resnames, npe_resids,
                        npe_atom_names):
                    fh.write(
                        f"{str(elem):<2} {x:12.6f} {y:12.6f} {z:12.6f}  "
                        f"{str(resn)}_npe  {int(resid)}  {str(atom_name)}\n")
                fh.write("@end\n\n")

                fh.write("@charges\n")
                for resn in pe_resname_set:
                    if resn not in pe_db:
                        raise KeyError(f"No PE parameters for residue name '{resn}'")
                    pattern_atoms = self._first_residue_atom_pattern(
                        pe_atom_names, pe_resindices, pe_resnames, resn
                    )
                    for atom in pattern_atoms:
                        resolved_atom = self._resolve_atom_name_for_pe_db(atom, pe_db[resn].keys())
                        if resolved_atom is None:
                            raise KeyError(
                                f"No PE params for {resn}/{atom}. Available: {sorted(pe_db[resn].keys())}"
                            )
                        p = pe_db[resn][resolved_atom]
                        fh.write(f"{p['element']:<2} {p['charge']:12.8f}  {resn}_pe\n")
                for resn in npe_resname_set:
                    pattern_atoms = self._first_residue_atom_pattern(
                        npe_atom_names, npe_resindices, npe_resnames, resn
                    )
                    for atom in pattern_atoms:
                        p = self._get_npe_params(atom, resn)
                        fh.write(f"{p['element']:<2} {p['charge']:12.8f}  {resn}_npe\n")
                fh.write("@end\n\n")

                fh.write("@polarizabilities\n")
                for resn in pe_resname_set:
                    pattern_atoms = self._first_residue_atom_pattern(
                        pe_atom_names, pe_resindices, pe_resnames, resn
                    )
                    for atom in pattern_atoms:
                        resolved_atom = self._resolve_atom_name_for_pe_db(atom, pe_db[resn].keys())
                        if resolved_atom is None:
                            raise KeyError(
                                f"No PE params for {resn}/{atom}. Available: {sorted(pe_db[resn].keys())}"
                            )
                        p = pe_db[resn][resolved_atom]
                        pol = p["polar"]
                        fh.write(
                            f"{p['element']:<2} {pol[0]:12.8f} {pol[1]:12.8f} {pol[2]:12.8f} "
                            f"{pol[3]:12.8f} {pol[4]:12.8f} {pol[5]:12.8f}  {resn}_pe\n"
                        )
                fh.write("@end\n")

        # Ensure all ranks proceed only after master finished writing
        self.comm.barrier()    

    @staticmethod
    def _apply_options_to_driver(driver, options: dict, skip_keys: set[str] | None = None):
        """
        Apply dictionary options to a driver using setattr with key validation.
        """
        if skip_keys is None:
            skip_keys = set()
        for key, val in options.items():
            if key in skip_keys:
                continue
            if not hasattr(driver, key):
                raise ValueError(
                    f"Unknown option '{key}' for {type(driver).__name__}."
                )
            setattr(driver, key, val)

    def _build_scf_driver(self, scf_options: dict):
        """
        Build SCF driver from scf_options['scf_type'].
        Allowed values:
            - 'restricted'           : ScfRestrictedDriver
            - 'unrestricted'         : ScfUnrestrictedDriver
            - 'restricted_openshell' : ScfRestrictedOpenDriver
        """
        scf_type = str(scf_options.get("scf_type", "restricted")).lower()
        if scf_type == "restricted":
            return ScfRestrictedDriver(self.comm, self.ostream)
        if scf_type == "unrestricted":
            return ScfUnrestrictedDriver(self.comm, self.ostream)
        if scf_type == "restricted_openshell":
            return ScfRestrictedOpenDriver(self.comm, self.ostream)
        raise ValueError(
            "Invalid scf_type in scf_options. "
            "Expected one of: 'restricted', 'unrestricted', 'restricted_openshell'."
        )

    def _build_property_driver(self, property_options: dict):
        """
        Build response driver from property_options using rules:
        - property='absorption' or 'ecd' + nstates      -> LinearResponseEigenSolver
        - property='absorption' or 'ecd' + frequencies  -> ComplexResponse
        """
        if not property_options:
            return None
 
        if "property" not in property_options:
            raise ValueError("property_options must contain key 'property'.")
 
        prop = str(property_options["property"]).lower().strip()
        has_nstates = "nstates" in property_options
        has_freqs = "frequencies" in property_options
 
        if prop not in {"absorption", "ecd"}:
            raise ValueError(
                "Invalid property in property_options. "
                "Expected 'absorption' or 'ecd'."
            )
 
        if has_nstates and has_freqs:
            raise ValueError(
                "property_options cannot contain both 'nstates' and 'frequencies'."
            )
 
        if has_freqs:
            drv = ComplexResponseSolver(self.comm, self.ostream)
            drv.set_cpp_property("absorption" if prop == "absorption" else "ecd")
            return drv
 
        if has_nstates:
            return LinearResponseEigenSolver(self.comm, self.ostream)
 
        raise ValueError(
            "property_options must define either 'nstates' (TD-DFT) "
            "or 'frequencies' (CPP)."
        )

    def compute(
        self,
        snapshots,
        basis_set: str,
        scf_options: dict | None = None,
        property_options: dict | None = None,
        potdir: str | Path = "pot_frames",
        write_pe_potfiles: bool = True,
        qm_charge: int | None = None,
        qm_multiplicity: int | None = None,
    ):
        """
        Drives the computation over the ensemble of snapshots.

        For each snapshot, an SCF calculation is performed for the QM subsystem, with
        optional PE/NPE environment terms. Optional response calculations are
        controlled via property_options.

        :param snapshots:
            A list of snapshot dictionaries (or a single snapshot dict).
        :param basis_set: (str)
            Basis set.
        :param scf_options: (dict)
            SCF/method options merged in one dictionary.
            Includes optional key 'scf_type' in {'restricted', 'unrestricted', 'restricted_openshell'}
            to select the SCF driver.
        :param property_options: (dict)
            Property options dictionary with keys:
            'property' ('absorption' or 'ecd') + ('nstates' or 'frequencies').
        :param potdir : (str or Path)
            Directory to store/read PE potfiles.
        :param write_pe_potfiles: (bool)
            If True, PE potfiles are (re)generated before the loop when needed.
        :param qm_charge:
            Optional override for the QM-region charge. If None, the value stored
            in each snapshot is used, defaulting to 0 when absent.
        :param qm_multiplicity:
            Optional override for the QM-region multiplicity. If None, the value stored
            in each snapshot is used, defaulting to 1 when absent.
        
        :return:
            Dictionary with keys:
            - scf_all: list of (frame, scf_results)
            - rsp_all: list of (frame, rsp_results), only present if response is run.

        :raises RuntimeError:
            If required PE/NPE models have not been selected.
        :raises ValueError:
            If snapshot fields required to build NPE point charges are missing.
        """
   
        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        if scf_options is None:
            scf_options = {}
        else:
            scf_options = dict(scf_options)

        if property_options is None:
            property_options = {}

        if qm_charge is not None:
            qm_charge = int(qm_charge)
        if qm_multiplicity is not None:
            qm_multiplicity = int(qm_multiplicity)
            if qm_multiplicity <= 0:
                raise ValueError("QM multiplicity must be a positive integer.")

        potdir = Path(potdir)

        # Detect whether PE / NPE models are needed
        has_any_pe = any(np.asarray(s.get("pe_coords", [])).size > 0 for s in snapshots)
        has_any_npe = any(np.asarray(s.get("npe_coords", [])).size > 0 for s in snapshots)

        # Only require environment models when corresponding environment atoms exist
        if has_any_pe and self.pe_model is None:
            raise RuntimeError(
                "Snapshots contain PE atoms but pe_model is not set. "
                "Call set_env_models(pe_model=...)."
            )

        if has_any_npe and self.npe_model is None:
            raise RuntimeError(
                "Snapshots contain NPE atoms but npe_model is not set. "
                "Call set_env_models(npe_model=...)."
            )

        for snap in snapshots:
            self._ensure_no_split_residues_between_pe_and_npe(snap)

        if write_pe_potfiles and has_any_pe:
            self.write_pot_files(snapshots, outdir=potdir)

        do_rsp = bool(property_options)

        # Guard against open-shell response requests for now:
        if do_rsp:
            scf_type = str(scf_options.get("scf_type", "restricted")).lower()
            if scf_type != "restricted":
                raise NotImplementedError(
                    "Open-shell response is not yet supported in EnsembleDriver.compute. "
                    "Use scf_type='restricted' for response calculations."
                )

        scf_all = []
        rsp_all = [] if do_rsp else None

        for snap in snapshots:
            scf_driver = self._build_scf_driver(scf_options)
            self._apply_options_to_driver(scf_driver, scf_options, skip_keys={"scf_type"})

            frame = int(snap["frame"])

            labels = [str(x) for x in snap["qm_elements"]]
            coords = np.asarray(snap["qm_coords"], dtype=float)

            snap_qm_charge = (
                qm_charge
                if qm_charge is not None 
                else int(snap.get("qm_charge", 0))
            )
            snap_qm_multiplicity = (
                qm_multiplicity
                if qm_multiplicity is not None 
                else int(snap.get("qm_multiplicity", 1))
            )

            molecule = Molecule(labels, coords)
            molecule.set_charge(snap_qm_charge)
            molecule.set_multiplicity(snap_qm_multiplicity)

            if not molecule.check_multiplicity():
                raise ValueError(
                    f"Incompatible QM charge ({snap_qm_charge}) and multiplicity"
                    f"({snap_qm_multiplicity}) for frame {frame}."
                )
            basis = MolecularBasis.read(molecule, basis_set)

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            npe_coords = np.asarray(snap.get("npe_coords", []), dtype=float)

            has_pe = pe_coords.size > 0
            has_npe = npe_coords.size > 0

            # Reset embedding inputs every frame
            scf_driver.potfile = None
            scf_driver.point_charges = None
            
            if has_pe:
                scf_driver.potfile = str(potdir / f"pe_frame_{frame:06d}.pot")

            if has_npe and not has_pe:
                npe_atom_names = snap.get("npe_atom_names", [])
                npe_resnames = snap.get("npe_resnames", [])
                if len(npe_atom_names) == 0:
                    raise ValueError("npe_atom_names missing in snapshots (required to read NPE charges from CSV)."
                    )
                scf_driver.point_charges = self._build_point_charges(
                    npe_coords, npe_atom_names, npe_resnames
                )

            scf_results = scf_driver.compute(molecule, basis)
            scf_all.append((frame, scf_results))

            if do_rsp:
                rsp_driver = self._build_property_driver(property_options)
                # property is routing metadata; do not setattr it on the driver
                self._apply_options_to_driver(
                    rsp_driver,
                    property_options,
                    skip_keys={"property"},
                )
                rsp_results = rsp_driver.compute(molecule, basis, scf_results)
                rsp_all.append((frame, rsp_results))

        results = {"scf_all": scf_all}
        if do_rsp:
            results["rsp_all"] = rsp_all
        return results

    def plot_uv_vis_spectra(
        self,
        results: dict,
        energy_min_ev=None,
        energy_max_ev=None,
        show_individual: bool = False,
        show_sticks: bool = True,
        show_std: bool = False,
        title: str = "Absorption Spectrum (Averaged)",
        ax=None,
        xlim_nm=None,
        save_averaged_spectra: bool = True,
        averaged_spectra_filename: str | Path = "averaged_spectra.csv",
    ):
        """
        Convenience wrapper to plot averaged UV/Vis spectra for an ensemble.

        This method expects the `results` dictionary returned by :meth:`compute`.
        It extracts `results["rsp_all"]` and passes it to
        :class:`SpectrumAverager`.

        :param results:
            Results dictionary returned by :meth:`compute`.
        :param energy_min_ev:
            Minimum photon energy in eV for the common grid. If None, it is set to
            (min excitation energy - padding).
        :param energy_max_ev:
            Maximum photon energy in eV for the common grid. If None, it is set to
            (max excitation energy + padding).
        :param show_individual:
            If True, plot individual broadened spectra for each snapshot.
        :param show_sticks:
            If True, plots oscillator strength sticks at transition wavelengths.
        :param show_std:
            If True, shows a shaded area corresponding to +/- one standard deviation.
        :param title:
            Title of the plot.
        :param ax:
            Matplotlib Axes object to plot on. If None, a new figure and axes are
            created.
        :param xlim_nm:
            Tuple (xmin, xmax) to set x-axis limits in nm. If None, automatic limits are used.
        :param save_averaged_spectra:
            If True, saves the averaged spectra data to a CSV.
        :param averaged_spectra_filename:
            Filename for saving the averaged spectra CSV if `save_averaged_spectra` is True.
            Default is "averaged_spectra.csv" in current working directory.
        :return:
            The Matplotlib Axes object containing the plot.
        :raises KeyError:
            If `rsp_all` is not present in the results dictionary.
        """
        if not isinstance(results, dict):
            raise TypeError(
                "results must be a dictionary as returned by EnsembleDriver.compute()."
            )

        if "rsp_all" not in results or results.get("rsp_all", None) is None:
            raise KeyError(
                "No 'rsp_all' found in results. Make sure you ran a response calculation "
                "(e.g., pass property_optios with 'property' and either "
                "'nstates' or 'frequencies' to EnsembleDriver.compute())"
            )

        rsp_all = results["rsp_all"]

        spec_avg = SpectrumAverager(comm=self.comm, ostream=self.ostream)
        ax = spec_avg.plot_uv_vis_spectra(
            rsp_all,
            energy_min_ev=energy_min_ev,
            energy_max_ev=energy_max_ev,
            show_individual=show_individual,
            show_sticks=show_sticks,
            show_std=show_std,
            title=title,
            ax=ax,
            xlim_nm=xlim_nm,
            save_averaged_spectra=save_averaged_spectra,
            averaged_spectra_filename=averaged_spectra_filename,
        )
        return ax

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
import numpy as np
import csv
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .lreigensolver import LinearResponseEigenSolver
from .sanitychecks import ensemble_driver_scf_sanity_check
from .sanitychecks import ensemble_driver_rsp_sanity_check
from .errorhandler import assert_msg_critical

from .veloxchemlib import (mpi_master, bohr_in_angstrom)

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

        # settings for SCF and choice of rsp driver
        self._is_restricted = True
        self._rsp_property = None

        # Other possible settings
        self.xcfun = None
        self.potfile = None
        self.embedding = None
        self.eri_thresh = None
        self.ri_coulomb = None
        self.grid_level = None

        # added scf_dict, method_dict and rsp_dict as
        # instance variables. Updated in update_settings
        # and used in compute
        self.scf_dict = None
        self.method_dict = None
        self.rsp_dict = None

        self.excited_states = False

        # response settings are also added here in a similar way
        self.nstates = None
        self.nto = None
        self.core_excitation = None
        self.num_core_orbitals = None

        db_dir = Path(__file__).resolve().parent / "database" / "environment_parameters"

        sep_parameters_file = db_dir / "pe_sep.csv"
        cp3_parameters_file = db_dir / "pe_cp3.csv"
        tip3p_parameters_file = db_dir / "npe_tip3p.csv"
        ff19sb_parameters_file = db_dir / "npe_ff19sb.csv"


        if self.rank == mpi_master():
            self._sep_db = self._load_pe_db(sep_parameters_file)
            self._cp3_db = self._load_pe_db(cp3_parameters_file)
            self._tip3p_db = self._load_npe_db(tip3p_parameters_file)
            self._ff19sb_db = self._load_npe_db(ff19sb_parameters_file)
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
    
    # A good routine to have to enable input-output runs later one.
    def update_settings(self, scf_dict=None, method_dict=None, rsp_dict=None):
        """
        Updates settings in the ensemble driver.

        :param scf_dict:
            The dictionary of SCF settings.
        :param method_dict:
            The dictionary of method settings.
        :param rsp_dict:
            The dictionary of response settings.
        """
        if scf_dict is None:
            scf_dict = {}
        if method_dict is None:
            method_dict = {}
        if rsp_dict is None:
            rsp_dict = {}

        self.scf_dict = dict(scf_dict)
        self.method_dict = dict(method_dict)
        self.rsp_dict = dict(rsp_dict)

        # Auto-enable response if response settings are provided
        if self.rsp_dict:
            self.excited_states = True

    @staticmethod
    def _parse_six_floats(field: str) -> list[float]:
        """
        Parse a six-component polarizability field from a CSV entry.

        :param field:
            Unoyt field containing six floating-point numbers.

        :return:
            List of six floats.
        """
        parts = [p for p in str(field).replace(",", " ").split() if p]
        if len(parts) != 6:
            raise ValueError(f"Expected 6 floats in P11, got {len(parts)} from: {field!r}")
        return [float(x) for x in parts]
    
    @staticmethod
    def _load_pe_db(csv_path: Path) -> dict:
        """
        Loads PE-like table:
            molecule,res_name,atom_name,element,M0,P11

        Returns:
            db[res_name][atom_name] = {"element": str, "charge": float, "polar": [6 floats]}
        """
        if not csv_path.is_file():
            raise FileNotFoundError(f"PE parameter file not found: {csv_path}")

        db: dict[str, dict[str, dict]] = {}

        with csv_path.open("r", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                resn = str(row["res_name"]).strip()
                atom = str(row["atom_name"]).strip()
                elem = str(row["element"]).strip()
                if not resn or not atom or not elem:
                    continue

                q = float(row["M0"])
                pol6 = EnsembleDriver._parse_six_floats(row["P11"])

                rdb = db.setdefault(resn, {})
                if atom in rdb:
                    old = rdb[atom]
                    if abs(old["charge"] - q) > 1e-12 or any(
                        abs(a - b) > 1e-12 for a, b in zip(old["polar"], pol6)
                    ):
                        raise ValueError(
                            f"Inconsistent duplicate entries for {resn}/{atom} in {csv_path}"
                        )
                rdb[atom] = {"element": elem, "charge": q, "polar": pol6}

        return db

    @staticmethod
    def _load_npe_db(csv_path: Path) -> dict:
        """
        Loads NPE-like table:
            molecule,res_name,atom_name,element,M0

        Returns:
            db[res_name][atom_name] = {"element": str, "charge": float}
        """
        if not csv_path.is_file():
            raise FileNotFoundError(f"NPE parameter file not found: {csv_path}")

        db: dict[str, dict[str, dict]] = {}

        with csv_path.open("r", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                resn = str(row["res_name"]).strip()
                atom = str(row["atom_name"]).strip()
                elem = str(row["element"]).strip()
                if not resn or not atom or not elem:
                    continue

                q = float(row["M0"])

                rdb = db.setdefault(resn, {})
                if atom in rdb:
                    old = rdb[atom]
                    if abs(old["charge"] - q) > 1e-12 or old["element"] != elem:
                        raise ValueError(
                            f"Inconsistent duplicate entries for {resn}/{atom} in {csv_path}"
                        )
                rdb[atom] = {"element": elem, "charge": q}
        return db

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
                                         water, and ions.

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
    def _first_residue_atom_pattern(atom_names, resids, resnames, target_resname: str) -> list[str]:
        """
        Return atom_name pattern for one residue instance of target_resname,
        preserving the order in the arrays.
        """
        atom_names = np.asarray(atom_names, dtype=object)
        resids = np.asarray(resids, dtype=int)
        resnames = np.asarray(resnames, dtype=object)

        mask = (resnames == target_resname)
        if not np.any(mask):
            return []

        first_resid = int(resids[mask][0])
        idx = np.where(mask & (resids == first_resid))[0]
        return [str(atom_names[i]) for i in idx]

    def _build_point_charges(self, coords_ang, atom_names, resnames) -> np.ndarray | None:
        """
        Build point charges array expected by SCF driver: shape (6, N), coords in bohr.
        Charges are taken from NPE database (db) by (resname, atom_name).
        """
        coords_ang = np.asarray(coords_ang, dtype=float)
        if coords_ang.size == 0:
            return None
        
        if self.npe_model is None:
            raise RuntimeError("Snapshot contains NPE atoms but npe_model is not set.")

        atom_names = np.asarray(atom_names, dtype=object)
        resnames = np.asarray(resnames, dtype=object)
        db = self.npe_model["db"]

        q = np.empty(atom_names.size, dtype=float)
        for i, (atom, resn) in enumerate(zip(atom_names, resnames)):
            resn = str(resn)
            atom = str(atom)
            if resn not in db:
                raise KeyError(f"No NPE parameters for residue name '{resn}'")
            if atom not in db[resn]:
                raise KeyError(f"No NPE charge for {resn}/{atom}. "
                               f"Available: {sorted(db[resn].keys())}")
            q[i] = db[resn][atom]["charge"]

        pc = np.zeros((6, coords_ang.shape[0]), dtype=float)
        pc[0:3, :] = coords_ang.T / bohr_in_angstrom()
        pc[3, :] = q
        return pc

    def write_pot_files(self, snapshots, outdir: str | Path):
        """
        Write PE environment snapshots to .pot files.

        Generates one .pot file per snapshot, named ``pe_frame_XXXXXX.pot`` 
        where XXXXXX is the snapshot ``frame`` index.

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

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        pe_db = self.pe_model["db"]

        for snap in snapshots:
            frame = int(snap["frame"])

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            pe_elements = np.asarray(snap.get("pe_elements", []), dtype=object)
            pe_resids = np.asarray(snap.get("pe_resids", []), dtype=int)
            pe_resnames = np.asarray(snap.get("pe_resnames", []), dtype=object)
            pe_atom_names = np.asarray(snap.get("pe_atom_names", []), dtype=object)

            if pe_coords.size == 0:
                continue

            if pe_atom_names.size != pe_coords.shape[0]:
                raise ValueError(
                    "pe_atom_names is missing or wrong length in snapshots. "
                    "Required for CP3/SEP when same element has different parameters."
                )

            # Unique residue names
            resname_set = []
            seen = set()
            for r in pe_resnames.tolist():
                r = str(r)
                if r not in seen:
                    seen.add(r)
                    resname_set.append(r)

            pot_path = outdir / f"pe_frame_{frame:06d}.pot"
            

            with pot_path.open("w") as fh:
                fh.write("@environment\n")
                fh.write("units: angstrom\n")
                fh.write("xyz:\n")
                for (x, y, z), elem, resn, resid in zip(pe_coords, pe_elements, pe_resnames, pe_resids):
                    # print (resid)
                    fh.write(f"{str(elem):<2} {x:12.6f} {y:12.6f} {z:12.6f}  {str(resn):>3}  {int(resid)}\n")
                fh.write("@end\n\n")

                fh.write("@charges\n")
                for resn in resname_set:
                    if resn not in pe_db:
                        raise KeyError(f"No PE parameters for residue name '{resn}'")
                    pattern_atoms = self._first_residue_atom_pattern(pe_atom_names, pe_resids, pe_resnames, resn)
                    for atom in pattern_atoms:
                        if atom not in pe_db[resn]:
                            raise KeyError(f"No PE params for {resn}/{atom}. Available: {sorted(pe_db[resn].keys())}")
                        p = pe_db[resn][atom]
                        fh.write(f"{p['element']:<2} {p['charge']:12.8f}  {resn}\n")
                fh.write("@end\n\n")

                fh.write("@polarizabilities\n")
                for resn in resname_set:
                    pattern_atoms = self._first_residue_atom_pattern(pe_atom_names, pe_resids, pe_resnames, resn)
                    for atom in pattern_atoms:
                        p = pe_db[resn][atom]
                        pol = p["polar"]
                        fh.write(
                            f"{p['element']:<2} {pol[0]:12.8f} {pol[1]:12.8f} {pol[2]:12.8f} "
                            f"{pol[3]:12.8f} {pol[4]:12.8f} {pol[5]:12.8f}  {resn}\n"
                        )
                fh.write("@end\n")

    def compute(
        self,
        snapshots,
        basis_label: str,
        potdir: str | Path = "pot_frames",
        write_pe_potfiles: bool = True,
    ):
        """
        Drives the computation over the ensemble of snapshots.

        For each snapshot, an SCF calculation is performed for the QM subsystem, with
        optional PE/NPE environment terms. Linear-response calculations are optional
        and are controlled by :attr:`excited_states` and the response options.

        :param snapshots:
            A list of snapshot dictionaries (or a single snapshot dict).
        :param basis_label: (str)
            Basis set label.
        :param potdir : (str or Path)
            Directory to store/read PE potfiles.
        :param write_pe_potfiles: (bool)
            If True, PE potfiles are (re)generated before the loop when needed.
        
        :return:
            Dictionary with keys:
            - scf_all: list of (frame, scf_results)
            - rsp_all: list of (frame, rsp_results), only present if response is run.

        :raises RuntimeError:
            If required PE/NPE models have not been selected.
        :raises ValueError:
            If snapshot fields required to build NPE point charges are missing.
        """
        if self.pe_model is None and self.npe_model is None:
            raise RuntimeError(
            "Models not set. Call set_env_models(pe_model=..., npe_model=...) first."
            )
    
        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        potdir = Path(potdir)

        # Detect whether we actually need PE / NPE
        has_any_pe = any(np.asarray(s.get("pe_coords", [])).size > 0 for s in snapshots)
        has_any_npe = any(np.asarray(s.get("npe_coords", [])).size > 0 for s in snapshots)

        if has_any_pe and self.pe_model is None:
            raise RuntimeError("Snapshots contain PE atoms but pe_model is not set. Call set_env_models(pe_model=...).")

        if has_any_npe and self.npe_model is None:
            raise RuntimeError("Snapshots contain NPE atoms but npe_model is not set. Call set_env_models(npe_model=...).")

        if write_pe_potfiles and has_any_pe:
            self.write_pot_files(snapshots, outdir=potdir)

        # Here starts the suggestion of how to handle 
        # SCF under the hood
        if self._is_restricted:
            scf_driver = ScfRestrictedDriver()
        else:
            scf_driver = ScfUnrestrictedDriver()

        rsp_driver = LinearResponseEigenSolver()
        
        # update settings for scf if necessary
        if self.scf_dict is not None or self.method_dict is not None:
            scf_driver.update_settings(self.scf_dict or {}, self.method_dict or {})

        do_rsp = bool(self.excited_states)
        if not do_rsp:
            rsp_opts_set = (
                (self.rsp_dict is not None and len(self.rsp_dict) > 0)
            or (self.nstates is not None)
            or (self.nto is not None)
            or (self.core_excitation is not None)
            or (self.num_core_orbitals is not None)
            )
            if rsp_opts_set:
                do_rsp = True
                self.excited_states = True
        
        rsp_driver = None
        if do_rsp:
            rsp_driver = LinearResponseEigenSolver()

            # update settings for rsp if necessary
            if self.rsp_dict is not None or self.method_dict is not None:
                rsp_driver.update_settings(self.rsp_dict or {}, self.method_dict or {})

        # sanity check -- update settings in scf and rsp according to
        # which variables the user set.
        ensemble_driver_scf_sanity_check(scf_driver, self)
        if do_rsp:
            ensemble_driver_rsp_sanity_check(rsp_driver, self)
        
        scf_all = []
        rsp_all = [] if do_rsp else None

        for snap in snapshots:
            frame = int(snap["frame"])

            labels = [str(x) for x in snap["qm_elements"]]
            coords = np.asarray(snap["qm_coords"], dtype=float)
            molecule = Molecule(labels, coords)
            basis = MolecularBasis.read(molecule, basis_label)

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            npe_coords = np.asarray(snap.get("npe_coords", []), dtype=float)

            has_pe = pe_coords.size > 0
            has_npe = npe_coords.size > 0

            # Reset embedding inputs every frame
            scf_driver.potfile = None
            scf_driver.point_charges = None
            
            if has_pe:
                scf_driver.potfile = str(potdir / f"pe_frame_{frame:06d}.pot")

            if has_npe:
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
                rsp_results = rsp_driver.compute(molecule, basis, scf_results)
                rsp_all.append((frame, rsp_results))

        results = {"scf_all": scf_all}
        if do_rsp:
            results["rsp_all"] = rsp_all
        return results

"""
redox_fetcher.py
================
Load and analyse results produced by RedoxCalculator from a folder of HDF5
files, and export a tidy summary CSV.

Design
------
- The HDF5 folder is scanned once on first use and the index is cached.
  Calling ``fetch()`` 360 times costs one disk scan, not 360.
- State matching uses graph isomorphism (not GED) for correctness and speed.
- Physical constants are imported from redox_calculator so there is a single
  source of truth.
- Output is wide-format: one row per molecule, one column per property.

Usage
-----
::

    from redox_fetcher import RedoxFetcher

    fetcher = RedoxFetcher("redox_results", pH=7.0)

    # Single molecule
    result = fetcher.fetch("c1ccccc1O", charge=0, mult=1)
    print(result)

    # Batch from CSV -> summary CSV
    fetcher.batch_summary("redox_molecules.csv", "redox_results_summary.csv")
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from typing import Optional

import h5py
import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms import isomorphism
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Import physical constants from the calculator so there is one source of truth.
# If redox_calculator is not on the path, fall back to inline values.
try:
    from redox_calculator import (
        G_PROTON_AQ, G_ELECTRON, CONV_EV, AU_TO_KCAL,
        E_REF_SHE, RT_LN10_KCAL, NERNST_FACTOR,
    )
except ImportError:
    G_PROTON_AQ  = -0.4390
    G_ELECTRON   = -0.0015
    CONV_EV      = 27.2114
    AU_TO_KCAL   = 627.509
    E_REF_SHE    =  4.28
    RT_LN10_KCAL =  1.364
    NERNST_FACTOR =  0.05916

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Index entry — one per HDF5 file
# ---------------------------------------------------------------------------

@dataclass
class _H5Entry:
    """Lightweight record for one HDF5 result file."""
    filename:          str
    g_total:           float
    g_corr:            float
    charge:            int
    mult:              int
    h_count:           int
    heavy_form:        str
    graph:             Optional[nx.Graph]
    skeleton:          str
    mm_corr_available: bool = True
    solvation_data:    dict = field(default_factory=dict)

    @property
    def full_formula(self) -> str:
        return f"{self.heavy_form}H{self.h_count}" if self.h_count else self.heavy_form


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _heavy_formula(formula_str: str) -> str:
    """
    Extract sorted heavy-atom counts from a molecular formula string.
    'C14H8N2O6' -> 'C14N2O6'
    """
    counts: dict[str, int] = {}
    for elem, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula_str):
        if elem == "H":
            continue
        counts[elem] = counts.get(elem, 0) + (int(count) if count else 1)
    return "".join(f"{k}{counts[k]}" for k in sorted(counts))


def _neutral_formula(mol) -> str:
    """
    Molecular formula from an RDKit Mol, with charge symbols stripped.
    'C6H5O2+' -> 'C6H5O2'
    """
    if mol is None:
        return ""
    return re.sub(r'[+\-]\d*$', '', CalcMolFormula(mol))


def _mol_to_graph(mol) -> Optional[nx.Graph]:
    """Heavy-atom connectivity graph from an RDKit Mol."""
    if mol is None:
        return None
    mol = Chem.RemoveHs(mol)
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), element=atom.GetSymbol())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G


def _graphs_isomorphic(g1: nx.Graph, g2: nx.Graph) -> bool:
    """True if g1 and g2 are graph-isomorphic (same elements, same topology)."""
    if len(g1) != len(g2):
        return False
    nm = isomorphism.categorical_node_match("element", None)
    return isomorphism.GraphMatcher(g1, g2, node_match=nm).is_isomorphic()


def _h_count_from_mol(mol) -> int:
    """Total implicit H count on the heavy-atom skeleton."""
    if mol is None:
        return 0
    return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())


def _read_solvation_data(h5file) -> dict[str, float]:
    """
    Read any named solvation values stored as 'solvation_data_<name>' keys.
    Returns a dict {name: value_au}.
    """
    result = {}
    for key in h5file.keys():
        if key.startswith("solvation_data_"):
            name = key[len("solvation_data_"):]
            try:
                result[name] = float(h5file[key][()])
            except Exception:
                pass
    return result


# ---------------------------------------------------------------------------
# Thermochemistry  (mirrors RedoxCalculator methods exactly)
# ---------------------------------------------------------------------------

def _reduction_potential(g_ox: float, g_red: float,
                          n: int = 1, m: int = 0, pH: float = 7.0) -> float:
    """E (V vs SHE) for  Ox + n e- + m H+ -> Red."""
    dg_au = g_red - (g_ox + m * G_PROTON_AQ + n * G_ELECTRON)
    e_abs = -dg_au * CONV_EV / n
    return (e_abs - E_REF_SHE) - (NERNST_FACTOR / n) * m * pH


def _pka(g_acid: float, g_base: float) -> float:
    """pKa for  HA -> A- + H+.  g_acid = G(HA), g_base = G(A-)."""
    dg_kcal = ((g_base + G_PROTON_AQ) - g_acid) * AU_TO_KCAL
    return dg_kcal / RT_LN10_KCAL


# ---------------------------------------------------------------------------
# Main fetcher class
# ---------------------------------------------------------------------------

class RedoxFetcher:
    """
    Load and query results from a RedoxCalculator output folder.

    The HDF5 folder is indexed once on first use (or on ``reload()``).
    Subsequent ``fetch()`` calls operate entirely in memory.

    Parameters
    ----------
    output_folder : str
        Path to the folder containing HDF5 result files.
    pH : float
        Solution pH for potential and pKa calculations.  Default 7.0.

    Examples
    --------
    ::

        fetcher = RedoxFetcher("redox_results", pH=7.0)

        # Single lookup
        result = fetcher.fetch("c1ccccc1O", charge=0, mult=1)

        # Batch export
        fetcher.batch_summary("redox_molecules.csv",
                              "redox_results_summary.csv")
    """

    def __init__(self, output_folder: str = "redox_results", pH: float = 7.0):
        self.output_folder = output_folder
        self.pH            = pH
        self._index: list[_H5Entry] = []   # populated lazily on first use

    # ------------------------------------------------------------------
    # Index management
    # ------------------------------------------------------------------

    def _ensure_index(self) -> None:
        """Build the index if it has not been built yet."""
        if not self._index:
            self.reload()

    def reload(self) -> None:
        """
        Re-scan ``output_folder`` and rebuild the in-memory index.

        Call this if new HDF5 files have been added since the last scan.
        """
        self._index = []
        if not os.path.isdir(self.output_folder):
            log.warning("Output folder not found: %s", self.output_folder)
            return

        n_loaded = n_skipped = 0
        for fname in os.listdir(self.output_folder):
            if not fname.endswith(".h5"):
                continue
            path = os.path.join(self.output_folder, fname)
            try:
                with h5py.File(path, "r") as f:
                    g_total   = float(f["g_total"][()])
                    g_corr    = float(f["g_corr"][()]) if "g_corr" in f else 0.0
                    charge    = int(f["charge"][()])
                    mult      = int(f["multiplicity"][()])
                    skeleton  = (
                        f["skeleton_smiles"][()].decode("utf-8")
                        if "skeleton_smiles" in f else "Unknown"
                    )
                    # mm_corr_available was added in a later version;
                    # older files default to True (assume correction was applied).
                    mm_ok = bool(f["mm_corr_available"][()]) if "mm_corr_available" in f else True
                    solv_data = _read_solvation_data(f)

                mol    = Chem.MolFromSmiles(skeleton) if skeleton != "Unknown" else None
                graph  = _mol_to_graph(mol)
                h_cnt  = _h_count_from_mol(mol)
                fname_stem = fname.split("_q")[0]
                heavy_form = _heavy_formula(fname_stem)

                self._index.append(_H5Entry(
                    filename          = fname,
                    g_total           = g_total,
                    g_corr            = g_corr,
                    charge            = charge,
                    mult              = mult,
                    h_count           = h_cnt,
                    heavy_form        = heavy_form,
                    graph             = graph,
                    skeleton          = skeleton,
                    mm_corr_available = mm_ok,
                    solvation_data    = solv_data,
                ))
                n_loaded += 1

            except Exception as exc:
                log.warning("Skipping %s: %s", fname, exc)
                n_skipped += 1

        log.info(
            "Index built: %d entries loaded, %d skipped from %s.",
            n_loaded, n_skipped, self.output_folder,
        )

    # ------------------------------------------------------------------
    # Core fetch
    # ------------------------------------------------------------------

    def fetch(
        self,
        query_smiles: str,
        charge: int,
        mult: int,
        solvation_name: Optional[str] = None,
    ) -> Optional[dict]:
        """
        Find the anchor state for ``query_smiles`` and assemble all related
        charge/protonation states into a thermodynamic result dict.

        Parameters
        ----------
        query_smiles : str
            SMILES of the molecule to look up.
        charge : int
            Formal charge of the anchor (reference) state.
        mult : int
            Spin multiplicity of the anchor state.
        solvation_name : str or None
            If given, use the named solvation free energy (stored in
            ``solvation_data``) instead of ``g_total`` when computing
            potentials and pKas.  The named value must have been registered
            via ``RedoxCalculator.register_solvation()`` before the
            calculation was run.

        Returns
        -------
        dict or None
            Dictionary with keys:

            ``"anchor_file"``  — filename of the matched anchor HDF5.
            ``"skeleton"``     — skeleton SMILES from the anchor file.
            ``"states"``       — dict mapping state label to ``_H5Entry``.
            ``"E_ox"``         — oxidation potential (V vs SHE) or None.
            ``"E_red"``        — reduction potential (V vs SHE) or None.
            ``"E_pcet_red"``   — PCET reduction potential or None.
            ``"E_pcet_ox"``    — PCET oxidation potential or None.
            ``"pKa_MH_plus"``  — pKa of MH+ -> M + H+ or None.
            ``"pKa_M"``        — pKa of M -> M_dep- + H+ or None.

            Returns None if no matching anchor is found.

        Notes
        -----
        If multiple anchor files match the query (can happen when several
        site-specific geometries were computed) the one with the lowest
        ``g_total`` is used.
        """
        self._ensure_index()

        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            log.error("Could not parse SMILES: %s", query_smiles)
            return None

        query_graph   = _mol_to_graph(query_mol)
        query_formula = _neutral_formula(query_mol)
        query_heavy   = _heavy_formula(query_formula)
        query_h       = _h_count_from_mol(query_mol)

        # ---- find anchor(s) -----------------------------------------
        # An anchor must match: heavy formula, H count, charge, mult,
        # and graph isomorphism.
        anchors = [
            e for e in self._index
            if (e.charge      == charge
                and e.mult    == mult
                and e.heavy_form == query_heavy
                and e.h_count == query_h
                and e.graph is not None
                and query_graph is not None
                and _graphs_isomorphic(e.graph, query_graph))
        ]

        if not anchors:
            log.debug("No anchor found for %s (q=%d, m=%d).", query_smiles, charge, mult)
            return None

        # If multiple anchors match take the lowest energy one
        if len(anchors) > 1:
            log.warning(
                "%d anchors matched %s — using lowest G_total.",
                len(anchors), query_smiles,
            )
        anchor = min(anchors, key=lambda e: e.g_total)

        # ---- find related states ------------------------------------
        # A related file must share the same heavy-atom skeleton (isomorphic
        # graph) and differ from the anchor only by ±1 H and/or ±1 charge.
        states: dict[str, _H5Entry] = {"M": anchor}

        for e in self._index:
            if e.filename == anchor.filename:
                continue
            if e.heavy_form != anchor.heavy_form:
                continue
            if e.graph is None or not _graphs_isomorphic(e.graph, anchor.graph):
                continue

            dh = e.h_count - anchor.h_count
            dq = e.charge  - anchor.charge

            key: Optional[str] = None
            if   dh ==  0 and dq == +1: key = "M_ox"
            elif dh ==  0 and dq == -1: key = "M_red"
            elif dh == +1 and dq == +1: key = "MH_plus"
            elif dh == +1 and dq ==  0: key = "MH_rad"
            elif dh == -1 and dq == -1: key = "M_deprot"
            elif dh == -1 and dq ==  0: key = "M_deprot_rad"

            if key:
                # Keep the lowest-energy candidate for each state
                if key not in states or e.g_total < states[key].g_total:
                    states[key] = e

        # ---- helper: get the effective G for a state ----------------
        def _g(state_key: str) -> Optional[float]:
            e = states.get(state_key)
            if e is None:
                return None
            if solvation_name:
                if solvation_name not in e.solvation_data:
                    log.warning(
                        "Solvation '%s' not found for state '%s' (%s). "
                        "Falling back to default g_total.",
                        solvation_name, state_key, e.filename,
                    )
                    return e.g_total
                # Reconstruct G using: scf_vacuum + dG_solv + g_corr
                # We only have g_total and solvation_data here; to recover
                # scf_vacuum we would need it stored separately. Since we
                # store g_total = scf_solvent + g_corr, the cleanest approach
                # for the fetcher is to store the named G directly.
                # Convention: solvation_data[name] IS the replacement g_total.
                return e.solvation_data[solvation_name]
            return e.g_total

        # ---- helper: normalised G for a pair -----------------------
        # Mirrors the thermal normalisation in RedoxCalculator.run():
        # if one state has a MM thermal correction and the other doesn't,
        # strip the correction from the one that has it so both are at
        # the pure E_SP(solvent) level for a fair comparison.
        def _normalise_pair(key_a: str, key_b: str) -> tuple:
            """
            Return (g_a, g_b, sp_only_states) normalised to the same level.
            sp_only_states is a list of state keys where g_corr was stripped.
            Returns (None, None, []) if either state is missing.
            """
            e_a = states.get(key_a)
            e_b = states.get(key_b)
            if e_a is None or e_b is None:
                return None, None, []

            g_a = _g(key_a)
            g_b = _g(key_b)
            stripped = []

            if e_a.mm_corr_available and not e_b.mm_corr_available:
                g_a -= e_a.g_corr
                stripped.append(key_a)
            elif e_b.mm_corr_available and not e_a.mm_corr_available:
                g_b -= e_b.g_corr
                stripped.append(key_b)

            return g_a, g_b, stripped

        # ---- compute thermodynamic quantities ----------------------
        result: dict = {
            "anchor_file":  anchor.filename,
            "skeleton":     anchor.skeleton,
            "states":       states,
            "E_ox":         None,
            "E_red":        None,
            "E_pcet_red":   None,
            "E_pcet_ox":    None,
            "pKa_MH_plus":  None,
            "pKa_M":        None,
            # sp_only tracks which quantities used asymmetric thermal levels
            "sp_only":      {},
        }

        if _g("M") is None:
            return result

        pairs = [
            ("E_ox",       "M_ox",        "M",           False),
            ("E_red",      "M",           "M_red",       False),
            ("E_pcet_red", "M",           "MH_rad",      True),
            ("E_pcet_ox",  "M_deprot_rad","M",           True),
            ("pKa_MH_plus","MH_plus",     "M",           None),
            ("pKa_M",      "M",           "M_deprot",    None),
        ]

        for prop, key_a, key_b, is_pcet in pairs:
            g_a, g_b, stripped = _normalise_pair(key_a, key_b)
            if g_a is None:
                continue
            if stripped:
                result["sp_only"][prop] = stripped
                log.warning(
                    "Fetch thermal normalisation for %s: stripped g_corr "
                    "from %s (Hessian had failed for partner state).",
                    prop, stripped,
                )
            if prop == "E_ox":
                result[prop] = _reduction_potential(g_a, g_b, pH=self.pH)
            elif prop == "E_red":
                result[prop] = _reduction_potential(g_a, g_b, pH=self.pH)
            elif prop == "E_pcet_red":
                result[prop] = _reduction_potential(g_a, g_b, m=1, pH=self.pH)
            elif prop == "E_pcet_ox":
                result[prop] = _reduction_potential(g_a, g_b, m=1, pH=self.pH)
            elif prop == "pKa_MH_plus":
                result[prop] = _pka(g_acid=g_a, g_base=g_b)
            elif prop == "pKa_M":
                result[prop] = _pka(g_acid=g_a, g_base=g_b)

        return result

    # ------------------------------------------------------------------
    # Batch export
    # ------------------------------------------------------------------

    def batch_summary(
        self,
        input_csv: str,
        output_csv: str = "redox_results_summary.csv",
        solvation_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Process every molecule in ``input_csv`` and write a wide-format
        summary CSV with one row per molecule.

        Parameters
        ----------
        input_csv : str
            Path to CSV with columns: ``smiles``, ``charge``,
            ``multiplicity``.
        output_csv : str
            Path for the output summary CSV.
        solvation_name : str or None
            Named solvation to use for G_total (see ``fetch()``).

        Returns
        -------
        pd.DataFrame  — the summary table (also saved to ``output_csv``).

        Output columns
        --------------
        smiles, charge, multiplicity, anchor_file, skeleton_smiles,
        E_ox_V, E_red_V, E_pcet_red_V, E_pcet_ox_V, pKa_MH_plus, pKa_M,
        sp_only_flags,
        mm_ok_M, mm_ok_M_ox, mm_ok_M_red, mm_ok_MH_plus, mm_ok_MH_rad,
        mm_ok_M_deprot, mm_ok_M_deprot_rad,
        file_M, file_M_ox, file_M_red, file_MH_plus, file_MH_rad,
        file_M_deprot, file_M_deprot_rad
        """
        self._ensure_index()

        df = pd.read_csv(input_csv)
        required = {"smiles", "charge", "multiplicity"}
        missing  = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Input CSV is missing required columns: {missing}"
            )

        rows = []
        n_found = n_missing = 0
        print(f"Processing {len(df)} molecules from {input_csv} ...")

        for _, row in df.iterrows():
            smi  = str(row["smiles"]).strip()
            q    = int(row["charge"])
            m    = int(row["multiplicity"])

            res = self.fetch(smi, charge=q, mult=m,
                             solvation_name=solvation_name)

            if res is None:
                log.info("No match: %s", smi)
                n_missing += 1
                rows.append({
                    "smiles": smi, "charge": q, "multiplicity": m,
                    "anchor_file": None, "skeleton_smiles": None,
                    "E_ox_V": None, "E_red_V": None,
                    "E_pcet_red_V": None, "E_pcet_ox_V": None,
                    "pKa_MH_plus": None, "pKa_M": None,
                    "sp_only_flags": None,
                    **{f"mm_ok_{k}": None for k in (
                        "M", "M_ox", "M_red", "MH_plus",
                        "MH_rad", "M_deprot", "M_deprot_rad"
                    )},
                    **{f"file_{k}": None for k in (
                        "M", "M_ox", "M_red", "MH_plus",
                        "MH_rad", "M_deprot", "M_deprot_rad"
                    )},
                })
                continue

            n_found += 1
            states   = res["states"]
            sp_only  = res.get("sp_only", {})

            # Summarise SP-ONLY flags as a semicolon-separated string
            sp_flags = "; ".join(
                f"{prop}(Hessian failed: {', '.join(v)})"
                for prop, v in sp_only.items()
            ) if sp_only else ""

            rows.append({
                "smiles":           smi,
                "charge":           q,
                "multiplicity":     m,
                "anchor_file":      res["anchor_file"],
                "skeleton_smiles":  res["skeleton"],
                "E_ox_V":           _round(res["E_ox"]),
                "E_red_V":          _round(res["E_red"]),
                "E_pcet_red_V":     _round(res["E_pcet_red"]),
                "E_pcet_ox_V":      _round(res["E_pcet_ox"]),
                "pKa_MH_plus":      _round(res["pKa_MH_plus"]),
                "pKa_M":            _round(res["pKa_M"]),
                "sp_only_flags":    sp_flags,
                # MM correction quality per state
                **{f"mm_ok_{k}": (states[k].mm_corr_available if k in states else None)
                   for k in ("M", "M_ox", "M_red", "MH_plus",
                             "MH_rad", "M_deprot", "M_deprot_rad")},
                # File provenance
                **{f"file_{k}": states[k].filename if k in states else None
                   for k in ("M", "M_ox", "M_red", "MH_plus",
                             "MH_rad", "M_deprot", "M_deprot_rad")},
            })

        out_df = pd.DataFrame(rows)
        out_df.to_csv(output_csv, index=False)
        print(
            f"Done — {n_found} matched, {n_missing} not found. "
            f"Summary saved to '{output_csv}'."
        )
        return out_df

    # ------------------------------------------------------------------
    # Convenience: pretty-print a single result
    # ------------------------------------------------------------------

    def print_result(self, result: dict) -> None:
        """Pretty-print the output of a single ``fetch()`` call."""
        if result is None:
            print("No result.")
            return

        sep = "-" * 56
        print(f"\n{sep}")
        print(f"  Anchor : {result['anchor_file']}")
        print(f"  SMILES : {result['skeleton']}")
        print(sep)

        sp_only = result.get("sp_only", {})

        def _flag(prop: str) -> str:
            if prop not in sp_only:
                return ""
            states = ", ".join(sp_only[prop])
            return f"  [!] SP-ONLY: Hessian failed for {states}"

        pot_map = [
            ("E_ox",       "Oxidation   M -> M+    "),
            ("E_red",      "Reduction   M -> M-    "),
            ("E_pcet_red", "PCET Red    M -> MH    "),
            ("E_pcet_ox",  "PCET Ox     M -> M_dep "),
        ]
        pka_map = [
            ("pKa_MH_plus", "pKa  MH+ -> M    + H+"),
            ("pKa_M",       "pKa  M   -> M_dep + H+"),
        ]

        any_pot = False
        for key, lbl in pot_map:
            val = result.get(key)
            if val is not None:
                if not any_pot:
                    print("\n  Reduction potentials (V vs SHE):")
                    any_pot = True
                print(f"    {lbl}  {val:+8.3f} V{_flag(key)}")

        any_pka = False
        for key, lbl in pka_map:
            val = result.get(key)
            if val is not None:
                if not any_pka:
                    print("\n  pKa values:")
                    any_pka = True
                print(f"    {lbl}  {val:7.2f}{_flag(key)}")

        # Per-state quality notes
        print("\n  State quality:")
        for state, entry in result["states"].items():
            note = ""
            if not entry.mm_corr_available:
                note = "  [!] Hessian failed — SP only"
            print(f"    {state:14s}  {entry.filename}{note}")
        print(sep + "\n")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _round(val: Optional[float], decimals: int = 3) -> Optional[float]:
    """Round a float or return None."""
    return round(val, decimals) if val is not None else None


# ---------------------------------------------------------------------------
# CLI convenience
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Redox results fetcher")
    parser.add_argument("--folder",  default="redox_results",
                        help="HDF5 output folder")
    parser.add_argument("--csv",     required=True,
                        help="Input CSV with smiles/charge/multiplicity")
    parser.add_argument("--out",     default="redox_results_summary.csv",
                        help="Output summary CSV")
    parser.add_argument("--pH",      type=float, default=7.0)
    parser.add_argument("--solvation", default=None,
                        help="Named solvation to use (e.g. 'MD')")
    args = parser.parse_args()

    fetcher = RedoxFetcher(args.folder, pH=args.pH)
    fetcher.batch_summary(args.csv, args.out, solvation_name=args.solvation)
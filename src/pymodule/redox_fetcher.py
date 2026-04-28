from __future__ import annotations

import ast
import logging
import os
import re
from dataclasses import dataclass, field
from typing import Optional

import h5py
import networkx as nx
import pandas as pd
from networkx.algorithms import isomorphism
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Import physical constants from the calculator so there is one source of truth.
# Support both package and local imports.
try:
    from veloxchem.redox_calculator import (
        G_PROTON_AQ,
        G_ELECTRON,
        CONV_EV,
        AU_TO_KCAL,
        E_REF_SHE,
        RT_LN10_KCAL,
        NERNST_FACTOR,
    )
except ImportError:
    try:
        from redox_calculator import (
            G_PROTON_AQ,
            G_ELECTRON,
            CONV_EV,
            AU_TO_KCAL,
            E_REF_SHE,
            RT_LN10_KCAL,
            NERNST_FACTOR,
        )
    except ImportError:
        G_PROTON_AQ = -0.43676
        G_ELECTRON = -0.0015
        CONV_EV = 27.2114
        AU_TO_KCAL = 627.509
        E_REF_SHE = 4.28
        RT_LN10_KCAL = 1.364
        NERNST_FACTOR = 0.05916

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Index entry — one per HDF5 file
# ---------------------------------------------------------------------------

@dataclass
class _H5Entry:
    """Lightweight record for one HDF5 result file."""
    filename: str
    g_total: float
    g_corr: float
    charge: int
    mult: int
    h_count: int
    heavy_form: str
    graph: Optional[nx.Graph]
    skeleton: str
    mm_corr_available: bool = True
    solvation_data: dict = field(default_factory=dict)

    @property
    def full_formula(self) -> str:
        return f"{self.heavy_form}H{self.h_count}" if self.h_count else self.heavy_form


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _heavy_formula(formula_str: str) -> str:
    """
    Extract sorted heavy-atom counts from a molecular formula string.
    Example: 'C14H8N2O6' -> 'C14N2O6'
    """
    counts: dict[str, int] = {}
    for elem, count in re.findall(r"([A-Z][a-z]*)(\d*)", formula_str):
        if elem == "H":
            continue
        counts[elem] = counts.get(elem, 0) + (int(count) if count else 1)
    return "".join(f"{k}{counts[k]}" for k in sorted(counts))


def _formula_h_count(formula_str: str) -> int:
    """
    Extract the hydrogen count from a formula string.
    Example: 'C2H4O2' -> 4
    """
    m = re.search(r"H(\d*)", formula_str)
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


def _neutral_formula(mol) -> str:
    """
    Molecular formula from an RDKit Mol, with charge suffix stripped.
    Example: 'C6H5O2+' -> 'C6H5O2'
    """
    if mol is None:
        return ""
    return re.sub(r"[+\-]\d*$", "", CalcMolFormula(mol))


def _mol_to_graph(mol) -> Optional[nx.Graph]:
    """
    Heavy-atom connectivity graph from an RDKit Mol.

    Bond order is intentionally ignored so the graph matches the HDF5-side
    representation built from atom_symbols + bond_edges.
    """
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


def _read_text_dataset(h5file, key: str) -> Optional[str]:
    """Read a scalar text dataset from HDF5."""
    if key not in h5file:
        return None
    try:
        value = h5file[key][()]
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return str(value)
    except Exception:
        return None


def _graph_from_symbols_and_edges(
    atom_symbols_text: str,
    bond_edges_text: str,
) -> tuple[Optional[nx.Graph], int]:
    """
    Build a heavy-atom graph from atom symbols and bond-edge list.

    Parameters
    ----------
    atom_symbols_text : str
        Comma-separated atom symbols in atom-index order, e.g.
        "C,C,O,O,H,H,H,H"
    bond_edges_text : str
        Python-literal list of 0-based bonded pairs, e.g.
        "[(0, 1), (1, 2), (1, 3)]"

    Returns
    -------
    (graph, h_count)
    """
    try:
        atom_symbols = [s.strip() for s in atom_symbols_text.split(",") if s.strip()]
        edges = ast.literal_eval(bond_edges_text)

        if not isinstance(edges, list):
            raise ValueError("bond_edges is not a list")

        G = nx.Graph()
        heavy_map: dict[int, int] = {}
        next_idx = 0

        for abs_idx, elem in enumerate(atom_symbols):
            if elem != "H":
                heavy_map[abs_idx] = next_idx
                G.add_node(next_idx, element=elem)
                next_idx += 1

        for pair in edges:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                continue

            i = int(pair[0])
            j = int(pair[1])

            if i < 0 or j < 0 or i >= len(atom_symbols) or j >= len(atom_symbols):
                continue

            if atom_symbols[i] != "H" and atom_symbols[j] != "H":
                G.add_edge(heavy_map[i], heavy_map[j])

        h_count = sum(1 for elem in atom_symbols if elem == "H")
        return G, h_count

    except Exception as exc:
        log.warning("Could not parse atom_symbols/bond_edges: %s", exc)
        return None, 0


def _formula_atom_sequence(formula_str: str) -> list[str]:
    """
    Expand a formula into an atom list in formula token order.

    Example:
        C2H4O2 -> ['C', 'C', 'H', 'H', 'H', 'H', 'O', 'O']

    This is only used as a legacy fallback and is not reliable if the edge-list
    atom ordering did not follow formula order.
    """
    atoms: list[str] = []
    for elem, count in re.findall(r"([A-Z][a-z]*)(\d*)", formula_str):
        n = int(count) if count else 1
        atoms.extend([elem] * n)
    return atoms


def _legacy_skeleton_to_graph_and_hcount(raw: str) -> tuple[Optional[nx.Graph], int]:
    """
    Legacy fallback for old files where skeleton_smiles contains:
        FORMULA|[(i, j), ...]

    Warning:
        This fallback assumes atom indices follow the formula-expanded order.
        That is not always true for mixed-element molecules, so old files may
        still be ambiguous without atom_symbols + bond_edges.
    """
    try:
        formula_part, edges_part = raw.split("|", 1)
        atom_seq = _formula_atom_sequence(formula_part)
        edges = ast.literal_eval(edges_part)

        if not isinstance(edges, list):
            raise ValueError("Legacy edge list is not a list")

        G = nx.Graph()
        heavy_map: dict[int, int] = {}
        next_idx = 0

        for abs_idx, elem in enumerate(atom_seq):
            if elem != "H":
                heavy_map[abs_idx] = next_idx
                G.add_node(next_idx, element=elem)
                next_idx += 1

        for pair in edges:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                continue

            i = int(pair[0])
            j = int(pair[1])

            if i < 0 or j < 0 or i >= len(atom_seq) or j >= len(atom_seq):
                continue

            if atom_seq[i] != "H" and atom_seq[j] != "H":
                G.add_edge(heavy_map[i], heavy_map[j])

        h_count = sum(1 for elem in atom_seq if elem == "H")
        return G, h_count

    except Exception as exc:
        log.warning("Could not parse legacy skeleton key '%s': %s", raw, exc)
        return None, 0


# ---------------------------------------------------------------------------
# Thermochemistry
# ---------------------------------------------------------------------------

def _reduction_potential(
    g_ox: float,
    g_red: float,
    n: int = 1,
    m: int = 0,
    pH: float = 7.0,
) -> float:
    """E (V vs SHE) for Ox + n e- + m H+ -> Red."""
    dg_au = g_red - (g_ox + m * G_PROTON_AQ + n * G_ELECTRON)
    e_abs = -dg_au * CONV_EV / n
    return (e_abs - E_REF_SHE) - (NERNST_FACTOR / n) * m * pH


def _pka(g_acid: float, g_base: float) -> float:
    """pKa for HA -> A- + H+. g_acid = G(HA), g_base = G(A-)."""
    dg_kcal = ((g_base + G_PROTON_AQ) - g_acid) * AU_TO_KCAL
    return dg_kcal / RT_LN10_KCAL


# ---------------------------------------------------------------------------
# Main fetcher class
# ---------------------------------------------------------------------------

class RedoxFetcher:
    """
    Load and query results from a RedoxCalculator output folder.

    Matching strategy
    -----------------
    Query side:
        SMILES -> RDKit Mol -> heavy-atom graph
    HDF5 side:
        atom_symbols + bond_edges -> heavy-atom graph

    Notes
    -----
    - Bond order is ignored in graph matching.
    - geometry_xyz is intentionally not used.
    - Old files that only contain legacy skeleton_smiles may still work for
      simple molecules, but mixed-element molecules can be ambiguous.
    """

    def __init__(self, output_folder: str = "redox_results", pH: float = 7.0):
        self.output_folder = output_folder
        self.pH = pH
        self._index: list[_H5Entry] = []

    # ------------------------------------------------------------------
    # Index management
    # ------------------------------------------------------------------

    def _ensure_index(self) -> None:
        if not self._index:
            self.reload()

    def reload(self) -> None:
        """
        Re-scan output_folder and rebuild the in-memory index.
        """
        self._index = []

        if not os.path.isdir(self.output_folder):
            log.warning("Output folder not found: %s", self.output_folder)
            return

        n_loaded = 0
        n_skipped = 0

        for fname in os.listdir(self.output_folder):
            if not fname.endswith(".h5"):
                continue

            path = os.path.join(self.output_folder, fname)

            try:
                with h5py.File(path, "r") as f:
                    g_total = float(f["g_total"][()])
                    g_corr = float(f["g_corr"][()]) if "g_corr" in f else 0.0
                    charge = int(f["charge"][()])
                    mult = int(f["multiplicity"][()])

                    raw_skeleton = _read_text_dataset(f, "skeleton_smiles")
                    atom_symbols = _read_text_dataset(f, "atom_symbols")
                    bond_edges = _read_text_dataset(f, "bond_edges")

                    mm_ok = bool(f["mm_corr_available"][()]) if "mm_corr_available" in f else True
                    solv_data = _read_solvation_data(f)

                fname_stem = fname.split("_q")[0]
                heavy_form = _heavy_formula(fname_stem)

                # Preferred path: typed connectivity
                if atom_symbols and bond_edges:
                    graph, h_cnt = _graph_from_symbols_and_edges(atom_symbols, bond_edges)
                else:
                    # Legacy fallback for older files
                    if raw_skeleton and "|" in raw_skeleton:
                        graph, h_cnt = _legacy_skeleton_to_graph_and_hcount(raw_skeleton)
                    else:
                        graph, h_cnt = None, 0

                if h_cnt == 0:
                    h_cnt = _formula_h_count(fname_stem)

                skeleton_display = raw_skeleton or fname_stem

                self._index.append(
                    _H5Entry(
                        filename=fname,
                        g_total=g_total,
                        g_corr=g_corr,
                        charge=charge,
                        mult=mult,
                        h_count=h_cnt,
                        heavy_form=heavy_form,
                        graph=graph,
                        skeleton=skeleton_display,
                        mm_corr_available=mm_ok,
                        solvation_data=solv_data,
                    )
                )
                n_loaded += 1

            except Exception as exc:
                log.warning("Skipping %s: %s", fname, exc)
                n_skipped += 1

        log.info(
            "Index built: %d entries loaded, %d skipped from %s.",
            n_loaded,
            n_skipped,
            self.output_folder,
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
        Find the anchor state for query_smiles and assemble all related
        charge/protonation states into a thermodynamic result dict.
        """
        self._ensure_index()

        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            log.error("Could not parse SMILES: %s", query_smiles)
            return None

        query_graph = _mol_to_graph(query_mol)
        query_formula = _neutral_formula(query_mol)
        query_heavy = _heavy_formula(query_formula)
        query_h = _h_count_from_mol(query_mol)

        anchors: list[_H5Entry] = []

        for e in self._index:
            if e.charge != charge:
                continue
            if e.mult != mult:
                continue
            if e.heavy_form != query_heavy:
                continue
            if e.h_count != query_h:
                continue
            if e.graph is None or query_graph is None:
                continue
            if not _graphs_isomorphic(e.graph, query_graph):
                continue
            anchors.append(e)

        if not anchors:
            log.debug("No anchor found for %s (q=%d, m=%d).", query_smiles, charge, mult)
            return None

        if len(anchors) > 1:
            log.warning(
                "%d anchors matched %s — using lowest G_total.",
                len(anchors),
                query_smiles,
            )

        anchor = min(anchors, key=lambda e: e.g_total)

        states: dict[str, _H5Entry] = {"M": anchor}

        for e in self._index:
            if e.filename == anchor.filename:
                continue
            if e.heavy_form != anchor.heavy_form:
                continue
            if e.graph is None or anchor.graph is None:
                continue
            if not _graphs_isomorphic(e.graph, anchor.graph):
                continue

            dh = e.h_count - anchor.h_count
            dq = e.charge - anchor.charge

            key: Optional[str] = None
            if dh == 0 and dq == +1:
                key = "M_ox"
            elif dh == 0 and dq == -1:
                key = "M_red"
            elif dh == +1 and dq == +1:
                key = "MH_plus"
            elif dh == +1 and dq == 0:
                key = "MH_rad"
            elif dh == -1 and dq == -1:
                key = "M_deprot"
            elif dh == -1 and dq == 0:
                key = "M_deprot_rad"

            if key:
                if key not in states or e.g_total < states[key].g_total:
                    states[key] = e

        def _g(state_key: str) -> Optional[float]:
            e = states.get(state_key)
            if e is None:
                return None
            if solvation_name:
                if solvation_name not in e.solvation_data:
                    log.warning(
                        "Solvation '%s' not found for state '%s' (%s). "
                        "Falling back to default g_total.",
                        solvation_name,
                        state_key,
                        e.filename,
                    )
                    return e.g_total
                return e.solvation_data[solvation_name]
            return e.g_total

        def _normalise_pair(key_a: str, key_b: str) -> tuple:
            """
            Return (g_a, g_b, stripped_states) normalised to the same level.
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

        result: dict = {
            "anchor_file": anchor.filename,
            "skeleton": anchor.skeleton,
            "states": states,
            "E_ox": None,
            "E_red": None,
            "E_pcet_red": None,
            "E_pcet_ox": None,
            "pKa_MH_plus": None,
            "pKa_M": None,
            "sp_only": {},
        }

        if _g("M") is None:
            return result

        pairs = [
            ("E_ox", "M_ox", "M", False),
            ("E_red", "M", "M_red", False),
            ("E_pcet_red", "M", "MH_rad", True),
            ("E_pcet_ox", "M_deprot_rad", "M", True),
            ("pKa_MH_plus", "MH_plus", "M", None),
            ("pKa_M", "M", "M_deprot", None),
        ]

        for prop, key_a, key_b, _is_pcet in pairs:
            g_a, g_b, stripped = _normalise_pair(key_a, key_b)
            if g_a is None:
                continue

            if stripped:
                result["sp_only"][prop] = stripped
                log.warning(
                    "Fetch thermal normalisation for %s: stripped g_corr from %s "
                    "(Hessian had failed for partner state).",
                    prop,
                    stripped,
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
        Process every molecule in input_csv and write a wide-format summary CSV.
        """
        self._ensure_index()

        df = pd.read_csv(input_csv)
        required = {"smiles", "charge", "multiplicity"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {missing}")

        rows = []
        n_found = 0
        n_missing = 0
        print(f"Processing {len(df)} molecules from {input_csv} ...")

        for _, row in df.iterrows():
            smi = str(row["smiles"]).strip()
            q = int(row["charge"])
            m = int(row["multiplicity"])

            res = self.fetch(smi, charge=q, mult=m, solvation_name=solvation_name)

            if res is None:
                log.info("No match: %s", smi)
                n_missing += 1
                rows.append(
                    {
                        "smiles": smi,
                        "charge": q,
                        "multiplicity": m,
                        "anchor_file": None,
                        "skeleton_smiles": None,
                        "E_ox_V": None,
                        "E_red_V": None,
                        "E_pcet_red_V": None,
                        "E_pcet_ox_V": None,
                        "pKa_MH_plus": None,
                        "pKa_M": None,
                        "sp_only_flags": None,
                        **{
                            f"mm_ok_{k}": None
                            for k in (
                                "M",
                                "M_ox",
                                "M_red",
                                "MH_plus",
                                "MH_rad",
                                "M_deprot",
                                "M_deprot_rad",
                            )
                        },
                        **{
                            f"file_{k}": None
                            for k in (
                                "M",
                                "M_ox",
                                "M_red",
                                "MH_plus",
                                "MH_rad",
                                "M_deprot",
                                "M_deprot_rad",
                            )
                        },
                    }
                )
                continue

            n_found += 1
            states = res["states"]
            sp_only = res.get("sp_only", {})

            sp_flags = (
                "; ".join(
                    f"{prop}(Hessian failed: {', '.join(v)})"
                    for prop, v in sp_only.items()
                )
                if sp_only
                else ""
            )

            rows.append(
                {
                    "smiles": smi,
                    "charge": q,
                    "multiplicity": m,
                    "anchor_file": res["anchor_file"],
                    "skeleton_smiles": res["skeleton"],
                    "E_ox_V": _round(res["E_ox"]),
                    "E_red_V": _round(res["E_red"]),
                    "E_pcet_red_V": _round(res["E_pcet_red"]),
                    "E_pcet_ox_V": _round(res["E_pcet_ox"]),
                    "pKa_MH_plus": _round(res["pKa_MH_plus"]),
                    "pKa_M": _round(res["pKa_M"]),
                    "sp_only_flags": sp_flags,
                    **{
                        f"mm_ok_{k}": (states[k].mm_corr_available if k in states else None)
                        for k in (
                            "M",
                            "M_ox",
                            "M_red",
                            "MH_plus",
                            "MH_rad",
                            "M_deprot",
                            "M_deprot_rad",
                        )
                    },
                    **{
                        f"file_{k}": states[k].filename if k in states else None
                        for k in (
                            "M",
                            "M_ox",
                            "M_red",
                            "MH_plus",
                            "MH_rad",
                            "M_deprot",
                            "M_deprot_rad",
                        )
                    },
                }
            )

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
        """Pretty-print the output of a single fetch() call."""
        if result is None:
            print("No result.")
            return

        sep = "-" * 56
        print(f"\n{sep}")
        print(f"  Anchor : {result['anchor_file']}")
        print(f"  Key    : {result['skeleton']}")
        print(sep)

        sp_only = result.get("sp_only", {})

        def _flag(prop: str) -> str:
            if prop not in sp_only:
                return ""
            states = ", ".join(sp_only[prop])
            return f"  [!] SP-ONLY: Hessian failed for {states}"

        pot_map = [
            ("E_ox", "Oxidation   M -> M+    "),
            ("E_red", "Reduction   M -> M-    "),
            ("E_pcet_red", "PCET Red    M -> MH    "),
            ("E_pcet_ox", "PCET Ox     M -> M_dep "),
        ]
        pka_map = [
            ("pKa_MH_plus", "pKa  MH+ -> M    + H+"),
            ("pKa_M", "pKa  M   -> M_dep + H+"),
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
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Redox results fetcher")
    parser.add_argument("--folder", default="redox_results", help="HDF5 output folder")
    parser.add_argument("--csv", required=True, help="Input CSV with smiles/charge/multiplicity")
    parser.add_argument(
        "--out",
        default="redox_results_summary.csv",
        help="Output summary CSV",
    )
    parser.add_argument("--pH", type=float, default=7.0)
    parser.add_argument(
        "--solvation",
        default=None,
        help="Named solvation to use (e.g. 'MD')",
    )
    args = parser.parse_args()

    fetcher = RedoxFetcher(args.folder, pH=args.pH)
    fetcher.batch_summary(args.csv, args.out, solvation_name=args.solvation)
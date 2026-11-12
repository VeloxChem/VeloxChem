from __future__ import annotations
import ast
import logging
import math
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
        G_PROTON_AQ = -0.43073
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
    scf_solvent: float
    scf_vacuum: float
    charge: int
    mult: int
    h_count: int
    heavy_form: str
    graph: Optional[nx.Graph]
    skeleton: str
    mm_corr_available: bool = True
    solvation_data: dict = field(default_factory=dict)

    def elec_energy(self, solvation_name: Optional[str] = None) -> float:
        """
        Electronic energy in solvent, without the Gibbs correction.
        Default: scf_solvent (SMD single point). With a named solvation:
            scf_vacuum + solvation_data[name]
        (as in GibbsResult.g_total_with(), minus g_corr). Falls back to
        g_total - g_corr for files that do not store scf_solvent.
        """
        if (
            solvation_name
            and solvation_name in self.solvation_data
            and not math.isnan(self.scf_vacuum)
        ):
            return self.scf_vacuum + self.solvation_data[solvation_name]
        if not math.isnan(self.scf_solvent):
            return self.scf_solvent
        return self.g_total - self.g_corr

    def energy(self, solvation_name: Optional[str] = None) -> float:
        """G at the requested solvation level: elec_energy() + g_corr."""
        return self.elec_energy(solvation_name) + self.g_corr

    def selection_key(self, solvation_name: Optional[str] = None) -> tuple:
        """
        Sort key for choosing between duplicate files of the same state:
        the lowest electronic energy in solvent (scf_solvent) wins,
        regardless of whether the Gibbs correction succeeded.
        """
        return (self.elec_energy(solvation_name),)

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
    return re.sub(r"[+-]\d*$", "", CalcMolFormula(mol))

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
        G.add_node(
            atom.GetIdx(),
            element=atom.GetSymbol(),
            nH=atom.GetTotalNumHs(),
            idx=atom.GetIdx(),
        )
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G

def _graphs_isomorphic(g1: nx.Graph, g2: nx.Graph) -> bool:
    """True if g1 and g2 are graph-isomorphic (same elements, same topology)."""
    if len(g1) != len(g2):
        return False
    nm = isomorphism.categorical_node_match("element", None)
    return isomorphism.GraphMatcher(g1, g2, node_match=nm).is_isomorphic()
# Upper bound on skeleton mappings tried when looking for the one with the
# fewest H differences (symmetric skeletons have several equivalent mappings).
_MAX_H_MAPPINGS = 2000

def _h_diff(g_ref: nx.Graph, g: nx.Graph) -> Optional[list[tuple]]:
    """
    Compare where the H atoms sit on two skeletons with the same heavy-atom
    connectivity.
    Returns None if the heavy-atom graphs are not isomorphic. Otherwise
    returns the smallest list of (ref_node, n_H_ref, n_H_other) over all
    skeleton mappings; [] means the H pattern is identical.
    """
    if len(g_ref) != len(g):
        return None
    nm_h = isomorphism.categorical_node_match(["element", "nH"], [None, None])
    if isomorphism.GraphMatcher(g_ref, g, node_match=nm_h).is_isomorphic():
        return []
    nm = isomorphism.categorical_node_match("element", None)
    best: Optional[list[tuple]] = None
    for n_tried, mapping in enumerate(
        isomorphism.GraphMatcher(g_ref, g, node_match=nm).isomorphisms_iter()
    ):
        diff = [
            (a, g_ref.nodes[a].get("nH", 0), g.nodes[b].get("nH", 0))
            for a, b in mapping.items()
            if g_ref.nodes[a].get("nH", 0) != g.nodes[b].get("nH", 0)
        ]
        if best is None or len(diff) < len(best):
            best = diff
        if len(best) <= 1 or n_tried + 1 >= _MAX_H_MAPPINGS:
            break
    return best

def _h_diff_matches(diff: Optional[list[tuple]], dh: int) -> bool:
    """
    True if the H difference is what the state's definition implies:
    identical pattern for dh == 0, one H added/removed at a single atom
    for dh == +1/-1.
    """
    if diff is None:
        return False
    if dh == 0:
        return not diff
    return len(diff) == 1 and diff[0][2] - diff[0][1] == dh

def _describe_h_diff(g_ref: nx.Graph, diff: Optional[list[tuple]]) -> str:
    """Human-readable H difference, e.g. 'N11: H 2->1; O31: H 0->1'."""
    if diff is None:
        return "skeleton mismatch"
    return "; ".join(
        f"{g_ref.nodes[a]['element']}{g_ref.nodes[a].get('idx', a)}: H {x}->{y}"
        for a, x, y in sorted(diff, key=lambda t: g_ref.nodes[t[0]].get("idx", t[0]))
    )
_STATE_KEYS = (
    "M",
    "M_ox",
    "M_red",
    "MH_plus",
    "MH_rad",
    "M_deprot",
    "M_deprot_rad",
)

def _h_check_keys() -> list[str]:
    """All '<state>_h_ok' / '<state>_h_diff' result keys, in a stable order."""
    return [f"{k}_{sfx}" for k in _STATE_KEYS for sfx in ("h_ok", "h_diff")]

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
                G.add_node(next_idx, element=elem, nH=0, idx=abs_idx)
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
            elif atom_symbols[i] == "H" and atom_symbols[j] != "H":
                G.nodes[heavy_map[j]]["nH"] += 1
            elif atom_symbols[j] == "H" and atom_symbols[i] != "H":
                G.nodes[heavy_map[i]]["nH"] += 1
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
                G.add_node(next_idx, element=elem, nH=0, idx=abs_idx)
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
            elif atom_seq[i] == "H" and atom_seq[j] != "H":
                G.nodes[heavy_map[j]]["nH"] += 1
            elif atom_seq[j] == "H" and atom_seq[i] != "H":
                G.nodes[heavy_map[i]]["nH"] += 1
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
# Column suffixes for the two files each property is computed from:
#   potentials: (oxidised state, reduced state)
#   pKa:        (acid, base)
_PAIR_FILE_SUFFIXES = {
    "E_ox": ("ox", "red"),
    "E_red": ("ox", "red"),
    "E_pcet_red": ("ox", "red"),
    "E_pcet_ox": ("ox", "red"),
    "pKa_MH_plus": ("acid", "base"),
    "pKa_M": ("acid", "base"),
}

def _pair_file_keys() -> list[str]:
    """All '<prop>_file_<suffix>' keys, in a stable order."""
    return [
        f"{prop}_file_{sfx}"
        for prop, sfxs in _PAIR_FILE_SUFFIXES.items()
        for sfx in sfxs
    ]
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
                    has_g_corr = "g_corr" in f
                    g_corr = float(f["g_corr"][()]) if has_g_corr else 0.0
                    scf_solvent = (
                        float(f["scf_solvent"][()])
                        if "scf_solvent" in f
                        else float("nan")
                    )
                    scf_vacuum = (
                        float(f["scf_vacuum"][()])
                        if "scf_vacuum" in f
                        else float("nan")
                    )
                    charge = int(f["charge"][()])
                    mult = int(f["multiplicity"][()])
                    raw_skeleton = _read_text_dataset(f, "skeleton_smiles")
                    atom_symbols = _read_text_dataset(f, "atom_symbols")
                    bond_edges = _read_text_dataset(f, "bond_edges")
                    # A file without g_corr cannot have been Gibbs-corrected.
                    if "mm_corr_available" in f:
                        mm_ok = bool(f["mm_corr_available"][()]) and has_g_corr
                    else:
                        mm_ok = has_g_corr
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
                        scf_solvent=scf_solvent,
                        scf_vacuum=scf_vacuum,
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

    def list_compound_files(
        self,
        smiles: str,
        exact_h: bool = False,
        charge: Optional[int] = None,
        mult: Optional[int] = None,
    ) -> pd.DataFrame:
        """List every indexed H5 file sharing the query's heavy-atom skeleton.
        Different charge states, protonation states, and conformers are included.
        If exact_h=True, files with the same H count but different H positions
        are omitted; states with other H counts remain included.
        This method never selects only the lowest-energy conformer.
        """
        self._ensure_index()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Cannot parse SMILES: {smiles}")
        query_graph = _mol_to_graph(mol)
        query_heavy = _heavy_formula(_neutral_formula(mol))
        query_h = _h_count_from_mol(mol)
        columns = [
            "filename", "charge", "multiplicity", "H_count",
            "scf_solvent_au", "scf_vacuum_au", "g_corr_au",
            "g_total_au", "gibbs_available", "H_difference",
        ]
        rows = []
        for entry in self._index:
            if charge is not None and entry.charge != charge:
                continue
            if mult is not None and entry.mult != mult:
                continue
            if entry.heavy_form != query_heavy or entry.graph is None:
                continue
            diff = _h_diff(query_graph, entry.graph)
            if diff is None:
                continue
            if exact_h and entry.h_count == query_h and diff:
                continue
            rows.append({
                "filename": entry.filename,
                "charge": entry.charge,
                "multiplicity": entry.mult,
                "H_count": entry.h_count,
                "scf_solvent_au": entry.scf_solvent,
                "scf_vacuum_au": entry.scf_vacuum,
                "g_corr_au": entry.g_corr if entry.mm_corr_available else float("nan"),
                "g_total_au": entry.g_total if entry.mm_corr_available else float("nan"),
                "gibbs_available": entry.mm_corr_available,
                "H_difference": _describe_h_diff(query_graph, diff),
            })
        result = pd.DataFrame(rows, columns=columns)
        if not result.empty:
            result = result.sort_values(
                ["charge", "multiplicity", "H_count", "scf_solvent_au", "filename"],
                na_position="last",
            ).reset_index(drop=True)
        return result

    # ------------------------------------------------------------------
    # Core fetch

    # ------------------------------------------------------------------
    def fetch(
        self,
        query_smiles: str,
        charge: int,
        mult: int,
        solvation_name: Optional[str] = None,
        exact_h: bool = True,
    ) -> Optional[dict]:
        """
        Find the anchor state for query_smiles and assemble all related
        charge/protonation states into a thermodynamic result dict.
        With exact_h=True (default) only files whose H positions match are
        used: the anchor must have exactly the H pattern of query_smiles, and
        each related state must differ from the anchor by the expected H
        change (none, or one H at a single atom). Tautomers are skipped, so
        fetch() returns None if the exact structure was never computed.
        With exact_h=False the closest tautomer is used and flagged in the
        '<state>_h_ok' / '<state>_h_diff' result keys.
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
        anchor_diffs: dict[str, list[tuple]] = {}
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
            diff = _h_diff(query_graph, e.graph)
            if diff is None:
                continue
            anchors.append(e)
            anchor_diffs[e.filename] = diff
        if exact_h and anchors:
            tautomers = [e for e in anchors if anchor_diffs[e.filename]]
            anchors = [e for e in anchors if not anchor_diffs[e.filename]]
            if not anchors:
                log.warning(
                    "No exact match for %s (q=%d, m=%d); skipped tautomer(s): %s",
                    query_smiles,
                    charge,
                    mult,
                    "; ".join(
                        f"{e.filename} ({_describe_h_diff(query_graph, anchor_diffs[e.filename])})"
                        for e in tautomers
                    ),
                )
        if not anchors:
            log.debug("No anchor found for %s (q=%d, m=%d).", query_smiles, charge, mult)
            return None
        if len(anchors) > 1:
            log.warning(
                "%d anchors matched %s — using lowest scf_solvent "
                "(exact H pattern preferred).",
                len(anchors),
                query_smiles,
            )
        # Exact H pattern first, then lowest scf_solvent.
        anchor = min(
            anchors,
            key=lambda e: (bool(anchor_diffs[e.filename]), e.selection_key(solvation_name)),
        )
        anchor_diff = anchor_diffs[anchor.filename]
        h_ok: dict[str, bool] = {"M": not anchor_diff}
        h_desc: dict[str, str] = {"M": _describe_h_diff(query_graph, anchor_diff)}
        if anchor_diff:
            log.warning(
                "No file has the exact H pattern of %s (q=%d, m=%d); using tautomer "
                "%s (%s).",
                query_smiles,
                charge,
                mult,
                anchor.filename,
                h_desc["M"],
            )
        states: dict[str, _H5Entry] = {"M": anchor}
        best_key: dict[str, tuple] = {}
        for e in self._index:
            if e.filename == anchor.filename:
                continue
            if e.heavy_form != anchor.heavy_form:
                continue
            if e.graph is None or anchor.graph is None:
                continue
            dh = e.h_count - anchor.h_count
            if abs(dh) > 1:
                continue
            diff = _h_diff(anchor.graph, e.graph)
            if diff is None:
                continue
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
                ok = _h_diff_matches(diff, dh)
                if exact_h and not ok:
                    continue
                # Expected H change first, then lowest scf_solvent.
                cand = (not ok, e.selection_key(solvation_name))
                if key not in states or cand < best_key[key]:
                    states[key] = e
                    best_key[key] = cand
                    h_ok[key] = ok
                    h_desc[key] = _describe_h_diff(anchor.graph, diff)
        for key in states:
            if key != "M" and not h_ok[key]:
                log.warning(
                    "State %s for %s uses %s, whose H pattern is not the expected "
                    "change from the anchor (%s).",
                    key,
                    query_smiles,
                    states[key].filename,
                    h_desc[key],
                )
        def _g(state_key: str) -> Optional[float]:
            e = states.get(state_key)
            if e is None:
                return None
            if solvation_name:
                if solvation_name not in e.solvation_data:
                    log.warning(
                        "Solvation '%s' not found for state '%s' (%s). "
                        "Falling back to scf_solvent.",
                        solvation_name,
                        state_key,
                        e.filename,
                    )
                    return e.elec_energy()
                if math.isnan(e.scf_vacuum):
                    log.warning(
                        "No scf_vacuum stored for state '%s' (%s); cannot apply "
                        "solvation '%s'. Falling back to scf_solvent.",
                        state_key,
                        e.filename,
                        solvation_name,
                    )
            return e.elec_energy(solvation_name)
        def _normalise_pair(key_a: str, key_b: str) -> tuple:
            """
            Return (g_a, g_b, no_corr_states) on a common level.
            If both states have a Gibbs correction, g = E_solv + g_corr for
            both. Otherwise g_corr is left out for both (g = E_solv), and
            no_corr_states lists the state(s) whose Hessian failed.
            """
            e_a = states.get(key_a)
            e_b = states.get(key_b)
            if e_a is None or e_b is None:
                return None, None, []
            g_a = _g(key_a)
            g_b = _g(key_b)
            if e_a.mm_corr_available and e_b.mm_corr_available:
                return g_a + e_a.g_corr, g_b + e_b.g_corr, []
            no_corr = [k for k, e in ((key_a, e_a), (key_b, e_b)) if not e.mm_corr_available]
            return g_a, g_b, no_corr
        result: dict = {
            "anchor_file": anchor.filename,
            "anchor_exact_match": not anchor_diff,
            **{f"{k}_h_ok": h_ok.get(k) for k in _STATE_KEYS},
            **{f"{k}_h_diff": h_desc.get(k) for k in _STATE_KEYS},
            "skeleton": anchor.skeleton,
            "states": states,
            "E_ox": None,
            "E_red": None,
            "E_pcet_red": None,
            "E_pcet_ox": None,
            "pKa_MH_plus": None,
            "pKa_M": None,
            "sp_only": {},
            **{k: None for k in _pair_file_keys()},
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
            sfx_a, sfx_b = _PAIR_FILE_SUFFIXES[prop]
            result[f"{prop}_file_{sfx_a}"] = states[key_a].filename
            result[f"{prop}_file_{sfx_b}"] = states[key_b].filename
            if stripped:
                result["sp_only"][prop] = stripped
                log.warning(
                    "%s computed without g_corr (electronic energies only): "
                    "Hessian failed for %s.",
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
        exact_h: bool = True,
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
            res = self.fetch(
                smi, charge=q, mult=m, solvation_name=solvation_name, exact_h=exact_h
            )
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
                        "anchor_exact_match": None,
                        **{k: None for k in _h_check_keys()},
                        **{k: None for k in _pair_file_keys()},
                        **{
                            f"g_total_{k}_au": None
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
                            f"g_corr_{k}_au": None
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
                            f"scf_solvent_{k}_au": None
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
                            f"scf_vacuum_{k}_au": None
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
                    "anchor_exact_match": res["anchor_exact_match"],
                    **{k: res[k] for k in _h_check_keys()},
                    **{k: res[k] for k in _pair_file_keys()},
                    **{
                        f"g_total_{k}_au": (states[k].g_total if k in states else None)
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
                        f"g_corr_{k}_au": (states[k].g_corr if k in states else None)
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
                        f"scf_solvent_{k}_au": (states[k].scf_solvent if k in states else None)
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
                        f"scf_vacuum_{k}_au": (states[k].scf_vacuum if k in states else None)
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
        if result.get("anchor_exact_match") is False:
            print(f"  [!] Not the exact SMILES structure: {result.get('M_h_diff')}")
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
                print(f"        ox:  {result.get(key + '_file_ox')}")
                print(f"        red: {result.get(key + '_file_red')}")
        any_pka = False
        for key, lbl in pka_map:
            val = result.get(key)
            if val is not None:
                if not any_pka:
                    print("\n  pKa values:")
                    any_pka = True
                print(f"    {lbl}  {val:7.2f}{_flag(key)}")
                print(f"        acid: {result.get(key + '_file_acid')}")
                print(f"        base: {result.get(key + '_file_base')}")
        print("\n  State energies (Hartree):")
        for state, entry in result["states"].items():
            note = ""
            if not entry.mm_corr_available:
                note = "  [!] Hessian failed — SP only"
            if result.get(f"{state}_h_ok") is False:
                note += f"  [!] H pattern: {result.get(f'{state}_h_diff')}"
            print(f"    {state:14s}  {entry.filename}{note}")
            print(f"        G_total      = {entry.g_total:.12f}")
            print(f"        g_corr       = {entry.g_corr:.12f}")
            print(f"        scf_solvent  = {entry.scf_solvent:.12f}")
            print(f"        scf_vacuum   = {entry.scf_vacuum:.12f}")
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
    parser.add_argument("--csv", help="Input CSV with smiles/charge/multiplicity")
    parser.add_argument("--smiles", help="List all matching H5 files for this SMILES")
    parser.add_argument("--charge", type=int, help="Optional charge filter for --smiles")
    parser.add_argument("--mult", type=int, help="Optional multiplicity filter for --smiles")
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
    parser.add_argument(
        "--allow-tautomers",
        action="store_true",
        help="Use the closest tautomer (flagged) when no exact H-pattern match exists",
    )
    args = parser.parse_args()
    fetcher = RedoxFetcher(args.folder, pH=args.pH)
    if args.smiles:
        df = fetcher.list_compound_files(
            args.smiles,
            exact_h=not args.allow_tautomers,
            charge=args.charge,
            mult=args.mult,
        )
        if df.empty:
            print("No matching H5 files found.")
        else:
            print(df.to_string(index=False))
            df.to_csv(args.out, index=False)
            print(f"\nSaved {len(df)} matches to {args.out}")
    elif args.csv:
        fetcher.batch_summary(
            args.csv, args.out,
            solvation_name=args.solvation,
            exact_h=not args.allow_tautomers,
        )
    else:
        parser.error("Provide either --csv or --smiles")

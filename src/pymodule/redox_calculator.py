"""
redox_calculator.py
====================
Automated redox potential and pKa calculator using VeloxChem.

Gibbs free energy composite protocol
--------------------------------------
G_total = E_SP(solvent, basis_sp, xcfun_sp)
        + G_corr(MM)

where the thermal correction is:

    G_corr(MM) = G_vib(MM) - E_elec(MM)

The MM force field is built automatically by MMForceFieldGenerator using
GAFF atom types and fast electronegativity-equalisation (EEM) charges
(resp=False), so no extra DFT charge calculation is needed for the
correction step.  The Hessian is evaluated numerically by MMHessianDriver.

This replaces the previous xTB dependency entirely.

Physical chemistry conventions
--------------------------------
* All internal energies are in atomic units (Hartree).
* Reduction potentials referenced to SHE via the absolute electrode scale
  (E_ref = 4.28 V -- Trasatti 1986; Reuter & Scheffler 2001).
* Proton free energy in aqueous solution: G(H+) = -0.4390 au
  (Tissandier et al., J. Phys. Chem. A 1998, 102, 7787).
* Electron free energy correction: G(e-) = -0.0015 au
  (Kelly et al., J. Phys. Chem. B 2006, 110, 16066).
* pKa = dG(deprotonation) / (RT ln10),  RT ln10 = 1.364 kcal/mol at 298 K.
* PCET Nernst shift: dE = -(0.05916 / n) * m * pH  (25 degC, aqueous).

Multiplicity assignment
------------------------
Closed-shell molecules -> singlet (mult = 1).
Open-shell molecules   -> doublet (mult = 2).
First-row TM complexes -> warning emitted; verify manually.
"""

from __future__ import annotations

import hashlib
import logging
import os
import time
from collections import Counter
from dataclasses import dataclass
from typing import Optional

import h5py
import networkx as nx
import numpy as np
import pandas as pd
import veloxchem as vlx
from networkx.algorithms import isomorphism


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
G_PROTON_AQ: float   = -0.43676   # au   G(H+, aq, 1M) = -274.1 kcal/mol
                                   # Tissandier et al. J.Phys.Chem.A 1998 (-265.9 kcal/mol solvation)
                                   # + gas-phase H+ free energy (-6.28 kcal/mol)
                                   # + 1 atm -> 1 M standard state correction (-1.89 kcal/mol)
                                   # used for DFT/SMD pKa calculations.
                                   # (Truhlar/Cramer convention; widely used with
                                   # implicit solvation models including SMD)
G_ELECTRON: float    = -0.0015   # au   G(e-) correction, Kelly 2006
CONV_EV: float       = 27.2114   # Ha -> eV
AU_TO_KCAL: float    = 627.509   # Ha -> kcal/mol
E_REF_SHE: float     =  4.28    # V    absolute SHE,     Trasatti 1986
RT_LN10_KCAL: float  =  1.364   # kcal/mol  (RT ln10 at 298.15 K)
NERNST_FACTOR: float =  0.05916  # V per pH unit per electron (25 degC)

FIRST_ROW_TM: frozenset = frozenset(
    {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
)

# Covalent cutoffs for protic-site detection (Angstrom)
BOND_CUTOFF: dict[str, float] = {
    "O": 1.10, "N": 1.15, "S": 1.45, "P": 1.50
}


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass
class GibbsResult:
    """
    All computed quantities for one charge/protonation state.

    The SMD and vacuum single-point energies are always stored.  Named
    solvation free energies from external computations (FEP, COSMO-RS, ...)
    are stored in ``solvation_data`` and can be used to recompute ``G_total``
    on demand via ``g_total_with(label)``.

    Default G_total formula
    -----------------------
    ::

        G_total = scf_solvent + g_corr

    With a named external solvation energy swapped in
    -------------------------------------------------
    ::

        G_total(label) = scf_vacuum + solvation_data[label] + g_corr

    Attributes
    ----------
    g_total : float
        Default composite Gibbs free energy using SMD (au).
    charge : int
        Formal charge.
    multiplicity : int
        Spin multiplicity.
    scf_solvent : float
        DFT single-point energy in implicit solvent (au).
    scf_vacuum : float
        DFT single-point energy in vacuum (au).
    g_corr : float
        MM Gibbs thermal correction G_vib(MM) - E_MM (au).
    geometry_xyz : str
        XYZ block of the optimised geometry (Angstrom).
    is_fragmented:   bool  = False
        True if the molecule broke apart during optimisation.
    skeleton_smiles : str
        Charge- and H-stripped canonical SMILES.
    label : str
        Short identifier, e.g. "M", "M+", "MH+_site3".
    solvation_data : dict[str, float]
        Named solvation free energies in Hartree, keyed by a user-supplied
        label such as "MD", "COSMO-RS", "FEP".  Populated via
        ``add_solvation()`` or ``RedoxCalculator.register_solvation()``.
        Each entry represents a *total* G_solv to replace scf_solvent with,
        not a correction on top of it.
    """
    g_total:         float = 0.0
    charge:          int   = 0
    multiplicity:    int   = 1
    scf_solvent:     float = 0.0
    scf_vacuum:      float = 0.0
    g_corr:          float = 0.0   # xTB thermal correction: G(xTB) - E_elec(xTB)
    geometry_xyz:    str   = ""
    is_fragmented:   bool  = False
    skeleton_smiles: str   = ""
    atom_symbols:    str   = ""
    bond_edges:      str   = ""
    label:           str   = ""
    mm_corr_available: bool = True   # False if MM Gibbs correction failed
    # Populated after construction; not a dataclass field with a mutable default
    # so we initialise it in __post_init__.
    solvation_data:  dict  = None   # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.solvation_data is None:
            self.solvation_data = {}

    def add_solvation(self, name: str, dg_solv: float) -> None:
        """
        Store a named solvation free energy (Hartree) for this state.

        Parameters
        ----------
        name : str
            Arbitrary label, e.g. ``"MD"``, ``"COSMO-RS"``, ``"FEP"``.
        dg_solv : float
            Total solvation free energy in Hartree.
            Sign convention: negative = stabilised by solvent.
            This value *replaces* ``scf_solvent`` in the G_total formula;
            it is not added on top of it.
        """
        self.solvation_data[name] = dg_solv
        log.debug(
            "GibbsResult '%s': stored solvation_data['%s'] = %.6f au",
            self.label, name, dg_solv,
        )

    def g_total_with(self, name: str) -> float:
        """
        Return G_total with the named solvation free energy swapped in.

        Formula::

            G_total(name) = scf_vacuum + solvation_data[name] + g_corr

        Parameters
        ----------
        name : str
            Key in ``solvation_data``.

        Returns
        -------
        float : G_total in Hartree.

        Raises
        ------
        KeyError if ``name`` is not in ``solvation_data``.
        """
        if name not in self.solvation_data:
            available = list(self.solvation_data.keys())
            raise KeyError(
                f"Solvation label {name!r} not found for state {self.label!r}. "
                f"Available: {available}"
            )
        return self.scf_vacuum + self.solvation_data[name] + self.g_corr

    def as_dict(self) -> dict:
        """Return a flat dictionary suitable for HDF5 storage."""
        d = {k: v for k, v in self.__dict__.items() if k != "solvation_data"}
        # Store each named solvation value as a separate key
        for name, val in self.solvation_data.items():
            d[f"solvation_data_{name}"] = val
        return d


@dataclass
class RedoxProfile:
    """
    Complete thermodynamic profile for one molecule.

    GibbsResult objects and scalar properties are named attributes,
    so IDE auto-complete and ``vars(profile)`` work naturally.

    Charge-state results
    --------------------
    M            : neutral reference state
    M_ox         : one-electron oxidised  (M+)
    M_red        : one-electron reduced   (M-)
    MH_plus      : protonated, charge +1  (MH+)
    MH_rad       : protonated neutral radical (MH-)
    M_deprot     : deprotonated anion     (M_dep-)
    M_deprot_rad : deprotonated neutral radical (M_dep-)

    Derived scalars  (None if the required states were not computed)
    ----------------------------------------------------------------
    E_ox         : M -> M+       oxidation potential (V vs SHE)
    E_red        : M -> M-       reduction potential (V vs SHE)
    E_pcet_red   : M -> MH-      PCET reduction      (V vs SHE)
    E_pcet_ox    : M -> M_dep-   PCET oxidation      (V vs SHE)
    pKa_MH_plus  : MH+ -> M  + H+
    pKa_MH_rad   : MH- -> M- + H+
    pKa_M        : M   -> M_dep- + H+
    pKa_M_plus   : M+  -> M_dep- + H+
    """
    # Gibbs results
    M:             Optional[GibbsResult] = None
    M_ox:          Optional[GibbsResult] = None
    M_red:         Optional[GibbsResult] = None
    MH_plus:       Optional[GibbsResult] = None
    MH_rad:        Optional[GibbsResult] = None
    M_deprot:      Optional[GibbsResult] = None
    M_deprot_rad:  Optional[GibbsResult] = None
    

    # Reduction potentials (V vs SHE)
    E_ox:        Optional[float] = None
    E_red:       Optional[float] = None
    E_pcet_red:  Optional[float] = None
    E_pcet_ox:   Optional[float] = None

    # pKa values
    pKa_MH_plus: Optional[float] = None
    pKa_MH_rad:  Optional[float] = None
    pKa_M:       Optional[float] = None
    pKa_M_plus:  Optional[float] = None

    def apply_solvation(
        self,
        name: str,
        calc: "RedoxCalculator",
    ) -> "RedoxProfile":
        """
        Return a new RedoxProfile with G_total recomputed for every state
        using the named solvation free energy, and all potentials and pKa
        values derived from those updated energies.

        The original profile is not modified.

        For states where ``solvation_data[name]`` is missing the default
        ``g_total`` (SMD-based) is kept unchanged, and a warning is logged.

        Formula applied per state::

            G_total(name) = scf_vacuum + solvation_data[name] + g_corr

        Parameters
        ----------
        name : str
            Key in each ``GibbsResult.solvation_data``.
        calc : RedoxCalculator
            Used to recompute potentials and pKas with the correct pH and
            physical constants.  Typically the same instance that produced
            this profile.

        Returns
        -------
        RedoxProfile with updated ``g_total`` values on each GibbsResult
        and freshly derived ``E_ox``, ``E_red``, ``pKa_M``, etc.

        Example
        -------
        ::

            profile = calc.run("c1ccccc1O")

            # Add FEP values after the fact
            profile.M.add_solvation("MD", -0.0145)
            profile.M_ox.add_solvation("MD", -0.0089)
            profile.M_red.add_solvation("MD", -0.0201)

            profile_md = profile.apply_solvation("MD", calc)
            print(profile_md.E_red)
        """
        import copy

        new = copy.deepcopy(self)

        # Update g_total for every state that has the named solvation entry
        state_attrs = (
            "M", "M_ox", "M_red", "MH_plus", "MH_rad",
            "M_deprot", "M_deprot_rad",
        )
        for attr in state_attrs:
            res: Optional[GibbsResult] = getattr(new, attr)
            if res is None:
                continue
            if name in res.solvation_data:
                res.g_total = res.g_total_with(name)
            else:
                log.warning(
                    "apply_solvation(%r): state '%s' has no entry for '%s', "
                    "keeping default G_total (SMD-based).",
                    name, res.label, name,
                )

        # Recompute all derived quantities from the updated g_total values
        new.E_ox = new.E_red = new.E_pcet_red = new.E_pcet_ox = None
        new.pKa_MH_plus = new.pKa_MH_rad = new.pKa_M = new.pKa_M_plus = None

        if new.M:
            g_M = new.M.g_total
            if new.M_ox:
                new.E_ox = calc.reduction_potential(
                    g_ox=new.M_ox.g_total, g_red=g_M
                )
            if new.M_red:
                new.E_red = calc.reduction_potential(
                    g_ox=g_M, g_red=new.M_red.g_total
                )
            if new.MH_rad:
                new.E_pcet_red = calc.reduction_potential(
                    g_ox=g_M, g_red=new.MH_rad.g_total, n=1, m=1
                )
            if new.M_deprot_rad:
                new.E_pcet_ox = calc.reduction_potential(
                    g_ox=new.M_deprot_rad.g_total, g_red=g_M, n=1, m=1
                )
            if new.MH_plus:
                new.pKa_MH_plus = calc.pka(
                    g_acid=new.MH_plus.g_total, g_base=new.M.g_total
                )
            if new.MH_rad and new.M_red:
                new.pKa_MH_rad = calc.pka(
                    g_acid=new.MH_rad.g_total, g_base=new.M_red.g_total
                )
            if new.M_deprot:
                new.pKa_M = calc.pka(
                    g_acid=new.M.g_total, g_base=new.M_deprot.g_total
                )
            if new.M_ox and new.M_deprot_rad:
                new.pKa_M_plus = calc.pka(
                    g_acid=new.M_ox.g_total, g_base=new.M_deprot_rad.g_total
                )

        return new


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _mol_formula(mol_obj) -> str:
    """Hill-ordered molecular formula from a VeloxChem Molecule."""
    counts = Counter(mol_obj.get_labels())
    order = [el for el in ("C", "H") if el in counts]
    order += sorted(k for k in counts if k not in ("C", "H"))
    return "".join(
        f"{el}{counts[el] if counts[el] > 1 else ''}" for el in order
    )


def _xyz_block(mol_obj) -> str:
    """Serialise a VeloxChem Molecule to an XYZ string."""
    labels = mol_obj.get_labels()
    coords = mol_obj.get_coordinates_in_angstrom()
    lines = [str(len(labels)), ""]
    for sym, c in zip(labels, coords):
        lines.append(f"{sym:4s} {c[0]:16.8f} {c[1]:16.8f} {c[2]:16.8f}")
    return "\n".join(lines)


    """
    Simple molecule identifier from VeloxChem Molecule: Hill formula.
    Used for HDF5 provenance only -- not for comparison or deduplication.
    """
    return _mol_formula(mol_obj)


def _skeleton_smiles(input_data) -> str:
    """
    Skeleton identifier for a molecule -- heavy-atom formula + sorted
    adjacency list from geometry-based bond perception.

    Accepts a VeloxChem Molecule or an XYZ string.
    Used for HDF5 storage and graph-based deduplication.
    Returns a stable string key, not a SMILES.
    """
    try:
        if isinstance(input_data, str):
            mol = vlx.Molecule.read_xyz_string(input_data)
        else:
            mol = input_data
        formula = _mol_formula(mol)
        bonds   = sorted(_geometry_bonds(mol))
        return f"{formula}|{bonds}"
    except Exception as exc:
        log.debug("Skeleton identifier failed: %s", exc)
        return "Unknown"
    
def _connectivity_payload(input_data) -> tuple[str, str]:
    """
    Return atom symbols and bond edges in a stable, index-preserving form.

    Output
    ------
    atom_symbols : str
        Comma-separated atom symbols in the exact atom index order used
        by the bond list, e.g. "C,C,O,O,H,H,H,H"
    bond_edges : str
        Python-literal list of 0-based bonded pairs, e.g.
        "[(0, 1), (1, 2), (1, 3)]"

    Accepts a VeloxChem Molecule or an XYZ string.
    """
    if isinstance(input_data, str):
        mol = vlx.Molecule.read_xyz_string(input_data)
    else:
        mol = input_data

    labels = list(mol.get_labels())
    bonds = sorted(_geometry_bonds(mol))
    atom_symbols = ",".join(labels)
    bond_edges = str(bonds)
    return atom_symbols, bond_edges

def _valid_multiplicity(mol_obj, charge: int) -> int:
    """Return 1 (singlet) or 2 (doublet); warn on TM detection."""
    n_elec = mol_obj.number_of_electrons()
    mult = 1 if n_elec % 2 == 0 else 2
    if any(sym in FIRST_ROW_TM for sym in mol_obj.get_labels()):
        log.warning(
            "Transition metal detected (charge %d -> mult %d). "
            "Higher spin states may be more stable -- verify manually.",
            charge, mult,
        )
    return mult


def _mol_from_smiles(smiles: str, charge: int = 0,
                     mult: Optional[int] = None):
    """
    Build a 3-D VeloxChem Molecule from SMILES using VeloxChem's own
    read_smiles() (which handles 3D embedding internally).

    Parameters
    ----------
    smiles : str
    charge : int   -- formal charge (default 0).
    mult   : int or None -- spin multiplicity. If None, inferred from
             electron count (1 = singlet, 2 = doublet).
    """
    vx_mol = vlx.Molecule.read_smiles(smiles)
    if vx_mol is None or vx_mol.number_of_atoms() == 0:
        raise ValueError(f"VeloxChem could not parse SMILES: {smiles!r}")
    vx_mol.set_charge(charge)
    resolved_mult = mult if mult is not None else _valid_multiplicity(vx_mol, charge)
    vx_mol.set_multiplicity(resolved_mult)
    return vx_mol


def _save_to_h5(path: str, result: GibbsResult) -> None:
    """Persist a GibbsResult to HDF5.  Warns if overwriting."""
    if os.path.exists(path):
        log.warning("Overwriting HDF5 file: %s", path)
    with h5py.File(path, "w") as f:
        for key, val in result.as_dict().items():
            if val is None:
                continue
            try:
                f.create_dataset(key, data=val)
            except TypeError:
                f.create_dataset(key, data=str(val))


def _run_sp(mol_obj, basis_name: str, xcfun: str, ri_jk: bool,
            dispersion: bool, solvation: Optional[str],
            mult: int, tmp_prefix: str) -> dict:
    """Run one DFT single-point and return the SCF result dictionary."""
    drv = vlx.ScfUnrestrictedDriver() if mult > 1 else vlx.ScfRestrictedDriver()
    drv.xcfun           = xcfun
    drv.conv_thresh     = 1e-6
    drv.ri_jk           = ri_jk
    drv.max_iter        = 200
    drv.dispersion      = dispersion
    drv.solvation_model = solvation   # None -> vacuum
    if mult > 1:
        drv.level_shifting = 0.3
    drv.filename = tmp_prefix
    drv.ostream.mute()
    basis = vlx.MolecularBasis.read(mol_obj, basis_name)
    return drv.compute(mol_obj, basis)


# ---------------------------------------------------------------------------
# Harmonic thermochemistry (self-contained, no VibrationalAnalysis needed)
# ---------------------------------------------------------------------------

def _thermal_gibbs(
    mol_obj,
    frequencies_cm: np.ndarray,
    temperature: float = 298.15,
    pressure: float = 101325.0,
) -> tuple:
    """
    Full ideal-gas harmonic-oscillator / rigid-rotor Gibbs free energy
    in Hartree, including translational, rotational, and vibrational
    contributions.

    Returns
    -------
    tuple : (G_vib_only, G_trans + G_rot + G_vib)
        G_vib_only    -- purely vibrational contribution (Hartree).
                        Use this for pKa and potentials so that
                        translational/rotational terms cancel correctly
                        between states.
        G_full        -- complete thermal correction including trans + rot.
                        Stored for reference but not used in thermochemistry.

    Why vib-only for thermochemistry
    ---------------------------------
    For isoatomic reactions (ET: M -> M+, M -> M-) the trans/rot terms are
    nearly identical for both states and cancel in the energy difference.
    For proton-transfer reactions (pKa: M -> M_dep- + H+) the proton's
    trans/rot contribution is already absorbed into G_PROTON_AQ = -0.43676 au
    (Tissandier 1998), so adding trans/rot again from the MM correction
    double-counts it and gives a wrong pKa.  Using only the vibrational
    correction avoids both problems.
        ------------------------------------------------
    Modes below 10 cm^-1 (translations, rotations, numerical noise) and
    imaginary modes (negative frequencies) are skipped.

    Parameters
    ----------
    mol_obj        : VeloxChem Molecule -- provides masses and coordinates.
    frequencies_cm : ndarray -- MM frequencies in cm^-1 (from MMHessianDriver).
    temperature    : float   -- temperature in Kelvin (default 298.15 K).
    pressure       : float   -- pressure in Pa (default 101325 Pa = 1 atm).

    Returns
    -------
    float : G in Hartree  (to be subtracted from E_MM to give G_corr).
    """
    # ---- physical constants ------------------------------------------
    h      = 6.62607015e-34     # J*s
    k_B    = 1.380649e-23      # J/K
    N_A    = 6.02214076e23     # mol^-1
    R      = k_B * N_A         # J/mol/K
    c      = 2.99792458e10     # cm/s
    pi     = np.pi
    amu    = 1.66053906660e-27  # kg
    J_to_Ha     = 1.0 / 4.3597447222e-18   # J/molecule -> Hartree
    JpMol_to_Ha = J_to_Ha / N_A            # J/mol      -> Hartree
    ang_to_m = 1e-10

    kT = k_B * temperature     # J/molecule (for partition function arguments)
    RT = R  * temperature      # J/mol      (for molar H and G)

    # ---- molecular mass ----------------------------------------------
    # Use VeloxChem's masses (amu), convert to kg
    masses_amu = np.array(mol_obj.get_masses(), dtype=float)  # (N,) in amu
    M_kg       = float(np.sum(masses_amu)) * amu  # total mass in kg

    # ---- moments of inertia ------------------------------------------
    coords_ang = mol_obj.get_coordinates_in_angstrom()   # (N, 3) A
    coords_m   = coords_ang * ang_to_m
    masses_kg  = masses_amu * amu

    # Centre of mass
    com = np.sum(masses_kg[:, None] * coords_m, axis=0) / M_kg
    r   = coords_m - com   # (N, 3) relative coords in m

    # Inertia tensor
    I = np.zeros((3, 3))
    for i, (m, ri) in enumerate(zip(masses_kg, r)):
        I[0, 0] += m * (ri[1]**2 + ri[2]**2)
        I[1, 1] += m * (ri[0]**2 + ri[2]**2)
        I[2, 2] += m * (ri[0]**2 + ri[1]**2)
        I[0, 1] -= m * ri[0] * ri[1]
        I[0, 2] -= m * ri[0] * ri[2]
        I[1, 2] -= m * ri[1] * ri[2]
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    I_eigs = np.sort(np.linalg.eigvalsh(I))   # (3,) ascending

    # Detect linear molecule: smallest moment negligible
    is_linear = (I_eigs[0] < 1e-44)   # ~1e-44 kg*m2 threshold
    sigma = 1   # symmetry number (conservative default)

    # ---- translational partition function (Sackur-Tetrode) -----------
    # q_trans = (2pi M kT / h2)^(3/2) * kT/P   (per molecule)
    # Molar entropy:  S_trans = R * (ln(q_trans) + 5/2)
    # Molar enthalpy: H_trans = 5/2 * RT  (3/2 kT kinetic + kT pV work)
    # Molar G_trans  = H_trans - T * S_trans
    lam_th  = h / np.sqrt(2 * pi * M_kg * kT)   # thermal de Broglie wavelength
    q_trans = (kT / pressure) / lam_th**3        # dimensionless, per molecule

    H_trans = 2.5 * RT                            # J/mol
    S_trans = R * (np.log(q_trans) + 2.5)         # J/mol/K  (Sackur-Tetrode)
    G_trans = H_trans - temperature * S_trans      # J/mol

    # ---- rotational partition function (rigid rotor) -----------------
    if mol_obj.number_of_atoms() == 1:
        G_rot = 0.0

    elif is_linear:
        I_rot = I_eigs[1]
        if I_rot < 1e-50:
            G_rot = 0.0
        else:
            q_rot = (8 * pi**2 * I_rot * kT) / (sigma * h**2)
            H_rot = RT                             # J/mol
            S_rot = R * (np.log(q_rot) + 1.0)     # J/mol/K
            G_rot = H_rot - temperature * S_rot    # J/mol

    else:
        I_A, I_B, I_C = I_eigs
        if I_A < 1e-50 or I_B < 1e-50 or I_C < 1e-50:
            G_rot = 0.0
        else:
            q_rot = (pi**0.5 / sigma) * np.sqrt(
                (8 * pi**2 * I_A * kT / h**2)
                * (8 * pi**2 * I_B * kT / h**2)
                * (8 * pi**2 * I_C * kT / h**2)
            )
            H_rot = 1.5 * RT                       # J/mol
            S_rot = R * (np.log(q_rot) + 1.5)      # J/mol/K
            G_rot = H_rot - temperature * S_rot     # J/mol

    # ---- vibrational partition function (quasi-RRHO, molar) ----------
    # Each mode contributes:
    #   U_i   = N_A * h*nu * (0.5 + 1/(exp(x)-1))    J/mol
    #   S_i   = R   * (x/(exp(x)-1) - ln(1-exp(-x))) J/mol/K
    #   G_i   = U_i - T * S_i                         J/mol
    # where x = h*nu / (k_B * T)
    cutoff       = 10.0   # cm^-1 -- quasi-RRHO floor
    sorted_freqs = np.sort(frequencies_cm)
    n_skip       = 6      # always skip 6 translations + rotations
    vib_freqs    = sorted_freqs[n_skip:]   # 3N-6 genuine modes

    G_vib_J = 0.0   # accumulate in J/mol
    for nu in vib_freqs:
        if nu <= 0.0:
            continue

        nu_eff = max(nu, cutoff)
        hnu    = h * nu_eff * c   # J per quantum (single molecule)

        if temperature < 1e-6:
            G_vib_J += N_A * 0.5 * hnu   # ZPE only, molar
            continue

        x = hnu / kT   # dimensionless

        if x > 700:
            G_vib_J += N_A * 0.5 * hnu
            continue

        exp_x  = np.exp(x)
        U_i    = N_A * hnu * (0.5 + 1.0 / (exp_x - 1.0))   # J/mol
        S_i    = R   * (x  / (exp_x - 1.0) - np.log(1.0 - np.exp(-x)))  # J/mol/K
        G_vib_J += U_i - temperature * S_i

    G_vib = G_vib_J * JpMol_to_Ha   # J/mol -> Hartree

    # ---- convert J/mol -> Hartree and return -------------------------
    G_trans_Ha  = G_trans * JpMol_to_Ha
    G_rot_Ha    = G_rot   * JpMol_to_Ha
    G_trans_rot = G_trans_Ha + G_rot_Ha
    # G_vib already in Hartree

    return G_vib, G_trans_rot, G_vib + G_trans_rot


# ---------------------------------------------------------------------------
# Bond perception from DFT geometry
# ---------------------------------------------------------------------------

def _geometry_bonds(mol_obj, scale: float = 1.2) -> set:
    """
    Detect bonded atom pairs using VeloxChem's own connectivity matrix.

    Returns
    -------
    set of (int, int) -- 0-based, always (min, max) ordered.
    """
    conn  = mol_obj.get_connectivity_matrix()
    n     = mol_obj.number_of_atoms()
    bonds = set()
    for i in range(n):
        for j in range(i + 1, n):
            if conn[i, j] == 1:
                bonds.add((i, j))
    return bonds


# ---------------------------------------------------------------------------
# xTB Gibbs thermal correction
# ---------------------------------------------------------------------------

def _xtb_gibbs_correction(
    mol_obj,
    temperature: float = 298.15,
) -> tuple:
    """
    Compute the Gibbs thermal correction using xTB + VeloxChem VibrationalAnalysis.

    VibrationalAnalysis.compute() returns a dict with:
      'gibbs_free_energy' -- G = E_elec(xTB) + ZPE + H_vib + H_trans + H_rot
                             - T*(S_vib + S_trans + S_rot)  (Hartree)

    The xTB electronic energy is available as vib.hessian_driver.elec_energy.

    The thermal correction added on top of the DFT single-point is:

        g_corr = gibbs_free_energy - elec_energy(xTB)
               = ZPE + H_vib + H_trans + H_rot - T*(S_vib + S_trans + S_rot)

    This is the standard composite protocol: DFT provides the electronic
    energy reference, xTB provides the thermal correction.

    Parameters
    ----------
    mol_obj     : VeloxChem Molecule at the DFT-optimised geometry.
    temperature : float -- temperature in Kelvin (default 298.15 K).

    Returns
    -------
    float : g_corr in Hartree  (gibbs_free_energy - elec_energy).
    """
    from veloxchem import XtbDriver, VibrationalAnalysis

    xtb_drv = XtbDriver()
    xtb_drv.ostream.mute()

    vib = VibrationalAnalysis(xtb_drv)
    vib.ostream.mute()
    vib.temperature = temperature

    import io
    from contextlib import redirect_stdout, redirect_stderr
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        vib_results = vib.compute(mol_obj)

    gibbs  = float(vib_results['gibbs_free_energy'])
    e_xtb  = float(vib.hessian_driver.elec_energy)
    g_corr = gibbs - e_xtb

    log.debug(
        "_xtb_gibbs_correction: G(xTB)=%.6f  E_elec(xTB)=%.6f  g_corr=%.6f au",
        gibbs, e_xtb, g_corr,
    )
    return g_corr


def _repair_ff_bonds(ff_gen, mol_obj, scale: float = 1.2) -> int:
    """
    Add any bonds present in the DFT geometry but missing from the GAFF
    force field, then rebuild the topology so angles/dihedrals are consistent.

    This is the primary defence against GAFF mis-typing radical cations and
    other unusual charge states: their bonding topology is taken as ground
    truth from the DFT-optimised geometry rather than from GAFF's atom-type
    perception.

    Force constant scaling
    ----------------------
    Missing bonds are added with a force constant scaled from a reference
    C-C aromatic value by the inverse square of the actual bond length::

        k = k_ref * (r_ref / r_actual)2

    where k_ref = 200 000 kJ/mol/nm2 and r_ref = 0.140 nm.  This gives
    a stiff but not unreasonable constraint that keeps the Hessian positive
    definite without requiring an accurate GAFF type assignment.

    Parameters
    ----------
    ff_gen  : MMForceFieldGenerator -- modified in place.
    mol_obj : VeloxChem Molecule -- the optimised geometry.
    scale   : float -- covalent radius tolerance (default 1.2).

    Returns
    -------
    int : number of bonds added.
    """
    geom_bonds = _geometry_bonds(mol_obj, scale)

    ff_bonds = set()
    for i, j in ff_gen.bonds:
        ff_bonds.add((min(i, j), max(i, j)))

    missing = geom_bonds - ff_bonds
    if not missing:
        return 0

    coords = mol_obj.get_coordinates_in_angstrom()
    k_ref  = 200_000.0   # kJ/mol/nm2
    r_ref  = 0.140       # nm

    for i, j in sorted(missing):
        r_nm = float(np.linalg.norm(coords[i] - coords[j])) * 0.1
        r_nm = max(r_nm, 0.05)
        log.info(
            "  Bond repair: adding bond %d-%d  r=%.3f nm  "
            "k=%.0f kJ/mol/nm2 [geometry-based]",
            i + 1, j + 1, r_nm, k_ref * (r_ref / r_nm) ** 2,
        )
        # add_bond updates connectivity_matrix in place (expects 1-based)
        ff_gen.add_bond([i + 1, j + 1])

    n_added = len(missing)

    # Rebuild topology so angles, dihedrals, and impropers are all consistent
    # with the corrected connectivity.
    ff_gen.create_topology(mol_obj, resp=False)

    # Override force constants for the repaired bonds with distance-scaled values
    for i, j in sorted(missing):
        key = (min(i, j), max(i, j))
        if key not in ff_gen.bonds:
            key = (max(i, j), min(i, j))
        if key in ff_gen.bonds:
            r_nm = float(np.linalg.norm(coords[i] - coords[j])) * 0.1
            r_nm = max(r_nm, 0.05)
            ff_gen.bonds[key]["force_constant"] = k_ref * (r_ref / r_nm) ** 2
            ff_gen.bonds[key]["comment"] += " [geometry-repaired]"

    return n_added


# ---------------------------------------------------------------------------
# MM Gibbs correction  (no VibrationalAnalysis dependency)
# ---------------------------------------------------------------------------

def _mm_gibbs_correction(
    mol_obj,
    hessian_step_size: float,
    openmm_platform: str,
    temperature: float = 298.15,
) -> float:
    """
    Compute G_corr = G_thermal(MM) - E_MM  in Hartree.

    After building the GAFF topology, ``_repair_ff_bonds`` validates the
    bonding against the DFT geometry and adds any missing bonds with
    distance-scaled force constants before rebuilding the topology.  This
    prevents radical cations and other unusual charge states from having
    broken connectivity that produces unphysical Hessian eigenvalues.

    All intermediate files are written to a TemporaryDirectory and deleted
    automatically on exit.

    Parameters
    ----------
    mol_obj           : VeloxChem Molecule at the optimised geometry.
    hessian_step_size : Adaptive step ratio for MMHessianDriver.
    openmm_platform   : OpenMM platform string ("CPU", "CUDA", ...).
    temperature       : Temperature in Kelvin (default 298.15 K).

    Returns
    -------
    float : G_corr in Hartree.
    """
    import tempfile
    from veloxchem import MMForceFieldGenerator, MMHessianDriver

    with tempfile.TemporaryDirectory() as tmpdir:

        # ---- 1. Force field (EEM charges, no QM) ---------------------
        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.molecule_name = os.path.join(tmpdir, "mol")
        ff_gen.create_topology(mol_obj, resp=False)

        # ---- 2. Bond repair ------------------------------------------
        # Check that every bond in the DFT geometry is present in the GAFF
        # topology.  Missing bonds (common for radical cations) are added
        # with geometry-derived force constants and the topology is rebuilt.
        n_repaired = _repair_ff_bonds(ff_gen, mol_obj)
        if n_repaired:
            log.warning(
                "  Bond repair: %d bond(s) missing from GAFF topology added "
                "from DFT geometry with scaled force constants.  "
                "Thermal correction for this state is approximate.",
                n_repaired,
            )

        # ---- 3. Numerical MM Hessian + frequencies -------------------
        # MMDriver is no longer needed for the energy -- E_MM is not used
        # in the composite protocol (see comment at step 5 below).
        # MMHessianDriver builds its own OpenMM simulation from ff_gen.
        hess_drv = MMHessianDriver()
        hess_drv.ostream.mute()
        hess_drv.step_size       = hessian_step_size
        hess_drv.openmm_platform = openmm_platform
        hess_drv.filename_prefix = os.path.join(tmpdir, "hess_res")
        hess_drv.compute(mol_obj, ff_gen)

        # Sanity check: after skipping the 6 lowest modes (translations +
        # rotations), count how many genuine vibrational modes are negative.
        # A well-behaved GAFF force field should give zero or at most 1-2
        # slightly negative modes from numerical noise in the finite differences.
        # More than 20% negative internal modes indicates a broken force field.
        n_atoms        = mol_obj.number_of_atoms()
        n_internal     = max(1, 3 * n_atoms - 6)
        sorted_freqs   = np.sort(hess_drv.frequencies)
        internal_freqs = sorted_freqs[6:]   # skip 6 rigid-body modes
        n_negative     = int(np.sum(internal_freqs < -10.0))   # < -10 cm^-1 is clearly unphysical
        if n_negative > max(1, n_internal // 5):
            raise ValueError(
                f"MM Hessian has {n_negative}/{n_internal} negative internal "
                f"frequencies (< -10 cm^-1) -- GAFF force field cannot describe "
                f"charge={mol_obj.get_charge()}, mult={mol_obj.get_multiplicity()}. "
                f"Most negative: {internal_freqs[0]:.1f} cm^-1."
            )

        # ---- 5. Vibrational and rigid-body thermal corrections -------
        # We return three quantities:
        #   g_vib       -- purely vibrational (ZPE + H_vib - T*S_vib)
        #   g_trans_rot -- translational + rotational free energy
        #   g_full      -- g_vib + g_trans_rot  (complete thermal correction)
        #
        # G_total = E_SP(DFT, solvent) + g_vib + g_trans_rot
        #
        # All three are stored in the HDF5 output for inspection.
        # The DFT single-point provides the electronic energy reference;
        # E_MM is never subtracted (it plays no role in the composite protocol).
        g_vib, g_trans_rot, g_full = _thermal_gibbs(
            mol_obj, hess_drv.frequencies, temperature
        )

    log.debug(
        "_mm_gibbs_correction: G_vib=%.6f  G_trans_rot=%.6f  G_full=%.6f au",
        g_vib, g_trans_rot, g_full,
    )
    return g_vib, g_trans_rot, g_full


# ---------------------------------------------------------------------------
# Main calculator class
# ---------------------------------------------------------------------------

class RedoxCalculator:
    """
    Compute reduction potentials and pKa values via a composite QM/xTB protocol.

    Thermodynamic cycle
    -------------------
    ::

        G_total = E_SP(solvent, basis_sp, xcfun_sp)
                + G_vib(xTB) + G_trans_rot(xTB)

    The thermal correction uses xTB (GFN2-xTB) via VeloxChem's XtbDriver
    and VibrationalAnalysis.  xTB correctly assigns X-H bond parameters
    (O-H, N-H, C-H) and produces reliable vibrational frequencies for all
    common organic and organometallic molecules without force-field typing
    failures.

    Parameters
    ----------
    output_folder : str
        Directory for HDF5 output files (created if absent).
        Default: ``"redox_results"``.
    basis_opt : str
        Basis set for DFT geometry optimisation. Default: ``"def2-SVP"``.
    basis_sp : str
        Basis set for DFT single-point energies. Default: ``"def2-SVPD"``.
    xcfun_opt : str
        XC functional for geometry optimisation. Default: ``"b3lyp"``.
    xcfun_sp : str
        XC functional for single-point energies. Default: ``"b3lyp"``.
    solvation_model : str or None
        Implicit solvent model. Default: ``"smd"``.
    pH : float
        Solution pH. Default: ``7.0``.
    temperature : float
        Temperature in Kelvin. Default: ``298.15``.
    ri_jk : bool
        Use RI-JK density fitting for single-points. Default: ``True``.
    sp_dispersion : bool
        Apply D3BJ dispersion to single-points. Default: ``True``.
    """

    def __init__(
        self,
        output_folder: str        = "redox_results",
        basis_opt: str            = "def2-SVP",
        basis_sp: str             = "def2-SVPD",
        xcfun_opt: str            = "b3lyp",
        xcfun_sp: str             = "b3lyp",
        solvation_model: str      = "smd",
        pH: float                 = 7.0,
        temperature: float        = 298.15,
        ri_jk: bool               = True,
        sp_dispersion: bool       = True,
    ) -> None:
        self.output_folder     = output_folder
        self.basis_opt         = basis_opt
        self.basis_sp          = basis_sp
        self.xcfun_opt         = xcfun_opt
        self.xcfun_sp          = xcfun_sp
        self.solvation_model   = solvation_model
        self.pH                = pH
        self.temperature       = temperature
        self.ri_jk             = ri_jk
        self.sp_dispersion     = sp_dispersion

        # Named solvation registry: { "MD": {"M": -0.0145, "M+": -0.0089, ...}, ... }
        # Populated via register_solvation(); values are attached to GibbsResult
        # objects during compute_gibbs() so they travel with the results.
        self._solvation_registry: dict[str, dict[str, float]] = {}

        os.makedirs(self.output_folder, exist_ok=True)

    # ------------------------------------------------------------------
    # Named solvation interface
    # ------------------------------------------------------------------

    def register_solvation(
        self,
        name: str,
        energies: dict[str, float],
    ) -> None:
        """
        Register a set of named solvation free energies to be attached to
        computed results during ``run()``.

        Each entry replaces ``scf_solvent`` in the G_total formula when you
        later call ``profile.apply_solvation(name, calc)``::

            G_total(name) = scf_vacuum + energies[state_label] + g_corr

        Parameters
        ----------
        name : str
            Arbitrary label for this solvation method, e.g. ``"MD"``,
            ``"COSMO-RS"``, ``"FEP"``.  This is the key you pass to
            ``apply_solvation()`` later.
        energies : dict[str, float]
            Mapping from state label to solvation free energy in Hartree.
            State labels are the base labels used internally: ``"M"``,
            ``"M+"``, ``"M-"``, ``"MH+"``, ``"MH_rad"``, ``"M_deprot"``,
            ``"M_deprot_rad"``.  You do not need to provide all states --
            missing ones keep their default SMD-based G_total.

        Examples
        --------
        Register FEP values before running::

            calc.register_solvation("MD", {
                "M":  -0.0145,
                "M+": -0.0089,
                "M-": -0.0201,
            })
            profile = calc.run("c1ccccc1O")
            profile_md = profile.apply_solvation("MD", calc)

        Load from a JSON file::

            import json
            with open("fep_solvation.json") as f:
                calc.register_solvation("FEP", json.load(f))
        """
        self._solvation_registry[name] = dict(energies)
        log.info(
            "Registered solvation '%s' for states: %s",
            name, list(energies.keys()),
        )

    def clear_solvation(self, name: Optional[str] = None) -> None:
        """
        Remove named solvation registrations from the calculator.

        Parameters
        ----------
        name : str or None
            If given, remove only that solvation set.
            If None (default), remove all registered solvation sets.

        Note: this does not affect ``solvation_data`` already attached to
        GibbsResult objects from a previous ``run()``.
        """
        if name is None:
            self._solvation_registry.clear()
            log.info("All solvation registrations cleared.")
        else:
            self._solvation_registry.pop(name, None)
            log.info("Solvation registration '%s' cleared.", name)

    # ------------------------------------------------------------------
    # Thermochemistry formulas
    # ------------------------------------------------------------------

    def reduction_potential(
        self,
        g_ox: float,
        g_red: float,
        n: int = 1,
        m: int = 0,
    ) -> float:
        """
        pH-corrected reduction potential (V vs SHE).

        Half-reaction:  Ox + n e- + m H+ -> Red

        Parameters
        ----------
        g_ox, g_red : float  -- Gibbs energies of Ox and Red states (au).
        n : int              -- electrons transferred.
        m : int              -- protons transferred (0 = pure ET, 1 = PCET).

        Returns
        -------
        float : E  (V vs SHE)
        """
        dg_au = g_red - (g_ox + m * G_PROTON_AQ + n * G_ELECTRON)
        e_abs = -dg_au * CONV_EV / n              # vacuum-scale absolute potential
        return (e_abs - E_REF_SHE) - (NERNST_FACTOR / n) * m * self.pH

    def pka(self, g_acid: float, g_base: float) -> float:
        """
        pKa for the deprotonation  HA -> A- + H+.

        Parameters
        ----------
        g_acid : float  -- G of the protonated (HA) form (au).
        g_base : float  -- G of the conjugate base (A-) form (au).

        Returns
        -------
        float : pKa  (positive = weak acid, negative = strong acid)
        """
        dg_kcal = ((g_base + G_PROTON_AQ) - g_acid) * AU_TO_KCAL
        return dg_kcal / RT_LN10_KCAL

    # ------------------------------------------------------------------
    # Protonation / deprotonation site handling
    # ------------------------------------------------------------------

    def find_protic_sites(
        self, mol_obj
    ) -> tuple[list[int], list[int]]:
        """
        Identify acidic H atoms and basic heavy atoms (O, N, S, P).

        Returns
        -------
        acidic_H : list[int]
            Indices of H atoms bonded to a heteroatom within the bond cutoff.
        basic_atoms : list[int]
            Indices of heteroatoms that could be protonated.
        """
        labels = mol_obj.get_labels()
        coords = mol_obj.get_coordinates_in_angstrom()
        acidic_H, basic_atoms = [], []
        for i, si in enumerate(labels):
            if si in BOND_CUTOFF:
                basic_atoms.append(i)
                cutoff = BOND_CUTOFF[si]
                for j, sj in enumerate(labels):
                    if sj == "H":
                        if float(np.linalg.norm(coords[i] - coords[j])) < cutoff:
                            acidic_H.append(j)
        return list(set(acidic_H)), list(set(basic_atoms))

    def protonate(self, mol_obj, atom_idx: int, new_charge: int):
        """
        Add one H atom 1.0 A from ``atom_idx``, directed radially outward.

        Parameters
        ----------
        mol_obj    : VeloxChem Molecule.
        atom_idx   : int -- index of the basic atom to protonate.
        new_charge : int -- formal charge of the resulting species.

        Returns
        -------
        VeloxChem Molecule (protonated).
        """
        labels = mol_obj.get_labels()
        coords = mol_obj.get_coordinates_in_angstrom()
        center = np.mean(coords, axis=0)
        vec    = coords[atom_idx] - center
        norm   = np.linalg.norm(vec)
        vec    = vec / norm if norm > 1e-6 else np.array([0.0, 0.0, 1.0])
        new_coords = np.vstack([coords, coords[atom_idx] + vec * 1.0])
        new_labels = list(labels) + ["H"]
        xyz = "\n".join(
            [str(len(new_labels)), ""]
            + [f"{s:4s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}"
               for s, c in zip(new_labels, new_coords)]
        )
        new_mol = vlx.Molecule.read_xyz_string(xyz)
        new_mol.set_charge(new_charge)
        new_mol.set_multiplicity(_valid_multiplicity(new_mol, new_charge))
        return new_mol

    def deprotonate(self, mol_obj, h_idx: int, new_charge: int):
        """
        Remove hydrogen atom ``h_idx`` from the molecule.

        Parameters
        ----------
        mol_obj    : VeloxChem Molecule.
        h_idx      : int -- index of the H atom to remove.
        new_charge : int -- formal charge of the resulting species.

        Returns
        -------
        VeloxChem Molecule (deprotonated).
        """
        labels     = mol_obj.get_labels()
        coords     = mol_obj.get_coordinates_in_angstrom()
        new_labels = [labels[i] for i in range(len(labels)) if i != h_idx]
        new_coords = np.delete(coords, h_idx, axis=0)
        xyz = "\n".join(
            [str(len(new_labels)), ""]
            + [f"{s:4s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}"
               for s, c in zip(new_labels, new_coords)]
        )
        new_mol = vlx.Molecule.read_xyz_string(xyz)
        new_mol.set_charge(new_charge)
        new_mol.set_multiplicity(_valid_multiplicity(new_mol, new_charge))
        return new_mol

    # ------------------------------------------------------------------
    # Core QC: composite Gibbs energy for one state
    # ------------------------------------------------------------------

    def compute_gibbs(
        self,
        mol_obj,
        label: str = "",
    ) -> Optional[GibbsResult]:
        """
        Compute the composite Gibbs free energy for one electronic state.

        Steps
        -----
        1. DFT geometry optimisation (``basis_opt`` / ``xcfun_opt``, gas-phase).

        3. MM Gibbs thermal correction (GAFF + EEM charges + numerical Hessian).
        4. DFT single-point in implicit solvent (always computed).
        5. DFT single-point in vacuum (always computed).
        6. G_total = scf_solvent + g_corr   (default, SMD-based).

        Any named solvation free energies pre-registered via
        ``register_solvation()`` are attached to the result's
        ``solvation_data`` dict so they can be used later via
        ``GibbsResult.g_total_with(name)`` or ``RedoxProfile.apply_solvation()``.

        Parameters
        ----------
        mol_obj      : VeloxChem Molecule with charge and multiplicity set.
        label        : str -- short identifier for output filenames.

        Returns
        -------
        GibbsResult on success, None on failure.

        Notes
        -----
        Fragmentation is detected by comparing the number of connected
        components in the input geometry against the optimised geometry,
        using RDKit covalent-radius bond perception.  This correctly handles
        charge/protonation states (MH+, M_deprot, etc.) whose formula
        legitimately differs from M -- they are not flagged as fragmented

        """
        charge = int(mol_obj.get_charge())
        mult   = _valid_multiplicity(mol_obj, charge)
        mol_obj.set_multiplicity(mult)

        formula        = _mol_formula(mol_obj)
        file_name_base = f"{formula}_q{charge}_m{mult}_{label}"

        log.info("compute_gibbs start: %s", file_name_base)

        # All VeloxChem SCF scratch files (*.h5, *_scf.h5) and geomeTRIC
        # optimisation files are written to a single TemporaryDirectory so
        # they are cleaned up automatically on success or failure.
        # Only the final result HDF5 is written to output_folder.
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            return self._compute_gibbs_in_tmpdir(
                mol_obj, label, file_name_base, charge, mult, tmpdir,
            )

    def _compute_gibbs_in_tmpdir(
        self, mol_obj, label, file_name_base, charge, mult, tmpdir,
    ):
        """Inner implementation of compute_gibbs running inside a temp dir."""

        # ---- DFT geometry optimisation --------------------------------
        # ri_jk is forced off here: ScfGradientDriver (called internally by
        # OptimizationDriver) does not yet support RI-JK in VeloxChem and will
        # raise an AssertionError if it is enabled.  RI-JK is still used for
        # the single-point steps below where no gradient is needed.
        scf_opt = (
            vlx.ScfUnrestrictedDriver() if mult > 1
            else vlx.ScfRestrictedDriver()
        )
        scf_opt.xcfun       = self.xcfun_opt
        scf_opt.ri_jk       = False
        scf_opt.max_iter    = 200
        scf_opt.conv_thresh = 1e-5
        scf_opt.dispersion  = True
        scf_opt.grid_level  = 3
        if mult > 1:
            scf_opt.level_shifting = 0.3
        scf_opt.filename = os.path.join(tmpdir, "opt")
        scf_opt.ostream.mute()

        try:
            log.info("  [%s] Step 1: reading basis %s", label, self.basis_opt)
            basis_opt            = vlx.MolecularBasis.read(mol_obj, self.basis_opt)
            log.info("  [%s] Step 2: geometry optimisation", label)
            opt_drv              = vlx.OptimizationDriver(scf_opt)
            opt_drv.conv_energy  = 1e-5
            opt_drv.conv_grms    = 3e-4
            opt_drv.conv_gmax    = 1.2e-3
            opt_drv.max_iter     = 50
            opt_drv.conv_maxiter = True
            opt_drv.ostream.mute()

            import io
            from contextlib import redirect_stdout, redirect_stderr
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                opt_res = opt_drv.compute(mol_obj, basis_opt)
            final_xyz = opt_res["final_geometry"]
            log.info("  [%s] Step 2 done: geometry converged", label)

            # ---- fragmentation check ----------------------------------
            # Fragmentation = the optimised geometry has more connected
            # components than the input geometry, OR the heavy-atom graph
            # of the optimised geometry is not isomorphic to the input.
            # We use geometry-based bond perception (covalent radii) via
            # _geometry_bonds(), consistent with the bond-repair code.
            # This correctly handles protonation/deprotonation states whose
            # formula legitimately differs from M.

            def _components_and_graph(xyz_or_mol):
                """
                Return (n_components, nx.Graph) for a geometry.
                Accepts either an XYZ string or a VeloxChem Molecule.
                """
                import networkx as nx
                if isinstance(xyz_or_mol, str):
                    tmp_mol = vlx.Molecule.read_xyz_string(xyz_or_mol)
                else:
                    tmp_mol = xyz_or_mol

                bonds = _geometry_bonds(tmp_mol)
                n     = tmp_mol.number_of_atoms()
                labels_tmp = tmp_mol.get_labels()

                G = nx.Graph()
                for i, sym in enumerate(labels_tmp):
                    G.add_node(i, element=sym)
                for i, j in bonds:
                    G.add_edge(i, j)

                return nx.number_connected_components(G), G

            n_comp_in,  g_in  = _components_and_graph(mol_obj)
            n_comp_opt, g_opt = _components_and_graph(final_xyz)

            # Check 1: more fragments after optimisation
            more_fragments = n_comp_opt > n_comp_in

            # Check 2: heavy-atom graph not isomorphic to input
            # (catches bond rearrangements / isomerisations)
            def _heavy_graph(G):
                import networkx as nx
                H = nx.Graph()
                mapping = {}
                new_idx = 0
                for node, data in G.nodes(data=True):
                    if data["element"] != "H":
                        H.add_node(new_idx, element=data["element"])
                        mapping[node] = new_idx
                        new_idx += 1
                for u, v in G.edges():
                    if u in mapping and v in mapping:
                        H.add_edge(mapping[u], mapping[v])
                return H

            from networkx.algorithms import isomorphism
            nm = isomorphism.categorical_node_match("element", None)
            hg_in  = _heavy_graph(g_in)
            hg_opt = _heavy_graph(g_opt)
            rearranged = (
                len(hg_in) == len(hg_opt) and
                not isomorphism.GraphMatcher(hg_in, hg_opt, node_match=nm).is_isomorphic()
            )

            is_fragmented = more_fragments or rearranged

            if is_fragmented:
                reason = "more fragments" if more_fragments else "bond rearrangement"
                log.warning(
                    "Fragmentation/rearrangement detected for '%s' (%s): "
                    "%d component(s) in -> %d out.",
                    label, reason, n_comp_in, n_comp_opt,
                )

            opt_mol = vlx.Molecule.read_xyz_string(final_xyz)
            opt_mol.set_charge(charge)
            opt_mol.set_multiplicity(mult)

            # ---- MM Gibbs thermal correction (replaces xTB) ----------
            log.info("  [%s] Step 3: xTB Gibbs correction", label)
            try:
                g_corr = _xtb_gibbs_correction(
                    opt_mol,
                    temperature = self.temperature,
                )
                # Sanity check: thermal correction should be small relative
                # to molecule size. A broken xTB run can give wildly wrong values.
                n_atoms_check = opt_mol.number_of_atoms()
                g_corr_limit  = 0.05 * n_atoms_check
                if abs(g_corr) > g_corr_limit:
                    raise ValueError(
                        f"g_corr = {g_corr:.4f} au exceeds the plausible "
                        f"limit of +/-{g_corr_limit:.3f} au for a "
                        f"{n_atoms_check}-atom molecule."
                    )
                log.info("  [%s] Step 3 done: g_corr = %.6f au", label, g_corr)
            except Exception as mm_exc:
                log.warning(
                    "  [%s] MM Gibbs correction failed (%s). "
                    "Setting g_corr = 0.0 and continuing with E_SP only.",
                    label, mm_exc,
                )
                g_corr = 0.0
            mm_corr_ok = (g_corr != 0.0)

            # ---- DFT single-points -----------------------------------
            log.info("  [%s] Step 4: solvent single-point (%s/%s)",
                     label, self.xcfun_sp, self.basis_sp)
            sp_solv = _run_sp(
                opt_mol, self.basis_sp, self.xcfun_sp,
                self.ri_jk, self.sp_dispersion, self.solvation_model,
                mult, os.path.join(tmpdir, "sp_solv"),
            )
            log.info("  [%s] Step 4 done: E_solvent = %.6f au",
                     label, sp_solv["scf_energy"])
            log.info("  [%s] Step 5: vacuum single-point", label)
            sp_vac = _run_sp(
                opt_mol, self.basis_sp, self.xcfun_sp,
                self.ri_jk, self.sp_dispersion, None,
                mult, os.path.join(tmpdir, "sp_vac"),
            )
            log.info("  [%s] Step 5 done: E_vacuum = %.6f au",
                     label, sp_vac["scf_energy"])

            # Sanity check: solvation energy = E_solvent - E_vacuum should
            # always be negative (solvent stabilises the molecule).
            # A positive value means the SMD calculation failed silently --
            # most common for anions with diffuse basis sets.
            dg_solv = sp_solv["scf_energy"] - sp_vac["scf_energy"]
            log.info("  [%s] Solvation energy = %.4f au", label, dg_solv)
            if dg_solv > 0.001:
                log.warning(
                    "  [%s] Solvation energy is POSITIVE (%.4f au = %.1f kJ/mol) "
                    "-- the SMD single-point likely did not converge properly. "
                    "The solvent SP energy will be replaced by the vacuum SP "
                    "energy for G_total to avoid a corrupted result. "
                    "Consider re-running with a tighter SCF threshold or a "
                    "different solvation model.",
                    label, dg_solv, dg_solv * 2625.5,
                )
                # Fall back to vacuum energy -- at least physically meaningful
                sp_solv_energy = sp_vac["scf_energy"]
                scf_solvent_stored = sp_solv["scf_energy"]  # store original for inspection
            else:
                sp_solv_energy     = sp_solv["scf_energy"]
                scf_solvent_stored = sp_solv["scf_energy"]

            # G_total = E_SP(solvent) + G_vib
            #
            # G_trans_rot is stored separately for inspection but NOT added
            # to G_total used in thermochemistry.  For pKa calculations the
            # proton reference G_PROTON_AQ already contains the proton's full
            # translational free energy, so adding G_trans_rot to each state
            # would double-count the proton contribution and shift all pKa
            # values by ~18 units.  For ET reactions (M -> M+, M -> M-)
            # G_trans_rot is nearly identical for both states and would cancel,
            # but including it risks numerical noise from the mass difference
            # of one electron.  G_corr_full (vib + trans + rot) is stored in
            # the HDF5 for reference and post-processing if needed.
            g_total = sp_solv_energy + g_corr

            atom_symbols, bond_edges = _connectivity_payload(final_xyz)

            result = GibbsResult(
                g_total           = g_total,
                charge            = charge,
                multiplicity      = mult,
                scf_solvent       = scf_solvent_stored,
                scf_vacuum        = sp_vac["scf_energy"],
                g_corr            = g_corr,
                geometry_xyz      = final_xyz,
                is_fragmented     = is_fragmented,
                skeleton_smiles   = _skeleton_smiles(final_xyz),
                atom_symbols      = atom_symbols,
                bond_edges        = bond_edges,
                label             = label,
                mm_corr_available = mm_corr_ok,
            )

            # Attach any pre-registered named solvation values
            base_label = label.split("_site")[0]
            for solv_name, state_map in self._solvation_registry.items():
                if base_label in state_map:
                    result.add_solvation(solv_name, state_map[base_label])
                    log.info(
                        "Attached solvation_data['%s'] = %.6f au to state '%s'.",
                        solv_name, state_map[base_label], label,
                    )

            h5_path = os.path.join(self.output_folder, f"{file_name_base}.h5")
            _save_to_h5(h5_path, result)
            log.info("compute_gibbs done: %s  G = %.6f au", label, g_total)
            return result

        except BaseException as exc:
            import traceback
            msg = traceback.format_exc()
            log.error("compute_gibbs failed for '%s': %s", label, exc)
            log.error("Full traceback:\n%s", msg)
            print(f"\n[ERROR] compute_gibbs failed for '{label}': {exc}\n{msg}",
                  flush=True)
            if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                raise
            return None

    # ------------------------------------------------------------------
    # Full thermodynamic profile
    # ------------------------------------------------------------------

    def run(self, input_data) -> RedoxProfile:
        """
        Compute the complete redox / pKa profile for a molecule.

        Parameters
        ----------
        input_data : str or VeloxChem Molecule
            A SMILES string or a pre-built VeloxChem Molecule.

        Returns
        -------
        RedoxProfile
            Dataclass with all GibbsResult objects (``profile.M``,
            ``profile.M_ox``, ``profile.MH_plus``, ...) and derived
            scalar quantities (``profile.E_red``, ``profile.pKa_M``, ...).
        """
        if isinstance(input_data, str):
            mol_obj   = _mol_from_smiles(input_data)
            mol_label = input_data
        else:
            mol_obj   = input_data
            mol_label = (
                f"{_mol_formula(mol_obj)} "
                f"(q={mol_obj.get_charge()}, m={mol_obj.get_multiplicity()})"
            )

        init_charge = int(mol_obj.get_charge())
        init_xyz    = _xyz_block(mol_obj)
        profile     = RedoxProfile()

        # ---- neutral, oxidised, reduced --------------------------------
        for attr, lbl, dq in [
            ("M",     "M",  0),
            ("M_ox",  "M+", 1),
            ("M_red", "M-", -1),
        ]:
            m_tmp = vlx.Molecule.read_xyz_string(init_xyz)
            m_tmp.set_charge(init_charge + dq)
            m_tmp.set_multiplicity(
                _valid_multiplicity(m_tmp, init_charge + dq)
            )
            setattr(
                profile, attr,
                self.compute_gibbs(m_tmp, lbl)
            )

        acidic_H, basic_atoms = self.find_protic_sites(mol_obj)

        # ---- protonation scan (basic sites) ----------------------------
        if basic_atoms:
            log.info(
                "Scanning %d basic site(s) for protonation.", len(basic_atoms)
            )
            best_plus = best_rad = None
            min_gp = min_gr = float("inf")

            for idx in basic_atoms:
                # MH+: protonated, charge +1
                res_p = self.compute_gibbs(
                    self.protonate(mol_obj, idx, init_charge + 1),
                    f"MH+_site{idx}",
                )
                if res_p and res_p.g_total < min_gp:
                    min_gp, best_plus = res_p.g_total, res_p

                res_r = self.compute_gibbs(
                    self.protonate(mol_obj, idx, init_charge),
                    f"MH_rad_site{idx}",
                )
                if res_r and res_r.g_total < min_gr:
                    min_gr, best_rad = res_r.g_total, res_r

            if best_plus:
                profile.MH_plus = best_plus
            if best_rad:
                profile.MH_rad = best_rad

        # ---- deprotonation scan (acidic sites) -------------------------
        if acidic_H:
            log.info(
                "Scanning %d acidic site(s) for deprotonation.", len(acidic_H)
            )
            best_dep = best_dep_rad = None
            min_gd = min_gdr = float("inf")

            for h_idx in acidic_H:
                # M_dep-: deprotonated anion, charge -1
                res_d = self.compute_gibbs(
                    self.deprotonate(mol_obj, h_idx, init_charge - 1),
                    f"M_deprot_site{h_idx}",
                )
                if res_d and res_d.g_total < min_gd:
                    min_gd, best_dep = res_d.g_total, res_d

                res_dr = self.compute_gibbs(
                    self.deprotonate(mol_obj, h_idx, init_charge),
                    f"M_deprot_rad_site{h_idx}",
                )
                if res_dr and res_dr.g_total < min_gdr:
                    min_gdr, best_dep_rad = res_dr.g_total, res_dr

            if best_dep:
                profile.M_deprot = best_dep
            if best_dep_rad:
                profile.M_deprot_rad = best_dep_rad

        # ---- thermal correction normalisation --------------------------
        # If one state in a pair has a MM thermal correction and the other
        # doesn't (mm_corr_available=False), the computed potential or pKa
        # would use asymmetric energies.  To ensure fair comparison, we
        # subtract g_corr from any state whose partner lacks it, effectively
        # putting both on the pure E_SP(solvent) level.
        # The original GibbsResult objects are not modified -- we build a
        # normalised g_total lookup used only for the thermochemistry below.

        def _g(res: Optional[GibbsResult]) -> Optional[float]:
            """Raw g_total for a result, or None."""
            return res.g_total if res is not None else None

        def _normalise_pair(
            res_a: Optional[GibbsResult],
            res_b: Optional[GibbsResult],
        ) -> tuple:
            """
            Return (g_a, g_b) normalised to the same thermal level.

            If one has a thermal correction and the other doesn't, strip the
            correction from the one that has it.  Returns (None, None) if
            either result is missing.
            """
            if res_a is None or res_b is None:
                return None, None
            g_a = res_a.g_total
            g_b = res_b.g_total
            stripped = []
            if res_a.mm_corr_available and not res_b.mm_corr_available:
                g_a -= res_a.g_corr
                stripped.append(res_a.label)
            elif res_b.mm_corr_available and not res_a.mm_corr_available:
                g_b -= res_b.g_corr
                stripped.append(res_b.label)
            if stripped:
                log.warning(
                    "Thermal normalisation: stripped g_corr from %s so that "
                    "both states are at the E_SP(solvent) level for a "
                    "consistent comparison.  This pair is flagged SP-ONLY.",
                    ", ".join(stripped),
                )
            return g_a, g_b

        # ---- reduction potentials --------------------------------------
        if profile.M:
            if profile.M_ox:
                g_ox, g_M = _normalise_pair(profile.M_ox, profile.M)
                if g_ox is not None:
                    profile.E_ox = self.reduction_potential(
                        g_ox=g_ox, g_red=g_M
                    )
            if profile.M_red:
                g_M, g_red = _normalise_pair(profile.M, profile.M_red)
                if g_M is not None:
                    profile.E_red = self.reduction_potential(
                        g_ox=g_M, g_red=g_red
                    )
            if profile.MH_rad:
                g_M, g_mhr = _normalise_pair(profile.M, profile.MH_rad)
                if g_M is not None:
                    profile.E_pcet_red = self.reduction_potential(
                        g_ox=g_M, g_red=g_mhr, n=1, m=1
                    )
            if profile.M_deprot_rad:
                g_mdr, g_M = _normalise_pair(
                    profile.M_deprot_rad, profile.M
                )
                if g_mdr is not None:
                    profile.E_pcet_ox = self.reduction_potential(
                        g_ox=g_mdr, g_red=g_M, n=1, m=1
                    )

        # ---- pKa values ------------------------------------------------
        if profile.MH_plus and profile.M:
            g_mhp, g_M = _normalise_pair(profile.MH_plus, profile.M)
            if g_mhp is not None:
                profile.pKa_MH_plus = self.pka(g_acid=g_mhp, g_base=g_M)

        if profile.MH_rad and profile.M_red:
            g_mhr, g_red = _normalise_pair(profile.MH_rad, profile.M_red)
            if g_mhr is not None:
                profile.pKa_MH_rad = self.pka(g_acid=g_mhr, g_base=g_red)

        if profile.M and profile.M_deprot:
            g_M, g_dep = _normalise_pair(profile.M, profile.M_deprot)
            if g_M is not None:
                profile.pKa_M = self.pka(g_acid=g_M, g_base=g_dep)

        if profile.M_ox and profile.M_deprot_rad:
            g_ox, g_mdr = _normalise_pair(
                profile.M_ox, profile.M_deprot_rad
            )
            if g_ox is not None:
                profile.pKa_M_plus = self.pka(g_acid=g_ox, g_base=g_mdr)

        self._print_summary(profile, label=mol_label)
        return profile

    # ------------------------------------------------------------------
    # Summary output
    # ------------------------------------------------------------------

    def _print_summary(self, profile: RedoxProfile, label: str = "") -> None:
        sep = "=" * 52
        print(f"\n{sep}")
        if label:
            print(f"  MOLECULE : {label}")
        print(f"  THERMODYNAMIC SUMMARY  (pH {self.pH:.1f})")
        print(sep)

        # Build a lookup of which states have quality issues so we can
        # annotate individual values inline rather than just at the bottom.
        def _flags(state_attrs: tuple) -> str:
            """Return a compact warning tag if any involved state is suspect."""
            tags = []
            for attr in state_attrs:
                res: Optional[GibbsResult] = getattr(profile, attr)
                if res is None:
                    continue
                if res.is_fragmented:
                    tags.append(f"FRAG:{attr}")
            missing = [
                attr for attr in state_attrs
                if getattr(profile, attr) is not None
                and not getattr(profile, attr).mm_corr_available
            ]
            if missing:
                states_str = ", ".join(missing)
                tags.append(f"SP-ONLY: Hessian failed for {states_str}")
            return f"  [!] {', '.join(tags)}" if tags else ""

        # ---- pKa values ------------------------------------------------
        # Map each pKa to the two states whose energies enter the formula
        pka_entries = [
            ("pKa_MH_plus", "MH+  ->  M   + H+",  ("MH_plus",  "M")),
            ("pKa_MH_rad",  "MH   ->  M-  + H+",  ("MH_rad",   "M_red")),
            ("pKa_M",       "M    ->  M_dep- + H+",("M",        "M_deprot")),
            ("pKa_M_plus",  "M+   ->  M_dep  + H+",("M_ox",     "M_deprot_rad")),
        ]
        any_pka = False
        for attr, lbl, states in pka_entries:
            val = getattr(profile, attr)
            if val is not None:
                if not any_pka:
                    print("\n  pKa values:")
                    any_pka = True
                flags = _flags(states)
                print(f"    {lbl:38s}  {val:7.2f}{flags}")

        # ---- reduction potentials --------------------------------------
        pot_entries = [
            ("E_ox",       "M  ->  M+         (ET oxidation)",   ("M_ox",        "M")),
            ("E_red",      "M  ->  M-         (ET reduction)",   ("M",           "M_red")),
            ("E_pcet_red", "M  ->  MH         (PCET reduction)", ("M",           "MH_rad")),
            ("E_pcet_ox",  "M  ->  M_dep      (PCET oxidation)", ("M_deprot_rad","M")),
        ]
        any_pot = False
        for attr, lbl, states in pot_entries:
            val = getattr(profile, attr)
            if val is not None:
                if not any_pot:
                    print("\n  Reduction potentials (V vs SHE):")
                    any_pot = True
                flags = _flags(states)
                print(f"    {lbl:42s}  {val:+8.3f} V{flags}")

        # ---- per-state warnings ----------------------------------------
        print("")
        for attr in (
            "M", "M_ox", "M_red", "MH_plus", "MH_rad",
            "M_deprot", "M_deprot_rad",
        ):
            res: Optional[GibbsResult] = getattr(profile, attr)
            if res and res.is_fragmented:
                print(
                    f"  [WARNING] {attr}: molecule fragmented during "
                    "optimisation -- energies are unreliable."
                )
            if res and not res.mm_corr_available:
                print(
                    f"  [WARNING] {attr}: MM Gibbs correction failed -- "
                    "G_total = E_SP(solvent) only, no thermal correction."
                )

        print(sep + "\n")


# ---------------------------------------------------------------------------
# pKa-only convenience subclass
# ---------------------------------------------------------------------------

class pKaCalculator(RedoxCalculator):
    """
    Convenience subclass that skips the ET charge states (M+, M-) and
    evaluates only pKa-relevant states (M, MH+, M_dep-).

    Use when only acid-base properties are needed; saves approximately half
    the DFT calculations compared to the full RedoxCalculator.

    All constructor parameters are inherited from RedoxCalculator.

    Example
    -------
    ::

        calc = pKaCalculator(xcfun_opt="pbe0", xcfun_sp="m06-2x",
                             sp_dispersion=False)
        profile = calc.run("CC(=O)O")
        print(f"pKa = {profile.pKa_M:.2f}")
    """

    def run(self, input_data) -> RedoxProfile:   # type: ignore[override]
        """Run only pKa-relevant states: M, MH+, M_dep-."""
        if isinstance(input_data, str):
            mol_obj = _mol_from_smiles(input_data)
        else:
            mol_obj = input_data

        init_charge = int(mol_obj.get_charge())
        profile     = RedoxProfile()

        profile.M = self.compute_gibbs(mol_obj, "M")

        acidic_H, basic_atoms = self.find_protic_sites(mol_obj)

        if basic_atoms and profile.M:
            best_plus, min_g = None, float("inf")
            for idx in basic_atoms:
                res = self.compute_gibbs(
                    self.protonate(mol_obj, idx, init_charge + 1),
                    f"MH+_site{idx}",
                )
                if res and res.g_total < min_g:
                    min_g, best_plus = res.g_total, res
            if best_plus:
                profile.MH_plus = best_plus
                profile.pKa_MH_plus = self.pka(
                    g_acid=profile.MH_plus.g_total,
                    g_base=profile.M.g_total,
                )

        if acidic_H and profile.M:
            best_dep, min_g = None, float("inf")
            for h_idx in acidic_H:
                res = self.compute_gibbs(
                    self.deprotonate(mol_obj, h_idx, init_charge - 1),
                    f"M_deprot_site{h_idx}",
                )
                if res and res.g_total < min_g:
                    min_g, best_dep = res.g_total, res
            if best_dep:
                profile.M_deprot = best_dep
                profile.pKa_M = self.pka(
                    g_acid=profile.M.g_total,
                    g_base=profile.M_deprot.g_total,
                )

        self._print_summary(profile)
        return profile


# ---------------------------------------------------------------------------
# CSV batch processor with graph-isomorphism deduplication
# ---------------------------------------------------------------------------

def _mol_to_graph(mol_obj) -> nx.Graph:
    """
    Build a heavy-atom NetworkX graph from a VeloxChem Molecule using
    geometry-based bond perception (_geometry_bonds / covalent radii).
    """
    labels = mol_obj.get_labels()
    bonds  = _geometry_bonds(mol_obj)

    G     = nx.Graph()
    heavy = {}   # atom index -> graph node index (heavy atoms only)
    idx   = 0
    for i, sym in enumerate(labels):
        if sym != "H":
            G.add_node(idx, element=sym)
            heavy[i] = idx
            idx += 1
    for i, j in bonds:
        if i in heavy and j in heavy:
            G.add_edge(heavy[i], heavy[j])
    return G


def _smiles_to_graph(smiles: str) -> Optional[nx.Graph]:
    """
    Build a heavy-atom NetworkX graph from a SMILES string.
    Uses _mol_from_smiles (VeloxChem first, RDKit fallback) for 3D structure,
    then _geometry_bonds for connectivity.
    """
    try:
        mol = _mol_from_smiles(smiles, charge=0)
        return _mol_to_graph(mol)
    except Exception as exc:
        log.debug("_smiles_to_graph failed for %s: %s", smiles, exc)
        return None


def _load_completed_graphs(output_folder: str) -> list:
    """
    Scan ``output_folder`` for completed HDF5 files (_M.h5) and return a
    list of heavy-atom molecular graphs for deduplication.

    Builds graphs from the stored geometry_xyz using VeloxChem +
    _geometry_bonds -- no SMILES parsing required.
    """
    graphs: list = []
    if not os.path.exists(output_folder):
        return graphs
    log.info("Building graph database from %s ...", output_folder)
    for fname in os.listdir(output_folder):
        if not fname.endswith("_M.h5"):
            continue
        try:
            with h5py.File(os.path.join(output_folder, fname), "r") as f:
                if "geometry_xyz" not in f:
                    continue
                xyz = f["geometry_xyz"][()].decode("utf-8")
            mol = vlx.Molecule.read_xyz_string(xyz)
            G   = _mol_to_graph(mol)
            graphs.append(G)
        except Exception as exc:
            log.debug("Could not read %s: %s", fname, exc)
    log.info("Graph database: %d unique scaffold(s).", len(graphs))
    return graphs


def _is_isomorphic(query_g: nx.Graph, existing: list) -> bool:
    """True if ``query_g`` is graph-isomorphic to any member of ``existing``."""
    nm = isomorphism.categorical_node_match("element", None)
    for target in existing:
        if len(target) != len(query_g):
            continue
        if isomorphism.GraphMatcher(target, query_g, node_match=nm).is_isomorphic():
            return True
    return False


def process_csv(df: pd.DataFrame, calc: RedoxCalculator) -> None:
    """
    Process a DataFrame of molecules using graph-isomorphism deduplication.

    Expected columns
    ----------------
    smiles       : SMILES string (required).
    charge       : integer formal charge (optional, default 0).
    multiplicity : spin multiplicity of the neutral reference state
                   (optional, default inferred from electron count).
                   Set this explicitly for radical ground states, e.g.
                   doublet radical cations or triplet biradicals.
    force_run    : "true"/"1"/"yes" to recompute known scaffolds (optional).
    """
    existing = _load_completed_graphs(calc.output_folder)

    for _, row in df.iterrows():
        smiles = str(row.get("smiles", "")).strip()
        if not smiles:
            continue
        charge    = int(row.get("charge", 0))
        # Read multiplicity from CSV if present; None means auto-infer
        mult_raw  = row.get("multiplicity", None)
        mult      = int(mult_raw) if pd.notna(mult_raw) else None
        force_run = (
            str(row.get("force_run", "False")).lower() in ("true", "1", "yes")
        )

        reason = "[FORCE]" if force_run else "[NEW]"
        try:
            mol = _mol_from_smiles(smiles, charge=charge, mult=mult)
        except Exception as exc:
            log.warning("Could not build molecule for %s: %s", smiles, exc)
            continue

        query_g = _mol_to_graph(mol)

        if _is_isomorphic(query_g, existing) and not force_run:
            log.info("Skip (graph exists): %.40s", smiles)
            continue

        log.info("%s Processing: %s (charge=%d, mult=%s)",
                 reason, smiles, charge, mult if mult is not None else "auto")
        try:
            calc.run(mol)
            existing.append(query_g)
        except Exception as exc:
            log.error("Failed for %s: %s", smiles, exc, exc_info=True)


def monitor_csv(
    file_path: str,
    calc: RedoxCalculator,
    interval_seconds: int = 30,
    max_iterations: Optional[int] = None,
) -> None:
    """
    Poll a CSV file for changes and trigger batch processing when it changes.

    Parameters
    ----------
    file_path        : str -- path to the CSV file.
    calc             : RedoxCalculator instance to use.
    interval_seconds : int -- polling interval (default 30 s).
    max_iterations   : int or None -- stop after N polls (None = run forever).

    CSV format: must have a ``smiles`` column; optional ``charge`` (int,
    default 0) and ``force_run`` (bool, default False) columns.
    """
    last_hash: Optional[str] = None
    iteration = 0
    print(f"Monitoring {file_path} every {interval_seconds}s. Ctrl-C to stop.")
    try:
        while max_iterations is None or iteration < max_iterations:
            if os.path.exists(file_path):
                try:
                    with open(file_path, "rb") as fh:
                        current_hash = hashlib.md5(fh.read()).hexdigest()
                    if current_hash != last_hash:
                        print(f"\n[!] Change detected at {time.ctime()}")
                        process_csv(pd.read_csv(file_path), calc)
                        last_hash = current_hash
                    else:
                        print(".", end="", flush=True)
                except Exception as exc:
                    log.error("Error reading CSV: %s", exc)
            else:
                log.warning("File not found: %s", file_path)
            time.sleep(interval_seconds)
            iteration += 1
    except KeyboardInterrupt:
        print("\nMonitoring stopped.")


# ---------------------------------------------------------------------------
# Demo entry point
# ---------------------------------------------------------------------------

def _demo() -> RedoxProfile:
    """
    Quick demonstration on phenol (PBE0/m06-2x).
    Run directly:  python redox_calculator.py
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    calc = RedoxCalculator(
        output_folder = "demo_results",
        xcfun_opt     = "pbe0",
        xcfun_sp      = "m06-2x",
        sp_dispersion = False,
        pH            = 7.0,
    )
    print("Running full redox/pKa profile for phenol ...")
    profile = calc.run("c1ccccc1O")

    print("\nRunning pKa-only profile for acetic acid ...")
    pka_calc = pKaCalculator(
        output_folder = "demo_results",
        xcfun_opt     = "pbe0",
        xcfun_sp      = "m06-2x",
        sp_dispersion = False,
    )
    pka_calc.run("CC(=O)O")

    return profile


if __name__ == "__main__":
    _demo()
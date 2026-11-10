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
"""Metal sites embedded in molecular or framework environments.

``FrameworkSystem`` prepares a VeloxChem baseline force field once for a
cluster or periodic cell, then builds one independent OpenMM system per
candidate. VeloxChem owns every non-target parameter: bonded terms, ligand and
hydroxide charges, and LJ parameters. Only the selected target metal site
changes. Metal-ligand contacts are nonbonded: the metal bonds, angles, and
torsions that VeloxChem infers from distances are removed before the OpenMM
system is created, so exclusions follow the remaining covalent graph.

Structural properties (``distance``, ``angle``, ``coordination_rmsd``,
``coordination_number``) are measured on a minimized geometry or averaged over
an MD trajectory, using minimum-image vectors in periodic cells (``geometry.py``).

The baseline section comes first: ``FrameworkBaseline`` is an immutable
snapshot (serialized OpenMM system, positions, charges, covalent graph) that
every candidate deserializes afresh; ``strip_metal_terms`` removes the
distance-inferred metal terms; ``compare_to_baseline`` lists anything a
candidate changed outside its site.
"""

from __future__ import annotations

from copy import copy, deepcopy
from dataclasses import dataclass, field
import xml.etree.ElementTree as ET
import hashlib
import inspect
import io
from pathlib import Path
import tempfile
import time
from typing import Any

import numpy as np

from ..atomtypeidentifier import AtomTypeIdentifier
from ..mmforcefieldgenerator import MMForceFieldGenerator
from ..errorhandler import assert_msg_critical
from ..molecule import Molecule
from ..outputstream import OutputStream
from ..veloxchemlib import bohr_in_angstrom
from .generator import new_c4_force, set_c4_pairs
from .geometry import (minimum_image, align_template, angle_degrees,
                       atom_labels, box_widths, count_within, distance_angstrom,
                       element_matched_rmsd, kabsch_rmsd, nearest_heavy_atoms,
                       resolve_atom)
from .models import all_finite, CalculationResult, default_ostream, ANGSTROM_TO_NM

# ==========================================================================
# framework baseline
# ==========================================================================


def metal_indices(molecule) -> tuple[int, ...]:
    """Indices of the atoms VeloxChem's atom typing treats as metals."""
    identifier = AtomTypeIdentifier(ostream=OutputStream(None))
    return tuple(i for i, elem_id in enumerate(molecule.get_element_ids())
                 if identifier.element_id_is_metal(elem_id))


WATER_MODELS_FOR_GROUPS = ("spc", "spce", "tip3p")


def covalent_groups(molecule) -> list[list[int]]:
    """Covalent fragments once metal-ligand bonds are removed.

    Uses VeloxChem's connectivity. Metals are left out; every other atom is in
    exactly one group (a ligand, a water, a counter-ion, ...), sorted by its
    first atom.
    """
    connectivity = np.asarray(molecule.get_connectivity_matrix(), dtype=bool)
    metals = set(metal_indices(molecule))
    unassigned = [i for i in range(molecule.number_of_atoms()) if i not in metals]
    seen, groups = set(), []
    for start in unassigned:
        if start in seen:
            continue
        group, stack = [], [start]
        seen.add(start)
        while stack:
            atom = stack.pop()
            group.append(atom)
            for other in np.flatnonzero(connectivity[atom]):
                if int(other) not in metals and int(other) not in seen:
                    seen.add(int(other))
                    stack.append(int(other))
        groups.append(sorted(group))
    return groups


_FORCE_FIELD_DATA = ("atoms", "bonds", "angles", "dihedrals", "impropers",
                     "pairs", "partial_charges", "atom_types",
                     "unique_atom_types", "connectivity_matrix")


@dataclass(frozen=True)
class FrameworkBaseline:
    """
    Immutable baseline snapshot; every candidate deserializes a fresh copy.

    Charges are in e. ``bonds`` is the covalent graph after metal-ligand bonds
    were removed.

    :param system_xml:
        The serialized OpenMM system.
    :param positions_nm:
        The read-only positions in nm.
    :param elements:
        The element labels.
    :param target_index:
        The index of the target metal.
    :param metal_indices:
        The indices of all metals.
    :param charges:
        The particle charges in e.
    :param bonds:
        The covalent bonds without metal-ligand bonds.
    :param removed_metal_terms:
        The number of removed metal terms by kind.
    :param box_vectors_nm:
        The periodic box vectors in nm, or None.
    :param source:
        The description of the force-field source.
    :param fingerprint:
        The fingerprint of system, positions and box.
    """

    system_xml: str
    positions_nm: np.ndarray
    elements: tuple[str, ...]
    target_index: int
    metal_indices: tuple[int, ...]
    charges: tuple[float, ...]
    bonds: tuple[tuple[int, int], ...]
    removed_metal_terms: dict[str, int]
    box_vectors_nm: np.ndarray | None
    source: str
    fingerprint: str

    @property
    def periodic(self) -> bool:
        """
        Checks whether the baseline has a periodic cell.

        :return:
            True for a periodic baseline.
        """

        return self.box_vectors_nm is not None

    @property
    def total_charge(self) -> float:
        """
        Gets the sum of the baseline particle charges.

        :return:
            The total charge in e.
        """

        return float(sum(self.charges))

    def new_system(self):
        """
        Creates an independent OpenMM System; the snapshot itself never changes.

        :return:
            The OpenMM System.
        """

        import openmm as mm
        return mm.XmlSerializer.deserialize(self.system_xml)


def copy_force_field(ff):
    """Copy a VeloxChem force field's mutable data without its native MPI/output objects."""
    cloned = copy(ff)
    for name in _FORCE_FIELD_DATA:
        if hasattr(ff, name):
            setattr(cloned, name, deepcopy(getattr(ff, name)))
    return cloned


def check_force_field_atoms(ff,
                            *,
                            exempt: int | None = None,
                            label=str) -> None:
    """Raise ValueError unless every atom has finite charge/LJ, nonnegative LJ, and the
    mixing is Lorentz-Berthelot (``comb_rule`` 2).

    ``exempt`` may have negative LJ (a target whose LJ each candidate replaces);
    ``label(index)`` names an atom in the messages.
    """
    for index, atom in ff.atoms.items():
        values = [atom.get(key) for key in ("charge", "sigma", "epsilon")]
        if not all_finite(values):
            raise ValueError(
                f"VeloxChem left atom {label(index)} without finite "
                "charge/LJ parameters")
        if index != exempt and (atom["sigma"] < 0 or atom["epsilon"] < 0):
            raise ValueError(f"negative LJ parameters for {label(index)}")
    if getattr(ff, "comb_rule", 2) != 2:
        raise ValueError("OpenMM NonbondedForce and native MM fitting require "
                         "Lorentz-Berthelot mixing (comb_rule=2)")


def strip_metal_terms(xml_text: str, metals: set[int]):
    """Remove VeloxChem's distance-inferred metal bonded terms from the XML.

    VeloxChem gives every atom its own class (``<index+1>_MOL``) and residue
    atom name, so terms touching a metal are identified exactly. Returns the
    filtered XML, residue atom names, the remaining covalent bonds, and counts.
    """
    root = ET.fromstring(xml_text)
    residue = root.find("Residues").find("Residue")
    names = [atom.get("name") for atom in residue.findall("Atom")]
    type_class = {t.get("name"): t.get("class") for t in root.find("AtomTypes")}
    class_of_index = {
        index: type_class[atom.get("type")]
        for index, atom in enumerate(residue.findall("Atom"))
    }
    metal_names = {names[i] for i in metals}
    metal_classes = {class_of_index[i] for i in metals}
    removed = {"bonds": 0, "angles": 0, "torsions": 0, "impropers": 0}
    bonds = []
    for bond in list(residue.findall("Bond")):
        pair = (bond.get("atomName1"), bond.get("atomName2"))
        if metal_names & set(pair):
            residue.remove(bond)
        else:
            bonds.append(
                tuple(sorted((names.index(pair[0]), names.index(pair[1])))))
    for tag, child, key in (("HarmonicBondForce", "Bond", "bonds"),
                            ("HarmonicAngleForce", "Angle", "angles"),
                            ("PeriodicTorsionForce", "Proper", "torsions"),
                            ("PeriodicTorsionForce", "Improper", "impropers")):
        force = root.find(tag)
        if force is None:
            continue
        for term in list(force.findall(child)):
            classes = {
                value for name, value in term.attrib.items()
                if name.startswith("class")
            }
            if classes & metal_classes:
                force.remove(term)
                removed[key] += 1
    return ET.tostring(root, encoding="unicode"), names, sorted(bonds), removed


def check_no_bond_crosses_cell(molecule, box_nm, bonds, metals):
    """VeloxChem connectivity is not periodic; refuse covalent bonds it would miss.

    Uses the bond rule of ``Molecule.get_connectivity_matrix`` (covalent radii
    times its default ``factor``, or ``H2_factor`` for H-H) on minimum-image
    distances.
    """
    defaults = inspect.signature(Molecule.get_connectivity_matrix).parameters
    factor = defaults["factor"].default
    h2_factor = defaults["H2_factor"].default
    positions = molecule.get_coordinates_in_bohr()
    box_bohr = np.asarray(box_nm) * 10.0 / bohr_in_angstrom()
    radii = molecule.covalent_radii_to_numpy()
    hydrogen = np.array([label == "H" for label in molecule.get_labels()])
    known, metal_set = set(bonds), set(metals)
    for i in range(len(positions) - 1):
        if i in metal_set:
            continue
        j = np.arange(i + 1, len(positions))
        delta = positions[j] - positions[i]
        direct = np.linalg.norm(delta, axis=1)
        image = np.linalg.norm(minimum_image(delta, box_bohr), axis=1)
        limit = (radii[i] + radii[j]) * np.where(hydrogen[i] & hydrogen[j],
                                                 h2_factor, factor)
        for k in np.flatnonzero((image <= limit) & (direct > limit)):
            if int(j[k]) not in metal_set and (i, int(j[k])) not in known:
                raise ValueError(
                    f"atoms {i} and {int(j[k])} are bonded across the periodic "
                    "boundary; VeloxChem connectivity is not periodic, so supply a "
                    "baseline whose covalent units are whole inside the cell")


def nonbonded_force(system):
    """The system's single NonbondedForce."""
    import openmm as mm

    forces = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError("framework systems require exactly one NonbondedForce")
    return forces[0]


def compare_to_baseline(system, baseline: FrameworkBaseline,
                        changed: set[int]) -> list[str]:
    """List every difference outside ``changed`` particles; empty when preserved.

    Checks nonbonded particles, exceptions, and bonded force terms that touch
    only baseline particles, so a candidate cannot alter its surroundings.
    """
    import openmm as mm

    reference = baseline.new_system()
    n_baseline = reference.getNumParticles()
    differences = []
    ref_nb, new_nb = nonbonded_force(reference), nonbonded_force(system)
    for index in range(n_baseline):
        if index in changed:
            continue
        if _values(ref_nb.getParticleParameters(index)) != _values(
                new_nb.getParticleParameters(index)):
            differences.append(f"nonbonded parameters of particle {index}")
        if reference.getParticleMass(index) != system.getParticleMass(index):
            differences.append(f"mass of particle {index}")
    ref_exceptions = {_pair(e): _values(e[2:]) for e in _exceptions(ref_nb)}
    new_exceptions = {
        _pair(e): _values(e[2:])
        for e in _exceptions(new_nb)
        if not set(_pair(e)) & changed and max(_pair(e)) < n_baseline
    }
    if ref_exceptions != new_exceptions:
        differences.append("nonbonded exceptions between baseline particles")
    for kind in (mm.HarmonicBondForce, mm.HarmonicAngleForce,
                 mm.PeriodicTorsionForce):
        ref_terms = _bonded_terms(reference, kind)
        new_terms = [
            term for term in _bonded_terms(system, kind)
            if max(term[0]) < n_baseline
        ]
        if sorted(ref_terms) != sorted(new_terms):
            differences.append(f"{kind.__name__} terms")
    return differences


def _exceptions(force):
    """Return the parameters of every exception of a ``NonbondedForce``."""
    return [
        force.getExceptionParameters(i) for i in range(force.getNumExceptions())
    ]


def _values(quantities) -> tuple[float, ...]:
    """Strip OpenMM units, returning plain floats in OpenMM's internal units."""
    return tuple(
        q._value if hasattr(q, "_value") else float(q) for q in quantities)


def _pair(exception) -> tuple[int, int]:
    """Return the exception's two particle indices in ascending order."""
    return tuple(sorted((exception[0], exception[1])))


def _bonded_terms(system, kind):
    """((atoms...), (parameters...)) for every term of one bonded force class."""
    getters = {
        "HarmonicBondForce": ("getNumBonds", "getBondParameters", 2),
        "HarmonicAngleForce": ("getNumAngles", "getAngleParameters", 3),
        "PeriodicTorsionForce": ("getNumTorsions", "getTorsionParameters", 4)
    }
    count, get, n_atoms = getters[kind.__name__]
    terms = []
    for force in system.getForces():
        if isinstance(force, kind):
            for i in range(getattr(force, count)()):
                parameters = getattr(force, get)(i)
                terms.append((tuple(parameters[:n_atoms]),
                              _values(parameters[n_atoms:])))
    return terms


# ==========================================================================
# framework
# ==========================================================================

FRAMEWORK_PROPERTIES = {
    "distance": "angstrom",
    "angle": "degree",
    "coordination_rmsd": "angstrom",
    "coordination_number": "count",
}
DEFAULT_PROTOCOL = {
    "method": "minimize",  # or "md": average over a short NVT trajectory
    "tolerance_kjmol_nm": 1.0,
    "max_iterations": 5000,
    "temperature_k": 300.0,
    "timestep_fs": 1.0,
    "friction_per_ps": 1.0,
    "equilibration_steps": 2000,
    "production_steps": 10000,
    "sample_interval": 100,
    "seed": 2024,
    "platform": None,  # None picks CPU, then Reference
    "threads": None,
}
DEFAULT_CN_CUTOFF_ANGSTROM = 2.5
DEFAULT_SHELL_SIZE = 6
_ATOMS_NEEDED = {"distance": 2, "angle": 3}


@dataclass
class FrameworkSystem:
    """
    A target metal site inside a cluster or periodic framework.

    ``charges`` selects the baseline charges: ``"resp"`` (VeloxChem RESP of the
    whole structure, runs an SCF), ``"groups"``, an explicit per-atom sequence
    in e, or ``"zero"``. It is required because VeloxChem otherwise assigns
    zero charges silently.

    ``"groups"`` is the scheme for a nonbonded metal model: each metal carries
    its formal charge (``metal_charges``, e.g. ``{"Al": 3}``, by element, atom
    label or index), and the structure is split into the covalent groups left
    when metal-ligand bonds are removed. Each water group (O with two H, no
    charge given) takes the charges of ``water_model`` (``"spc"``, ``"spce"``
    or ``"tip3p"``). Every other group is fitted by RESP as an isolated species
    at its geometry (``MMForceFieldGenerator``, HF/6-31G*), with the integer
    net charge from ``group_charges`` (keyed by any atom label or index in the
    group, e.g. ``{"O1": -1}`` for a hydroxide); unlisted groups are neutral.
    Charge transfer to the metal is then left to the fitted CDAM parameters
    instead of being counted in the charges. The metal and group charges must
    add up to the structure's charge.
    ``baseline_force_field`` accepts an already generated
    ``MMForceFieldGenerator`` for the same structure; it is copied, never changed.
    ``box_vectors`` (3x3, angstrom) makes the system periodic with PME.
    ``representation`` is ``"atom"`` (fit the real metal's LJ) or ``"cdam"``.

    In the CDAM representation the target becomes a central particle plus
    charge-only dummies. ``metal_site_charge`` is the total site charge kept
    while dummy charges change; by default the baseline charge of the target,
    which preserves the system charge. ``metal_dummy_bond`` is ``"constraint"``
    (exact distance) or ``"harmonic"`` (finite flexibility, using the
    candidate's bond force constant). D-M-D harmonic angles hold the geometry,
    as in the aqueous CDAM model. ``c4_partners`` selects the oxygens that feel
    the 12-6-4 C4 term: ``"water"`` (O with exactly two H) or ``"oxygen"``.

    :param structure:
        The structure: a Molecule, an xyz file, or xyz text.
    :param target_metal:
        The target metal (index or atom label).
    :param charges:
        The baseline charges: "resp", "groups", "zero", or a per-atom sequence in e.
    :param charge:
        The molecular charge (required for an xyz structure).
    :param multiplicity:
        The spin multiplicity.
    :param representation:
        The representation of the target, "atom" or "cdam".
    :param baseline_force_field:
        An already generated force field of the same structure.
    :param box_vectors:
        The 3x3 periodic box vectors in angstrom.
    :param protocol:
        The protocol settings (see ``DEFAULT_PROTOCOL``).
    :param metal_site_charge:
        The total charge of the CDAM site.
    :param metal_dummy_bond:
        The metal-dummy bond, "constraint" or "harmonic".
    :param c4_partners:
        The oxygens that feel the C4 term, "water" or "oxygen".
    :param metal_charges:
        The formal metal charges for charges="groups", by element, label or index.
    :param group_charges:
        The net ligand-group charges for charges="groups", by a member atom.
    :param water_model:
        The water model whose charges water groups take for charges="groups".
    :param name:
        The system name.
    :param ostream:
        The output stream.
    """

    structure: Any
    target_metal: int | str
    charges: Any = None
    charge: int | None = None
    multiplicity: int = 1
    representation: str = "atom"
    baseline_force_field: Any = None
    box_vectors: Any = None
    protocol: dict[str, Any] = field(default_factory=dict)
    metal_site_charge: float | None = None
    metal_dummy_bond: str = "constraint"
    c4_partners: str = "water"
    metal_charges: dict[Any, float] | None = None
    group_charges: dict[Any, int] | None = None
    water_model: str = "spce"
    name: str = "framework"
    ostream: Any = None

    def __post_init__(self):
        """Validate options, read the structure and resolve the target metal.

        It is a critical error if the target is not a metal, no baseline
        charge source is chosen, or ``box_vectors`` (angstrom) is not a finite,
        nonsingular 3x3 matrix. The VeloxChem baseline itself is built lazily
        by ``prepare_baseline``.
        """
        self._check_options()
        self.molecule = self._read_structure()
        self.elements = tuple(self.molecule.get_labels())
        self.labels = tuple(atom_labels(self.molecule))
        self.target_index = resolve_atom(self.target_metal, self.labels)
        assert_msg_critical(
            self.target_index in metal_indices(self.molecule),
            f"target {self.target_metal!r} is {self.elements[self.target_index]}, "
            "not a metal")
        assert_msg_critical(
            self.charges is not None or self.baseline_force_field is not None,
            "choose baseline charges: charges='resp' (VeloxChem RESP), "
            "'groups' (formal metal charges with per-group ligand RESP), an "
            "explicit per-atom list, or 'zero'; VeloxChem would otherwise "
            "assign zero charges without RESP")
        if self.box_vectors is not None:
            box = np.asarray(self.box_vectors, dtype=float)
            assert_msg_critical(
                box.shape == (3, 3) and bool(np.isfinite(box).all()) and
                abs(np.linalg.det(box)) >= 1e-6,
                "box_vectors must be a finite, nonsingular 3x3 matrix in angstrom")
            self.box_vectors = box
        self.protocol = dict(self.protocol)
        self.ostream = default_ostream(self.ostream)
        self.baseline_builds = 0
        self._baseline: FrameworkBaseline | None = None

    def _check_options(self):
        """Critical error for any option value outside the supported set."""
        assert_msg_critical(self.representation in ("atom", "cdam"),
                            "representation must be 'atom' or 'cdam'")
        assert_msg_critical(
            self.metal_dummy_bond in ("constraint", "harmonic"),
            "metal_dummy_bond must be 'constraint' or 'harmonic'")
        assert_msg_critical(self.c4_partners in ("water", "oxygen"),
                            "c4_partners must be 'water' or 'oxygen'")
        assert_msg_critical(
            self.water_model in WATER_MODELS_FOR_GROUPS,
            "water_model must be one of " + ", ".join(WATER_MODELS_FOR_GROUPS) +
            " (a 4-site model needs a virtual site the framework cannot add)")
        if isinstance(self.charges, str) and self.charges == "groups":
            assert_msg_critical(
                bool(self.metal_charges),
                "charges='groups' needs metal_charges, e.g. {'Al': 3}")

    def _read_structure(self):
        """A VeloxChem Molecule from the given object, xyz file, or xyz text."""

        molecule = self.structure
        if isinstance(molecule, (str, Path)):
            text = Path(molecule).read_text() if Path(
                str(molecule)).suffix == ".xyz" else str(molecule)
            molecule = Molecule.read_xyz_string(text)
            assert_msg_critical(
                self.charge is not None,
                "give the molecular charge for an xyz structure")
        if self.charge is not None:
            molecule.set_charge(self.charge)
        molecule.set_multiplicity(self.multiplicity)
        return molecule

    # ================================================================ baseline
    def prepare_baseline(self) -> FrameworkBaseline:
        """
        Generates the VeloxChem baseline once and keeps an immutable snapshot.

        :return:
            The baseline snapshot.
        """

        if self._baseline is None:
            self._baseline = self._build_baseline()
            self.baseline_builds += 1
        return self._baseline

    def _build_baseline(self) -> FrameworkBaseline:
        """VeloxChem force field -> stripped XML -> OpenMM system -> snapshot."""
        import openmm as mm
        import openmm.unit as unit

        ff, source = self._generate_force_field()
        # The target may have negative LJ here because each candidate replaces it.
        check_force_field_atoms(ff,
                                exempt=self.target_index,
                                label=self.labels.__getitem__)
        metals = metal_indices(self.molecule)
        xml_text, atom_names, bonds, removed = strip_metal_terms(
            self._openmm_xml(ff), set(metals))
        topology = self._topology(atom_names, bonds)
        positions = self.molecule.get_coordinates_in_angstrom() * ANGSTROM_TO_NM
        box_nm = None
        if self.box_vectors is not None:
            box_nm = self.box_vectors * ANGSTROM_TO_NM
            check_no_bond_crosses_cell(self.molecule, box_nm, bonds, metals)
        system = self._create_system(xml_text, topology, box_nm)
        nonbonded = nonbonded_force(system)
        charges = [
            nonbonded.getParticleParameters(i)[0].value_in_unit(
                unit.elementary_charge) for i in range(len(self.elements))
        ]
        system_xml = mm.XmlSerializer.serialize(system)
        positions = np.array(positions, dtype=float)
        positions.setflags(write=False)
        return FrameworkBaseline(system_xml=system_xml,
                                 positions_nm=positions,
                                 elements=self.elements,
                                 target_index=self.target_index,
                                 metal_indices=metals,
                                 charges=tuple(charges),
                                 bonds=tuple(bonds),
                                 removed_metal_terms=removed,
                                 box_vectors_nm=box_nm,
                                 source=source,
                                 fingerprint=self._fingerprint(
                                     system_xml, positions, box_nm))

    def _generate_force_field(self):
        """(force field, source description) with the chosen baseline charges."""

        n_atoms = len(self.elements)
        if self.baseline_force_field is not None:
            ff = copy_force_field(self.baseline_force_field)
            if len(ff.atoms) != n_atoms:
                raise ValueError(
                    "baseline_force_field does not match the structure's atom count"
                )
            return ff, "supplied VeloxChem force field"
        charges = self.charges
        if isinstance(charges, str) and charges == "groups":
            charges = self._group_scheme_charges()
        ff = MMForceFieldGenerator()
        ff.ostream.mute()
        if isinstance(charges, str):
            if charges not in ("resp", "zero"):
                raise ValueError("charges must be 'resp', 'groups', 'zero', "
                                 "or a per-atom sequence")
            resp = charges == "resp"
            source = (
                "VeloxChem MMForceFieldGenerator, RESP charges" if resp else
                "VeloxChem MMForceFieldGenerator, zero charges (explicitly requested)"
            )
        else:
            values = np.asarray(charges, dtype=float)
            if values.shape != (n_atoms,) or not np.isfinite(values).all():
                raise ValueError(
                    f"charges must be {n_atoms} finite values in e")
            if not np.isclose(
                    values.sum(), self.molecule.get_charge(), atol=1e-4):
                raise ValueError(
                    f"charges sum to {values.sum():.4f} e but the structure's "
                    f"charge is {self.molecule.get_charge()}")
            ff.partial_charges = values.copy()
            resp = True  # preset charges are kept; VeloxChem skips its own RESP
            source = ("VeloxChem MMForceFieldGenerator, formal metal charges "
                      f"with per-group RESP ligands and {self.water_model} water"
                      if isinstance(self.charges, str) else
                      "VeloxChem MMForceFieldGenerator, user-supplied charges")
        try:
            # RESP writes SCF files named after molecule_name; keep them out
            # of the working directory.
            with tempfile.TemporaryDirectory() as directory:
                name = ff.molecule_name
                ff.molecule_name = str(Path(directory) / name)
                try:
                    ff.create_topology(self.molecule, resp=resp)
                finally:
                    ff.molecule_name = name
        except Exception as exc:
            raise ValueError(
                f"VeloxChem could not parameterize the {self.name} structure: "
                f"{type(exc).__name__}: {exc}") from exc
        if not isinstance(charges, str):
            kept = np.array([ff.atoms[i]["charge"] for i in range(n_atoms)])
            if not np.allclose(kept, np.round(charges, 6), atol=1e-5):
                raise ValueError("VeloxChem changed the supplied charges")
        return ff, source

    def _group_scheme_charges(self) -> list[float]:
        """Per-atom charges of the ``"groups"`` scheme (see the class docstring)."""
        from ..waterparameters import get_water_parameters

        n_atoms = len(self.elements)
        metals = metal_indices(self.molecule)
        groups = covalent_groups(self.molecule)
        group_of = {atom: g for g, atoms in enumerate(groups) for atom in atoms}
        charges = np.zeros(n_atoms)

        # formal metal charges: by element first, then by label or index
        metal_set, assigned = set(metals), {}
        for key, value in self.metal_charges.items():
            if isinstance(key, str) and any(
                    self.elements[m].upper() == key.strip().upper() for m in metals):
                targets = [m for m in metals
                           if self.elements[m].upper() == key.strip().upper()]
            else:
                targets = [resolve_atom(key, self.labels)]
                if targets[0] not in metal_set:
                    raise ValueError(f"metal_charges key {key!r} is not a metal")
            for m in targets:
                assigned[m] = float(value)
        missing = [self.labels[m] for m in metals if m not in assigned]
        if missing:
            raise ValueError("metal_charges has no formal charge for "
                             + ", ".join(missing))
        for m, value in assigned.items():
            charges[m] = value

        # integer net charge of each ligand group, named by any member atom
        totals = {}
        for key, value in (self.group_charges or {}).items():
            atom = resolve_atom(key, self.labels)
            if atom in metal_set:
                raise ValueError(f"group_charges key {key!r} is a metal; "
                                 "use metal_charges")
            if float(value) != int(value):
                raise ValueError(f"group charge {value!r} of {key!r} is not an integer")
            group = group_of[atom]
            if group in totals and totals[group] != int(value):
                raise ValueError(f"group_charges gives two charges to the group of {key!r}")
            totals[group] = int(value)

        water = get_water_parameters()[self.water_model]
        coordinates = self.molecule.get_coordinates_in_angstrom()
        for g, atoms in enumerate(groups):
            elements = [self.elements[i] for i in atoms]
            total = totals.get(g, 0)
            names = " ".join(self.labels[i] for i in atoms)
            if sorted(elements) == ["H", "H", "O"] and g not in totals:
                for i in atoms:
                    charges[i] = water["ow" if self.elements[i] == "O" else "hw"]["charge"]
                how = f"{self.water_model} water"
            elif len(atoms) == 1:
                charges[atoms[0]] = float(total)
                how = "single atom"
            else:
                charges[atoms] = self._isolated_group_resp(elements,
                                                           coordinates[atoms],
                                                           total, names)
                how = "RESP of the isolated group"
            self.ostream.print_info(f"Charge group {names}: net {total:+d}, {how}")

        target = self.molecule.get_charge()
        if not np.isclose(charges.sum(), target, atol=1e-4):
            raise ValueError(
                f"metal and group charges sum to {charges.sum():.4f} e but the "
                f"structure's charge is {target}; set group_charges for the "
                "charged ligands")
        self.ostream.flush()
        return charges.tolist()

    @staticmethod
    def _isolated_group_resp(elements, coordinates, total, names) -> np.ndarray:
        """RESP charges of one ligand group as an isolated closed-shell species."""
        group = Molecule(list(elements), np.asarray(coordinates), "angstrom")
        group.set_charge(total)
        electrons = sum(group.get_element_ids()) - total
        if electrons % 2:
            raise ValueError(
                f"group {names} with net charge {total:+d} has an odd number of "
                "electrons; give its closed-shell charge in group_charges")
        group.set_multiplicity(1)
        ff = MMForceFieldGenerator()
        ff.ostream.mute()
        with tempfile.TemporaryDirectory() as directory:
            ff.molecule_name = str(Path(directory) / "group")
            try:
                ff.create_topology(group, resp=True)
            except Exception as exc:
                raise ValueError(
                    f"RESP failed for the group {names}: "
                    f"{type(exc).__name__}: {exc}") from exc
        return np.array([ff.atoms[i]["charge"] for i in range(len(elements))])

    @staticmethod
    def _openmm_xml(ff) -> str:
        """Return the force field's OpenMM XML text, written via a temporary directory."""
        with tempfile.TemporaryDirectory() as directory:
            ff.write_openmm_files(str(Path(directory) / "baseline"))
            return (Path(directory) / "baseline.xml").read_text()

    def _topology(self, atom_names, bonds):
        """Build a one-residue (``MOL``) OpenMM topology with the given atom names and bonds."""
        import openmm.app as app

        topology = app.Topology()
        residue = topology.addResidue("MOL", topology.addChain("A"))
        atoms = [
            topology.addAtom(atom_names[i],
                             app.element.get_by_symbol(self.elements[i]),
                             residue) for i in range(len(self.elements))
        ]
        for left, right in bonds:
            topology.addBond(atoms[left], atoms[right])
        return topology

    def _create_system(self, xml_text, topology, box_nm):
        """OpenMM system with no constraints; PME in a cell, no cutoff otherwise."""
        import openmm.app as app
        import openmm.unit as unit

        options = {
            "constraints": None,
            "rigidWater": False,
            "removeCMMotion": False
        }
        if box_nm is not None:
            topology.setPeriodicBoxVectors(box_nm * unit.nanometer)
            cutoff = min(1.0, 0.5 * box_widths(box_nm).min() - 1e-3)
            options.update(nonbondedMethod=app.PME,
                           nonbondedCutoff=cutoff * unit.nanometer)
        else:
            options.update(nonbondedMethod=app.NoCutoff)
        system = app.ForceField(io.StringIO(xml_text)).createSystem(
            topology, **options)
        nonbonded = nonbonded_force(system)
        for index in range(nonbonded.getNumExceptions()):
            left, right, *_ = nonbonded.getExceptionParameters(index)
            if self.target_index in (left, right):
                raise ValueError(
                    "baseline still couples the target metal through bonded "
                    "exceptions; metal-ligand terms were not removed")
        return system

    def _fingerprint(self, system_xml, positions, box_nm) -> str:
        """Return a short SHA-256 digest of system, positions, target and box."""
        digest = hashlib.sha256()
        for part in (system_xml, positions.tobytes(), str(self.target_index),
                     "" if box_nm is None else np.asarray(box_nm).tobytes()):
            digest.update(part if isinstance(part, bytes) else part.encode())
        return digest.hexdigest()[:16]

    # ============================================================== candidates
    def build(self, candidate):
        """
        Builds a fresh system for one CDAM candidate.

        ``candidate`` is a ``CDAMCandidate`` point. In the atom representation
        only its metal sigma/epsilon (and C4) are used; the metal charge stays
        at its baseline value. In the CDAM representation dummies are added
        around the target. Every other particle, exception, and bonded term is
        the baseline's; a fresh System per candidate lets constraints change.

        :param candidate:
            The CDAM candidate.

        :return:
            The tuple (system, positions_nm, site).
        """

        baseline = self.prepare_baseline()
        system = baseline.new_system()
        positions = np.array(baseline.positions_nm, dtype=float)
        target = baseline.target_index
        nonbonded = nonbonded_force(system)
        openmm_site = candidate.openmm_parameters()
        sigma, epsilon = openmm_site["sigma"], openmm_site["epsilon"]
        if candidate.metal != self.elements[target]:
            raise ValueError(
                f"candidate metal {candidate.metal} does not match target "
                f"{self.elements[target]}")
        if self.representation == "atom":
            nonbonded.setParticleParameters(target, baseline.charges[target],
                                            sigma, epsilon)
            site = {
                "metal": target,
                "dummies": [],
                "metal_site_charge": baseline.charges[target],
                "central_charge": baseline.charges[target]
            }
        else:
            site = self._add_cdam_site(system, nonbonded, candidate, sigma,
                                       epsilon)
            positions = np.vstack([positions, site.pop("dummy_positions_nm")])
        if candidate.potential_form == "lj12_6_4":
            self._add_c4(system, nonbonded, candidate, site)
        return system, positions, site

    def _metal_site_charge(self) -> float:
        """Total CDAM site charge in e: ``metal_site_charge`` or the target's baseline charge."""
        baseline = self.prepare_baseline()
        if self.metal_site_charge is None:
            return baseline.charges[baseline.target_index]
        return float(self.metal_site_charge)

    def first_shell(self, size: int) -> list[int]:
        """
        Finds the nearest non-hydrogen, non-metal atoms to the target.

        :param size:
            The number of atoms.

        :return:
            The atom indices.
        """

        baseline = self.prepare_baseline()
        return nearest_heavy_atoms(baseline.positions_nm, baseline.target_index,
                                   self.elements, baseline.metal_indices,
                                   baseline.box_vectors_nm, size)

    def dummy_directions(self, candidate) -> np.ndarray:
        """
        Rotates the template dummy directions onto the first coordination shell.

        The rotation minimizing the RMSD between template directions and the
        observed metal-ligand unit vectors (with the best atom assignment) is
        used, so dummies point at the coordinating atoms of the input structure.

        :param candidate:
            The CDAM candidate (its geometry gives the template).

        :return:
            The rotated unit directions.
        """

        template = np.asarray(candidate.geometries[candidate.geometry],
                              dtype=float)
        baseline = self.prepare_baseline()
        target = baseline.target_index
        shell = self.first_shell(len(template))
        observed = minimum_image(
            baseline.positions_nm[shell] - baseline.positions_nm[target],
            baseline.box_vectors_nm)
        observed /= np.linalg.norm(observed, axis=1)[:, None]
        return align_template(template, observed)

    def _add_cdam_site(self, system, nonbonded, candidate, sigma, epsilon):
        """Replace the target by a central particle plus charge-only dummies."""
        import openmm as mm

        baseline = self.prepare_baseline()
        target = baseline.target_index
        directions = self.dummy_directions(candidate)
        n_dummies = len(directions)
        site_charge = self._metal_site_charge()
        dummy_charge = float(candidate.dummy_charge)
        central_charge = site_charge - n_dummies * dummy_charge
        metal_mass = system.getParticleMass(target)._value
        if metal_mass <= n_dummies * candidate.dummy_mass:
            raise ValueError("metal mass must exceed the total dummy mass")
        system.setParticleMass(target,
                               metal_mass - n_dummies * candidate.dummy_mass)
        nonbonded.setParticleParameters(target, central_charge, sigma, epsilon)
        dummies = []
        for _ in range(n_dummies):
            dummies.append(system.addParticle(candidate.dummy_mass))
            nonbonded.addParticle(dummy_charge, 0.0, 0.0)
        if any(
                isinstance(other, mm.CustomNonbondedForce)
                for other in system.getForces()):
            raise ValueError(
                "baseline contains a CustomNonbondedForce of unknown meaning")
        site = [target, *dummies]
        for i, left in enumerate(site):  # internal Coulomb/LJ fully excluded
            for right in site[i + 1:]:
                nonbonded.addException(left, right, 0.0, 1.0, 0.0)
        distance = candidate.openmm_parameters()["distance"]
        self._add_site_restraints(system, candidate, target, dummies,
                                  directions, distance)
        dummy_positions = (baseline.positions_nm[target] +
                           distance * directions)
        return {
            "metal": target,
            "dummies": dummies,
            "metal_site_charge": site_charge,
            "central_charge": central_charge,
            "dummy_charge": dummy_charge,
            "metal_dummy_distance_nm": distance,
            "dummy_positions_nm": dummy_positions
        }

    def _add_site_restraints(self, system, candidate, target, dummies,
                             directions, distance):
        """Metal-dummy constraints (or harmonic bonds) and D-M-D harmonic angles."""
        import openmm as mm

        if self.metal_dummy_bond == "constraint":
            for dummy in dummies:
                system.addConstraint(target, dummy, distance)
        else:
            bonds = mm.HarmonicBondForce()
            bonds.setName("CDAM_metal_dummy")
            k = candidate.openmm_parameters()["bond_k"]
            for dummy in dummies:
                bonds.addBond(target, dummy, distance, k)
            system.addForce(bonds)
        angles = mm.HarmonicAngleForce()
        angles.setName("CDAM_geometry")
        k_angle = candidate.openmm_parameters()["angle_k"]
        for i in range(len(dummies)):
            for j in range(i + 1, len(dummies)):
                theta = np.arccos(
                    np.clip(directions[i] @ directions[j], -1.0, 1.0))
                angles.addAngle(dummies[i], target, dummies[j], float(theta),
                                k_angle)
        system.addForce(angles)

    def _add_c4(self, system, nonbonded, candidate, site):
        """Add -C4/r^4 between the metal and selected oxygens, honoring exclusions."""
        partners = self.c4_oxygens()
        c4 = candidate.openmm_parameters()["c4"]
        force = new_c4_force(system, nonbonded, system.getNumParticles())
        set_c4_pairs(force, system.getNumParticles(), site["metal"], partners,
                     c4)

    def c4_oxygens(self) -> set[int]:
        """
        Finds the oxygens that feel the C4 term: water O (two H neighbours) or
        every O.

        :return:
            The oxygen indices.
        """

        baseline = self.prepare_baseline()
        neighbours = {i: [] for i in range(len(self.elements))}
        for left, right in baseline.bonds:
            neighbours[left].append(right)
            neighbours[right].append(left)
        oxygens = {
            i for i, element in enumerate(self.elements) if element == "O"
        }
        if self.c4_partners == "water":
            oxygens = {
                i for i in oxygens if len(neighbours[i]) == 2 and all(
                    self.elements[j] == "H" for j in neighbours[i])
            }
        return oxygens

    def ineffective_parameters(self, specs) -> list[str]:
        """
        Finds adjustable parameters that cannot change this system's interactions.

        :param specs:
            The parameter specifications.

        :return:
            The names of the ineffective optimizable parameters.
        """

        unused = set()
        if self.representation == "atom":
            unused |= {
                "dummy_charge", "metal_dummy_distance", "bond_force_constant",
                "angle_force_constant"
            }
        elif self.metal_dummy_bond == "constraint":
            unused.add("bond_force_constant")
        return [
            spec.name
            for spec in specs
            if spec.optimizable and spec.name in unused
        ]

    # ============================================================== properties
    def check_targets(self, targets) -> None:
        """
        Validates property names and atom selections before any calculation.

        Raises ValueError for an unsupported property, a wrong number of atoms,
        or an unknown method.

        :param targets:
            The targets of this system.
        """

        for target in targets:
            prop = target.property_name
            if prop not in FRAMEWORK_PROPERTIES:
                raise ValueError(
                    f"target {target.name!r}: {self.name} supports "
                    f"{sorted(FRAMEWORK_PROPERTIES)}, not {prop!r}")
            atoms = target.metadata.get("atoms")
            needed = _ATOMS_NEEDED.get(prop)
            if needed is not None:
                if atoms is None or len(atoms) != needed:
                    raise ValueError(
                        f"target {target.name!r}: {prop} needs {needed} atoms")
                [resolve_atom(atom, self.labels) for atom in atoms]
            method = target.metadata.get(
                "method", self.protocol.get("method", "minimize"))
            if method not in ("minimize", "md"):
                raise ValueError(
                    f"target {target.name!r}: method must be 'minimize' or 'md'"
                )

    def calculate(self, candidate, targets, protocol=None) -> CalculationResult:
        """
        Builds the candidate system, relaxes or samples it, and measures targets.

        Targets may ask for ``method`` "minimize" or "md"; each requested
        method runs once. ``protocol`` temporarily overrides protocol settings
        (used for cheap screening). Failures return a failed result.

        :param candidate:
            The CDAM candidate.
        :param targets:
            The targets of this system.
        :param protocol:
            The temporary protocol overrides.

        :return:
            The calculation result.
        """

        targets = list(targets)
        started = time.monotonic()
        try:
            self.check_targets(targets)
            baseline = self.prepare_baseline()
            system, positions, site = self.build(candidate)
            preserved = compare_to_baseline(system,
                                            baseline,
                                            changed={baseline.target_index})
            if preserved:
                raise RuntimeError(
                    f"non-target parameters changed: {preserved}")
            protocol = {**DEFAULT_PROTOCOL, **self.protocol, **(protocol or {})}
            methods = {
                target.metadata.get("method", protocol["method"])
                for target in targets
            }
            frames = {
                method: self._sample(system, positions, protocol, method)
                for method in sorted(methods)
            }
        except Exception as exc:
            return CalculationResult(
                failed=True,
                failure_reason=f"{self.name}: {type(exc).__name__}: {exc}",
                metadata={"system": self.name})
        values, metadata = {}, {}
        for target in targets:
            method = target.metadata.get("method", protocol["method"])
            samples = [
                self.measure(target, xyz, site) for xyz in frames[method]
            ]
            values[target.name] = float(np.mean(samples))
            metadata[target.name] = self._target_metadata(
                target, method, samples, baseline)
        raw_data = {
            "system": self.name,
            "site": site,
            "runtime_seconds": time.monotonic() - started
        }
        return CalculationResult(values=values,
                                 metadata=metadata,
                                 raw_data=raw_data)

    def _target_metadata(self, target, method, samples,
                         baseline) -> dict[str, Any]:
        """Return the per-target metadata dict (unit, method, frame count, spread, fingerprint)."""
        minimized = method == "minimize"
        spread = float(np.std(samples)) if len(samples) > 1 else 0.0
        return {
            "system": self.name,
            "property": target.property_name,
            "unit": FRAMEWORK_PROPERTIES[target.property_name],
            "method": "OpenMM minimization" if minimized else "OpenMM NVT MD",
            "geometry": "minimized" if minimized else "md_average",
            "n_frames": len(samples),
            "standard_deviation": spread,
            "baseline_fingerprint": baseline.fingerprint,
            "representation": self.representation,
            "non_target_preserved": True,
            "quality_status": "usable",
        }

    def _sample(self, system, positions, protocol, method) -> list[np.ndarray]:
        """Minimize, then optionally run NVT MD; returns frames in nm."""
        import openmm as mm
        import openmm.unit as unit

        integrator = mm.LangevinMiddleIntegrator(
            protocol["temperature_k"] * unit.kelvin,
            protocol["friction_per_ps"] / unit.picosecond,
            protocol["timestep_fs"] * unit.femtosecond)
        integrator.setRandomNumberSeed(int(protocol["seed"]))
        platform, properties = _platform(protocol)
        context = mm.Context(system, integrator, platform, properties)
        baseline = self.prepare_baseline()
        if baseline.box_vectors_nm is not None:
            context.setPeriodicBoxVectors(*(baseline.box_vectors_nm *
                                            unit.nanometer))
        context.setPositions(positions * unit.nanometer)
        mm.LocalEnergyMinimizer.minimize(context,
                                         protocol["tolerance_kjmol_nm"],
                                         int(protocol["max_iterations"]))
        state = context.getState(getPositions=True, getEnergy=True)
        if not np.isfinite(state.getPotentialEnergy()._value):
            raise RuntimeError("minimized energy is not finite")
        if method == "minimize":
            return [state.getPositions(asNumpy=True)._value]
        context.setVelocitiesToTemperature(
            protocol["temperature_k"] * unit.kelvin, int(protocol["seed"]))
        integrator.step(int(protocol["equilibration_steps"]))
        frames = []
        interval = int(protocol["sample_interval"])
        for _ in range(max(1, int(protocol["production_steps"]) // interval)):
            integrator.step(interval)
            state = context.getState(getPositions=True, getEnergy=True)
            if not np.isfinite(state.getPotentialEnergy()._value):
                raise RuntimeError("MD energy became nonfinite")
            frames.append(state.getPositions(asNumpy=True)._value)
        return frames

    def measure(self, target, xyz_nm: np.ndarray, site) -> float:
        """
        Measures one property value for one geometry.

        :param target:
            The target (its property and atom selection).
        :param xyz_nm:
            The coordinates in nm.
        :param site:
            The site returned by ``build``.

        :return:
            The value in angstrom, degree, or count.
        """

        baseline = self.prepare_baseline()
        box = baseline.box_vectors_nm
        prop = target.property_name
        atoms = [
            resolve_atom(atom, self.labels)
            for atom in target.metadata.get("atoms", [])
        ]
        if prop == "distance":
            return distance_angstrom(xyz_nm, atoms, box)
        if prop == "angle":
            return angle_degrees(xyz_nm, atoms, box)
        if prop == "coordination_number":
            metal = baseline.target_index
            cutoff = float(
                target.metadata.get("cutoff_angstrom",
                                    DEFAULT_CN_CUTOFF_ANGSTROM))
            heavy = [
                i for i, element in enumerate(self.elements)
                if element != "H" and i != metal and
                i not in baseline.metal_indices
            ]
            return count_within(xyz_nm, metal, heavy, cutoff * ANGSTROM_TO_NM,
                                box)
        return self._coordination_rmsd(target, xyz_nm)

    def _coordination_rmsd(self, target, xyz_nm) -> float:
        """Kabsch RMSD of the metal and its shell against a reference (angstrom).

        ``selection`` is "first_shell" (the ``shell_size`` nearest heavy atoms
        of the input structure; default 6) or explicit atom selectors. The
        reference defaults to the input structure itself (identical atoms);
        ``reference_structure`` (``[[element, x, y, z], ...]`` in angstrom,
        metal first) is matched by element with optimal assignment, so
        symmetry-equivalent ligands in any order are handled.
        """

        baseline = self.prepare_baseline()
        metal = baseline.target_index
        selection = target.metadata.get("selection", "first_shell")
        if selection == "first_shell":
            shell = self.first_shell(
                int(target.metadata.get("shell_size", DEFAULT_SHELL_SIZE)))
        else:
            shell = [resolve_atom(atom, self.labels) for atom in selection]
        atoms = [metal, *shell]
        box = baseline.box_vectors_nm
        local = minimum_image(xyz_nm[atoms] - xyz_nm[metal],
                              box) / ANGSTROM_TO_NM
        reference = target.metadata.get("reference_structure")
        if reference is None:
            ref = minimum_image(
                baseline.positions_nm[atoms] - baseline.positions_nm[metal],
                box) / ANGSTROM_TO_NM
            return kabsch_rmsd(local, ref)
        ref_elements = [str(row[0]).capitalize() for row in reference]
        ref = np.array([[float(v) for v in row[1:4]] for row in reference])
        elements = [self.elements[i] for i in atoms]
        if sorted(ref_elements) != sorted(elements):
            raise ValueError(
                f"reference_structure elements {ref_elements} do not match the "
                f"selected atoms {elements}")
        return element_matched_rmsd(local, elements, ref - ref[0], ref_elements)


def _platform(protocol):
    """The protocol's OpenMM platform, else CPU, else Reference."""
    import openmm as mm

    if protocol.get("platform"):
        names = [protocol["platform"]]
    else:
        names = ["CPU", "Reference"]
    for name in names:
        try:
            platform = mm.Platform.getPlatformByName(name)
        except Exception:
            continue
        properties = {}
        if name == "CPU" and protocol.get("threads"):
            properties["Threads"] = str(protocol["threads"])
        return platform, properties
    raise ValueError(f"no OpenMM platform among {names}")

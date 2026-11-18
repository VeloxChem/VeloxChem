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
"""A metal node (a capped MOF node or a metal cluster) as a fit property.

``MetalNode`` reads the node geometry and classifies it the way a chemist
would: the metal atoms, and the residues left once the metal-ligand bonds are
removed (VeloxChem's connectivity). Each residue is typed automatically --
oxo, hydroxide, water, carboxylate (formate and acetate caps included),
halide -- which gives its formal charge, its donor atoms and whether it
bridges several metals. Anything the rules do not recognize is set with
``set_residue_type``.

It is not part of the AutoCDAM driver. Given as the ``system`` of one
property, it hands the driver two ordinary inputs (see
``multisystem.expand_system_objects``): ``system()``, the framework system
specification, and ``targets()``, the grouped geometry targets of the fitted
metal -- the metal-donor distance per residue type and the RMSD of the
donor-metal-donor angles (cis and trans) against the input geometry.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import json

import numpy as np

from ..molecule import Molecule
from ..errorhandler import assert_msg_critical
from .framework import WATER_MODELS_FOR_GROUPS, covalent_groups, metal_indices
from .geometry import atom_labels, resolve_atom

# formal charge of each automatically recognized residue type; a carboxylate
# residue carries -1 per COO group
RESIDUE_CHARGES = {
    "oxo": -2,
    "hydroxide": -1,
    "water": 0,
    "carboxylate": -1,
    "fluoride": -1,
    "chloride": -1,
    "bromide": -1,
    "iodide": -1,
}
_HALIDES = {"F": "fluoride", "Cl": "chloride", "Br": "bromide", "I": "iodide"}
_CAPS = {"CHO2": "formate", "C2H3O2": "acetate"}
TRANS_ANGLE_DEGREES = 135.0
DEFAULT_TOLERANCES = {"distance": 0.05, "angle": 10.0}
_TARGET_KEYS = {"system", "reference", "tolerance", "weight", "role", "method"}


class Residue:
    """
    One residue of a metal node: a covalent group once metal bonds are removed.

    :param number:
        The residue number (1-based, in order of the first atom).
    :param atoms:
        The atom indices.
    :param formula:
        The Hill formula.
    :param kind:
        The residue type, or None when it is not recognized.
    :param charge:
        The formal charge, or None when it is not known.
    :param donors:
        The atoms bonded to a metal.
    :param metals:
        The metal atoms the residue binds.
    """

    def __init__(self, number, atoms, formula, kind, charge, donors, metals):
        self.number = number
        self.atoms = atoms
        self.formula = formula
        self.kind = kind
        self.charge = charge
        self.donors = donors
        self.metals = metals

    @property
    def bridging(self) -> int:
        """The number of metals the residue binds when more than one, else 0."""
        return len(self.metals) if len(self.metals) > 1 else 0

    @property
    def label(self) -> str:
        """The type with its bridging prefix, e.g. 'mu2-hydroxide'."""
        if self.kind is None:
            return "unknown"
        return f"mu{self.bridging}-{self.kind}" if self.bridging else self.kind


class MetalNode:
    """
    A metal node, typed into residues, to fit a metal inside a framework.

    Instance variables
        - molecule: The VeloxChem Molecule of the node.
        - metal: The index of the fitted metal.
        - residues: The residues (see ``Residue``).

    :param molecule:
        The node: a Molecule (its charge is the node charge).
    :param metal:
        The fitted metal (element, atom label or index); default the only
        metal atom.
    :param charge:
        The node charge; default the molecule's charge.
    :param name:
        The system name.
    :param charges:
        The framework charges: "groups" (default; formal metal charges and
        per-residue RESP with the residue charges), "resp", a per-atom list in
        e (used as given), or a charge file (json with "charges_e",
        "metal_charges", "group_charges" and "water_model") that is checked
        against the group scheme.
    :param metal_charges:
        The formal metal charges by element, label or index; default from the
        node charge minus the residue charges.
    """

    def __init__(self,
                 molecule,
                 metal=None,
                 charge=None,
                 name="node",
                 charges="groups",
                 metal_charges=None):
        assert_msg_critical(
            isinstance(molecule, Molecule),
            "MetalNode: give a VeloxChem Molecule (or use read_xyz_file)")
        if charge is not None:
            molecule.set_charge(charge)
        self.molecule = molecule
        self.name = str(name)
        self.charges = charges
        self.metal_charges = dict(metal_charges) if metal_charges else None
        self.elements = tuple(molecule.get_labels())
        self.labels = tuple(atom_labels(molecule))
        self.metals = metal_indices(molecule)
        assert_msg_critical(bool(self.metals), "MetalNode: no metal atom found")
        self.metal = self._resolve_metal(metal)
        self.residues = self._classify()

    @classmethod
    def read_xyz_file(cls, path, metal=None, charge=0, **kwargs):
        """
        Reads a node from an xyz file.

        :param path:
            The xyz file.
        :param metal:
            The fitted metal (element, atom label or index).
        :param charge:
            The node charge.

        :return:
            The metal node.
        """

        molecule = Molecule.read_xyz_file(str(path))
        return cls(molecule, metal=metal, charge=charge, **kwargs)

    # ------------------------------------------------------------ classifying

    def _resolve_metal(self, metal) -> int:
        """The index of the fitted metal."""
        if metal is None:
            assert_msg_critical(
                len(self.metals) == 1,
                "MetalNode: the node has several metals "
                f"({', '.join(self.labels[m] for m in self.metals)}); choose "
                "one with metal=")
            return self.metals[0]
        if isinstance(metal, str):
            same = [m for m in self.metals
                    if self.elements[m].upper() == metal.strip().upper()]
            if len(same) == 1:
                return same[0]
        index = resolve_atom(metal, self.labels)
        assert_msg_critical(index in self.metals,
                            f"MetalNode: {metal!r} is not a metal atom")
        return index

    def _classify(self) -> list[Residue]:
        """The typed residues of the node."""
        connectivity = np.asarray(self.molecule.get_connectivity_matrix(),
                                  dtype=bool)
        self._connectivity = connectivity
        metals = set(self.metals)
        residues = []
        for number, atoms in enumerate(covalent_groups(self.molecule), 1):
            donors = [a for a in atoms
                      if any(connectivity[a, m] for m in metals)]
            bound = sorted({m for a in donors for m in metals
                            if connectivity[a, m]})
            kind, charge = self._residue_type(atoms, connectivity)
            residues.append(
                Residue(number, atoms, self._formula(atoms), kind, charge,
                        donors, bound))
        return residues

    def _formula(self, atoms) -> str:
        """The Hill formula of some atoms."""
        counts = Counter(self.elements[a] for a in atoms)
        order = (["C", "H"] if "C" in counts else []) + sorted(
            e for e in counts if "C" not in counts or e not in ("C", "H"))
        return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}"
                       for e in order)

    def _residue_type(self, atoms, connectivity):
        """(type, formal charge) of a residue, or (None, None)."""
        elements = sorted(self.elements[a] for a in atoms)
        if len(atoms) == 1:
            if elements == ["O"]:
                return "oxo", RESIDUE_CHARGES["oxo"]
            if elements[0] in _HALIDES:
                kind = _HALIDES[elements[0]]
                return kind, RESIDUE_CHARGES[kind]
            return None, None
        if elements == ["H", "O"]:
            return "hydroxide", RESIDUE_CHARGES["hydroxide"]
        if elements == ["H", "H", "O"]:
            return "water", RESIDUE_CHARGES["water"]
        inside = set(atoms)

        def neighbours(atom):
            return [int(b) for b in np.flatnonzero(connectivity[atom])
                    if int(b) in inside]

        # carboxylate: a C with exactly two terminal O (no H on either O)
        groups = 0
        for atom in atoms:
            if self.elements[atom] != "C":
                continue
            oxygens = [b for b in neighbours(atom) if self.elements[b] == "O"]
            if len(oxygens) == 2 and all(neighbours(o) == [atom]
                                         for o in oxygens):
                groups += 1
        if groups:
            return "carboxylate", groups * RESIDUE_CHARGES["carboxylate"]
        return None, None

    def set_residue_type(self, residue, kind, charge=None) -> None:
        """
        Sets the type (and formal charge) of a residue the rules missed.

        :param residue:
            The residue number (see ``print_residues``) or one of its atoms
            (label or index).
        :param kind:
            The residue type, e.g. "phosphonate".
        :param charge:
            The formal charge; default the charge of a known type.
        """

        target = self._residue(residue)
        if charge is None:
            assert_msg_critical(
                kind in RESIDUE_CHARGES,
                f"MetalNode: give the formal charge of a {kind!r} residue")
            charge = RESIDUE_CHARGES[kind]
        assert_msg_critical(
            float(charge) == int(charge),
            f"MetalNode: the residue charge {charge!r} is not an integer")
        target.kind, target.charge = str(kind), int(charge)

    def _residue(self, residue) -> Residue:
        """A residue by number or by one of its atoms."""
        if isinstance(residue, Residue):
            return residue
        if isinstance(residue, (int, np.integer)) and not isinstance(
                residue, bool) and 1 <= residue <= len(self.residues):
            return self.residues[int(residue) - 1]
        atom = resolve_atom(residue, self.labels)
        return next(r for r in self.residues if atom in r.atoms)

    # ------------------------------------------------------------- reporting

    def print_residues(self, ostream=None) -> None:
        """Prints the metals and the typed residues."""
        lines = [
            f"Metal node {self.name!r}: {self._formula(range(len(self.elements)))}, "
            f"charge {self.molecule.get_charge():+.0f}, fitted metal "
            f"{self.labels[self.metal]}",
            f"{'Residue':>7}  {'Type':<23}{'Formula':<10}{'Charge':>6}  "
            f"{'Donors':<14}{'Binds':<12}{'M-donor (A)'}",
        ]
        xyz = self.molecule.get_coordinates_in_angstrom()
        for r in self.residues:
            kind = r.label + (f" ({_CAPS[r.formula]})"
                              if r.formula in _CAPS else "")
            distances = " ".join(
                f"{np.linalg.norm(xyz[d] - xyz[m]):.3f}" for d in r.donors
                for m in r.metals if self._bonded(d, m))
            charge = "?" if r.charge is None else f"{r.charge:+d}"
            lines.append(
                f"{r.number:>7}  {kind:<23}{r.formula:<10}{charge:>6}  "
                f"{' '.join(self.labels[d] for d in r.donors) or '-':<14}"
                f"{' '.join(self.labels[m] for m in r.metals) or '-':<12}"
                f"{distances}")
        unknown = [r for r in self.residues if r.kind is None]
        if unknown:
            lines.append("Unrecognized residues: set them with "
                         "set_residue_type(number, type, charge)")
        if ostream is None:
            print("\n".join(lines))
            return
        for line in lines:
            ostream.print_info(line)
        ostream.print_blank()
        ostream.flush()

    def _bonded(self, a, b) -> bool:
        return bool(self._connectivity[a, b])

    # ------------------------------------------------- inputs for the driver

    def check(self) -> None:
        """Critical error for unrecognized residues or charges that do not add up."""
        unknown = [r for r in self.residues if r.kind is None]
        assert_msg_critical(
            not unknown,
            f"MetalNode {self.name!r}: residues not recognized: " + "; ".join(
                f"{r.number} {r.formula} ({' '.join(self.labels[a] for a in r.atoms)})"
                for r in unknown) +
            ". Set each with node.set_residue_type(number, type, charge)")
        metal_charges = self._metal_charges()
        total = sum(r.charge for r in self.residues) + sum(
            metal_charges.values())
        assert_msg_critical(
            np.isclose(total, self.molecule.get_charge()),
            f"MetalNode {self.name!r}: metal and residue charges add up to "
            f"{total:+g}, not the node charge "
            f"{self.molecule.get_charge():+g}; check the residue types and "
            "metal_charges")

    def _metal_charges(self) -> dict[int, float]:
        """Formal charge of each metal atom."""
        if self.metal_charges:
            charges = {}
            for key, value in self.metal_charges.items():
                if isinstance(key, str) and not key[-1].isdigit():
                    for m in self.metals:
                        if self.elements[m].upper() == key.strip().upper():
                            charges[m] = float(value)
                else:
                    charges[resolve_atom(key, self.labels)] = float(value)
            missing = [self.labels[m] for m in self.metals if m not in charges]
            assert_msg_critical(
                not missing, f"MetalNode {self.name!r}: metal_charges has no "
                f"formal charge for {', '.join(missing)}")
            return charges
        rest = self.molecule.get_charge() - sum(
            r.charge for r in self.residues if r.charge is not None)
        share = rest / len(self.metals)
        assert_msg_critical(
            len({self.elements[m] for m in self.metals}) == 1 and
            float(share).is_integer(),
            f"MetalNode {self.name!r}: cannot split the metal charge {rest:+g} "
            "over the metals; give metal_charges")
        return {m: share for m in self.metals}

    def system(self, water_model=None, ostream=None) -> dict:
        """
        The framework system specification of the node (CDAM representation).

        :param water_model:
            The water model of the study; water residues take its charges.
        :param ostream:
            The output stream for the residue table.

        :return:
            The keyword settings of a framework system.
        """

        self.check()
        if ostream is not None:
            self.print_residues(ostream)
        water_model = water_model.lower() if water_model else None
        waters = any(r.kind == "water" for r in self.residues)
        if water_model is None or water_model not in WATER_MODELS_FOR_GROUPS:
            assert_msg_critical(
                not waters or water_model is None,
                f"MetalNode {self.name!r}: the node has water residues, which "
                f"take {water_model!r} charges; the framework supports "
                + ", ".join(WATER_MODELS_FOR_GROUPS))
            water_model = WATER_MODELS_FOR_GROUPS[1]
        spec = {
            "type": "framework",
            "structure": self.molecule,
            "target_metal": self.labels[self.metal],
            "representation": "cdam",
            "charge": int(round(self.molecule.get_charge())),
            "water_model": water_model,
        }
        metal_charges = {self.labels[m]: q
                         for m, q in self._metal_charges().items()}
        # water residues take the water model's charges; other residues
        # are listed by their first donor (or first atom)
        group_charges = {
            self.labels[(r.donors or r.atoms)[0]]: r.charge
            for r in self.residues if r.kind != "water" and r.charge != 0
        }
        charges = self.charges
        if isinstance(charges, (str, Path)) and str(charges).endswith(".json"):
            # a saved group-scheme result: used as is and checked against
            # the scheme with the water model it was made for
            data = json.loads(Path(charges).read_text())
            spec.update(charges=data["charges_e"],
                        metal_charges=data.get("metal_charges", metal_charges),
                        group_charges=data.get("group_charges", group_charges),
                        water_model=data.get("water_model", water_model))
        elif isinstance(charges, str):
            # "groups": formal metal charges and RESP of each residue
            spec.update(charges=charges,
                        metal_charges=metal_charges,
                        group_charges=group_charges)
        else:
            spec.update(charges=list(charges))  # a per-atom list, as given
        return spec

    def targets(self, name, values) -> dict[str, dict]:
        """
        The grouped geometry targets of the fitted metal.

        One distance target per residue type that binds the metal (the mean
        metal-donor distance) and the RMSD of the cis and trans
        donor-metal-donor angles against the input geometry. ``reference``
        maps a residue type (``"hydroxide"``, or ``"mu2-hydroxide"`` for a
        bridging one) to a distance in angstrom; by default the input
        geometry's. ``tolerance`` is ``{"distance": 0.05, "angle": 10.0}``.

        :param name:
            The property name, the prefix of the target names.
        :param values:
            The property settings: system, reference, tolerance, weight,
            role, method.

        :return:
            The target settings by target name.
        """

        unknown = set(values) - _TARGET_KEYS
        assert_msg_critical(
            not unknown,
            f"property {name!r} (a metal node) has unsupported keys "
            f"{sorted(unknown)}; use {sorted(_TARGET_KEYS)}")
        self.check()
        tolerance = values.get("tolerance", {})
        if not isinstance(tolerance, dict):
            tolerance = {"distance": tolerance}
        tolerance = {**DEFAULT_TOLERANCES, **tolerance}
        references = dict(values.get("reference") or {})
        common = {k: values[k] for k in ("weight", "role", "method")
                  if k in values}
        metal = self.labels[self.metal]
        xyz = self.molecule.get_coordinates_in_angstrom()

        # metal-donor distances by residue type
        by_type = {}
        for r in self.residues:
            if self.metal not in r.metals:
                continue
            for d in r.donors:
                if self._bonded(d, self.metal):
                    by_type.setdefault((r.label, r.kind, self.elements[d]),
                                       []).append(d)
        targets = {}
        for (label, kind, element), donors in by_type.items():
            reference = references.get(label, references.get(kind))
            if reference is None:
                reference = round(float(np.mean(
                    [np.linalg.norm(xyz[d] - xyz[self.metal])
                     for d in donors])), 4)
            targets[f"{name}: {self.elements[self.metal]}-{element}({label})"] = {
                "system": self.name, "property": "distance",
                "atoms": [[metal, self.labels[d]] for d in donors],
                "value": float(reference), "tolerance": tolerance["distance"],
                "unit": "angstrom", **common}
        leftover = set(references) - {k for key in by_type for k in key[:2]}
        assert_msg_critical(
            not leftover,
            f"property {name!r}: no donor of type {sorted(leftover)} binds "
            f"{metal}; types are {sorted({l for l, _, _ in by_type})}")

        # donor-metal-donor angles, split into cis and trans by the input
        donors = [d for ds in by_type.values() for d in ds]
        cis, trans = [], []
        for i, a in enumerate(donors):
            for b in donors[i + 1:]:
                u, v = xyz[a] - xyz[self.metal], xyz[b] - xyz[self.metal]
                angle = np.degrees(np.arccos(np.clip(
                    u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
                (trans if angle > TRANS_ANGLE_DEGREES else cis).append(
                    [self.labels[a], metal, self.labels[b]])
        elements = {self.elements[d] for d in donors}
        ligand = elements.pop() if len(elements) == 1 else "L"
        shape = f"{ligand}-{self.elements[self.metal]}-{ligand}"
        for kind, triples in (("cis ", cis), ("trans ", trans)):
            if not triples:
                continue
            kind = kind if trans else ""
            targets[f"{name}: {kind}{shape} RMSD"] = {
                "system": self.name, "property": "angle_rmsd",
                "atoms": triples, "value": 0.0,
                "tolerance": tolerance["angle"], "unit": "degree", **common}
        return targets

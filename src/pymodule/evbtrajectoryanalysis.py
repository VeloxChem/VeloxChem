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

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from mpi4py import MPI

from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .reactionsystembuilder import EvbForceGroup

try:
    import openmm as mm
    import openmm.unit as mmunit
except ImportError:
    pass

# OpenMM's Coulomb constant (kJ nm / (mol e^2)), as its NonbondedForce uses it
ONE_4PI_EPS0 = 138.93545764438198

HAMILTONIANS = ("E1", "E2", "E_int")
# The fragment class is split per topology residue ('fragment:A:124:ASP'),
# in the order of the residues; see EvbTrajectoryAnalysis.class_names
ATOM_CLASSES = ("protein", "fragment", "water", "ions", "other")
BONDED_CATEGORIES = ("bond", "angle", "torsion", "other")
# Label of the class-pair columns, and the separator of the two classes
CLASS_PAIR_LABEL = "classes"
CLASS_PAIR_SEPARATOR = "~"

_WATER_RESIDUES = {"HOH", "SOL", "WAT", "TIP3", "TIP3P", "T3P", "SPC", "TIP4"}
_ION_RESIDUES = {
    "ION", "NA", "CL", "K", "NA+", "CL-", "K+", "LI", "RB", "CS", "BR", "IOD",
    "F"
}
# Forces that hold no interaction energy between particles
_IGNORED_FORCE_TYPES = {
    "CMMotionRemover", "MonteCarloBarostat", "MonteCarloAnisotropicBarostat",
    "MonteCarloMembraneBarostat", "MonteCarloFlexibleBarostat",
    "AndersenThermostat"
}
_RESTRAINT_GROUPS = {"CONSTRAINT", "RESTRAINT", "PES_CONSTRAINT", "FROZEN"}
_COULOMB_GROUPS = {"SYSCOUL_STATIC", "SYSCOUL_DYNAMIC"}
_LJ_GROUPS = {"SYSLJ_STATIC", "SYSLJ_DYNAMIC"}
_TERM_CONTAINERS = ("Bonds", "Angles", "Torsions")
_MAX_FORCE_GROUPS = 32


class EvbTrajectoryAnalysis:
    """
    Recomputes, for every frame of a sampled EVB trajectory, the energies of
    chosen force-field terms and residues, in the reactant (E1) and product
    (E2) PES systems and in the integration system of the frame's own lambda
    (E_int).

    Tracked terms are tuples of 2, 3 or 4 atoms. Each gets its geometry
    (distance in nm, angle or dihedral in degrees) and, per Hamiltonian, the
    energy of every bonded term on exactly those atoms ('bonded', restraints
    apart in 'restraint'); an atom pair also gets its nonbonded Coulomb and
    Lennard-Jones energy as the Hamiltonian counts it (exclusions, 1-4
    scaling, soft-core pair forces of the reacting fragment). The Coulomb
    energy of a pair in the NonbondedForce is plain Coulomb at the minimum
    image distance: under PME a single pair has no other definition.

    Tracked residues get, per Hamiltonian:
        - bonded/{bond,angle,torsion,other}: terms with every atom inside the
          residue. Terms shared with another residue (the peptide bond and
          the angles and torsions across it) are left out entirely, so no
          term is counted twice.
        - restraint: restraint terms with every atom inside the residue.
        - nb_internal/{coul,lj}: nonbonded pairs inside the residue, pairwise.
        - nb_rest/{coul,lj}: nonbonded interaction with every atom outside
          the residue, consistent with the simulation's electrostatics:
          E(all) - E(residue switched off) - E(residue alone), evaluated with
          the system's own NonbondedForce (PME included), plus the soft-core
          pair forces that cross the residue boundary.
        - nb_rest/<class>/{coul,lj}: the same interaction split over the
          atom classes (class_names: protein, one fragment class per residue
          holding reacting atoms, water, ions and other), pairwise: plain
          Coulomb at the minimum image of each residue (no cutoff),
          Lennard-Jones with the system's cutoff and switching function. The
          Coulomb class parts therefore do not sum to the PME total in
          nb_rest/coul; the Lennard-Jones parts sum to nb_rest/lj up to the
          long-range dispersion correction.

    The interaction of two tracked residues with each other is part of both
    residues' nb_rest, so nb_rest must not be summed over residues. Instead
    of a residue, an atom group can be tracked the same way: 'fragment' (the
    reacting atoms, which can reach into protein residues) or a dict with
    'label' and 'atoms'. Groups may overlap; for a group spanning several
    residues, 'internal' means inside the group.

    A selected residue that holds reacting atoms and other atoms (a residue
    cut by pdb_active_res 'cut', e.g. DhlA Asp124 CG/OD1/OD2) is tracked as
    two groups, '<residue>[fragment]' (its reacting atoms) and
    '<residue>[rest]' (the others; the cut's host atom belongs here). The
    bonded terms with atoms on both sides go to '<residue>[cut]/<H>/bonded/
    <category>', so no bonded term of the residue is lost.

    With class_pairs, the NonbondedForce energy is also split over pairs of
    atom classes, consistently with the simulation's electrostatics:
    classes/<H>/<a>~<b>/{coul,lj} = E(a and b on) - E(a on) - E(b on), and
    classes/<H>/<a>~<a>/{coul,lj} = E(a on), every other class switched off,
    exceptions between classes added pairwise, the long-range dispersion
    correction split analytically. Summed over all pairs (a <= b) they give
    the NonbondedForce energy exactly. Soft-core pair forces (SYS*) are not
    in these columns; they are fragment-internal and appear in the nb_rest
    of the fragment groups.

    Instance variables:
        - comm: The MPI communicator.
        - ostream: The output stream.
        - topology: The OpenMM topology of the sampled system.
        - force_group_map: Force group name -> force group number of the
          systems, as saved with them (force_group_name.json).
        - hamiltonians: The Hamiltonians evaluated, a subset of E1, E2, E_int.
        - terms: The tracked terms, dicts with label, kind and atoms.
        - residues: The tracked residues and groups, dicts with label, kind
          ('residue', 'fragment part', 'rest part', 'group') and atoms.
        - cuts: The residues split into a fragment and a rest part, dicts
          with label ('<residue>[cut]'), atoms, fragment and rest.
        - fragment_atoms: Particle indices of the reacting fragment.
        - class_names: The atom classes: protein, fragment:<residue> per
          residue holding reacting atoms, water, ions and other.
        - atom_classes: Class index (into class_names) of every particle.
        - class_pairs: Whether the class-pair columns are computed.
        - columns: The output column names, in order.
        - periodic: Whether the frames carry a periodic box.
        - platform: OpenMM platform name, or None for OpenMM's choice.
        - platform_properties: OpenMM platform properties.
    """

    def __init__(self,
                 topology,
                 reference_system,
                 force_group_map=None,
                 terms=None,
                 residues=None,
                 hamiltonians=HAMILTONIANS,
                 platform=None,
                 platform_properties=None,
                 class_pairs=True,
                 comm=None,
                 ostream=None):
        """
        Initializes the analysis and resolves the selections.

        :param topology:
            The OpenMM topology of the sampled system.
        :param reference_system:
            A system of the run (the reactant PES system), used to find the
            reacting fragment.
        :param force_group_map:
            Force group name -> number as the systems were built with. None
            uses the current EvbForceGroup definition.
        :param terms:
            Tracked terms: a list of atom selections, or a dict label -> atom
            selection. An atom selection is a list of 2, 3 or 4 atoms, each a
            particle index, a 'chain:resid:name' string or a (chain, resid,
            name) tuple.
        :param residues:
            Tracked residues: 'chain:resid' strings or (chain, resid) tuples;
            also 'fragment' (the reacting atoms) or an atom group {'label':
            name, 'atoms': [atom selections]}. Groups may overlap. A residue
            holding reacting atoms and other atoms is split (see above).
        :param hamiltonians:
            The Hamiltonians to evaluate, a subset of ('E1', 'E2', 'E_int').
        :param platform:
            OpenMM platform name, or None.
        :param platform_properties:
            OpenMM platform properties, or None.
        :param class_pairs:
            Whether to split the NonbondedForce energy over pairs of atom
            classes (PME-consistent; n(n+1)/2 extra energy evaluations per
            kind for n non-empty classes).
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbTrajectoryAnalysis.')

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == 0:
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.ostream = ostream

        self.topology = topology
        if force_group_map is None:
            force_group_map = {fg.name: fg.value for fg in EvbForceGroup}
        self.force_group_map = dict(force_group_map)
        self._group_names = {v: k for k, v in self.force_group_map.items()}

        hamiltonians = tuple(hamiltonians)
        for ham in hamiltonians:
            assert_msg_critical(
                ham in HAMILTONIANS,
                f"Unknown Hamiltonian '{ham}'; choose from {HAMILTONIANS}")
        assert_msg_critical(
            len(hamiltonians) > 0 and len(set(hamiltonians)) == len(hamiltonians),
            "hamiltonians must name each of E1, E2, E_int at most once")
        self.hamiltonians = hamiltonians

        self.platform = platform
        self.platform_properties = platform_properties

        # From the system: a vacuum topology read from a PDB can carry a cell
        self.periodic = reference_system.usesPeriodicBoundaryConditions()
        self._atoms = list(topology.atoms())
        self.n_atoms = len(self._atoms)
        assert_msg_critical(
            reference_system.getNumParticles() == self.n_atoms,
            f"The topology has {self.n_atoms} atoms but the system has "
            f"{reference_system.getNumParticles()} particles")

        self.class_pairs = bool(class_pairs)
        self.fragment_atoms = self._fragment_atoms(reference_system)
        self.terms = self._resolve_terms(terms)
        self.cuts = []
        self.residues = self._resolve_residues(residues)
        assert_msg_critical(
            len(self.terms) + len(self.residues) > 0 or self.class_pairs,
            "EvbTrajectoryAnalysis needs at least one term or residue to "
            "track, or class_pairs")

        self.class_names, self.atom_classes = self._classify_atoms()

        # Per-residue images: every topology residue is moved as a whole to
        # the image nearest the tracked residue
        self._atom_group = np.array([a.residue.index for a in self._atoms])
        first = {}
        for atom in self._atoms:
            first.setdefault(atom.residue.index, atom.index)
        self._group_ref = np.array(
            [first[g] for g in range(max(first) + 1)] if first else [])
        for residue in self.residues:
            idx = residue["atoms"]
            mask = np.ones(self.n_atoms, dtype=bool)
            mask[idx] = False
            residue["others"] = np.where(mask)[0]
            residue["other_classes"] = self.atom_classes[residue["others"]]

        self.columns = self._make_columns()

    # ------------------------------------------------------------------ #
    # Selections
    # ------------------------------------------------------------------ #

    @staticmethod
    def _atom_label(atom):
        return f"{atom.residue.chain.id}:{atom.residue.id}:{atom.name}"

    def _resolve_atom(self, selection):
        """Particle index of one atom selection."""

        if isinstance(selection, (int, np.integer)):
            index = int(selection)
            assert_msg_critical(
                0 <= index < self.n_atoms,
                f"Atom index {index} is outside the system (0 to "
                f"{self.n_atoms - 1})")
            return index

        if isinstance(selection, str):
            parts = selection.split(":")
        else:
            parts = list(selection)
        assert_msg_critical(
            len(parts) == 3,
            f"Atom selection {selection!r} must be an index, "
            "'chain:resid:name' or (chain, resid, name)")
        chain, resid, name = (str(p).strip() for p in parts)
        matches = [
            a.index for a in self._atoms
            if a.residue.chain.id == chain and str(a.residue.id) == resid
            and a.name == name
        ]
        assert_msg_critical(
            len(matches) == 1,
            f"Atom selection {selection!r} matches {len(matches)} atoms; it "
            "must match exactly one")
        return matches[0]

    def _resolve_residue(self, selection):
        """The topology residue of one residue selection."""

        if isinstance(selection, str):
            parts = selection.split(":")
        else:
            parts = list(selection)
        assert_msg_critical(
            len(parts) == 2,
            f"Residue selection {selection!r} must be 'chain:resid' or "
            "(chain, resid)")
        chain, resid = (str(p).strip() for p in parts)
        matches = [
            r for r in self.topology.residues()
            if r.chain.id == chain and str(r.id) == resid
        ]
        assert_msg_critical(
            len(matches) == 1,
            f"Residue selection {selection!r} matches {len(matches)} "
            "residues; it must match exactly one")
        return matches[0]

    def _resolve_terms(self, terms):
        if terms is None:
            return []
        if isinstance(terms, dict):
            items = list(terms.items())
        else:
            items = [(None, sel) for sel in terms]

        resolved = []
        kinds = {2: "distance", 3: "angle", 4: "dihedral"}
        for label, selection in items:
            atoms = [self._resolve_atom(sel) for sel in selection]
            assert_msg_critical(
                len(atoms) in kinds,
                f"A tracked term needs 2, 3 or 4 atoms, got {selection!r}")
            assert_msg_critical(
                len(set(atoms)) == len(atoms),
                f"Tracked term {selection!r} names an atom twice")
            if label is None:
                label = "-".join(
                    self._atom_label(self._atoms[i]) for i in atoms)
            resolved.append({
                "label": str(label),
                "kind": kinds[len(atoms)],
                "atoms": atoms,
            })
        labels = [t["label"] for t in resolved]
        assert_msg_critical(
            len(set(labels)) == len(labels),
            f"Tracked terms must have unique labels: {labels}")
        return resolved

    def _resolve_residues(self, residues):
        if residues is None:
            return []
        resolved = []
        fragment = set(self.fragment_atoms.tolist())

        def add(label, kind, atoms):
            resolved.append({
                "label": label,
                "kind": kind,
                "atoms": np.array(sorted(atoms), dtype=int),
            })

        for selection in residues:
            if isinstance(selection, str) and selection == "fragment":
                assert_msg_critical(
                    len(fragment) > 0, "The systems have no reacting fragment")
                add("fragment", "group", fragment)
            elif isinstance(selection, dict):
                assert_msg_critical(
                    set(selection) == {"label", "atoms"},
                    "An atom group is a dict with 'label' and 'atoms'")
                label = str(selection["label"])
                atoms = [self._resolve_atom(a) for a in selection["atoms"]]
                assert_msg_critical(
                    len(atoms) > 0 and len(set(atoms)) == len(atoms),
                    f"Atom group {label!r} is empty or names an atom twice")
                add(label, "group", atoms)
            else:
                residue = self._resolve_residue(selection)
                label = f"{residue.chain.id}:{residue.id}:{residue.name}"
                atoms = {a.index for a in residue.atoms()}
                inside = atoms & fragment
                if not inside or inside == atoms:
                    add(label, "residue", atoms)
                    continue
                # A residue cut by the reacting fragment (pdb_active_res
                # 'cut'): its reacting atoms sit in REA_*/SYS* force groups,
                # the rest in the protein's, so it is tracked as two parts
                rest = atoms - inside
                add(f"{label}[fragment]", "fragment part", inside)
                add(f"{label}[rest]", "rest part", rest)
                self.cuts.append({
                    "label": f"{label}[cut]",
                    "residue": label,
                    "atoms": np.array(sorted(atoms), dtype=int),
                    "fragment": np.array(sorted(inside), dtype=int),
                    "rest": np.array(sorted(rest), dtype=int),
                })
        labels = ([r["label"] for r in resolved] +
                  [c["label"] for c in self.cuts])
        assert_msg_critical(
            len(set(labels)) == len(labels),
            f"A residue or group is selected twice: {labels}")
        assert_msg_critical(
            CLASS_PAIR_LABEL not in labels,
            f"'{CLASS_PAIR_LABEL}' is reserved for the class-pair columns")
        return resolved

    def _fragment_atoms(self, system):
        """Particles in the reaction's own force groups (REA_*, SYS*)."""

        groups = {
            value for name, value in self.force_group_map.items()
            if name.startswith("REA_") or name.startswith("SYS")
        }
        atoms = set()
        for force in system.getForces():
            if force.getForceGroup() not in groups:
                continue
            parsed = _parse_force_terms(force)
            if parsed is None:
                continue
            for term_atoms in parsed["atoms"]:
                atoms.update(term_atoms)
        return np.array(sorted(atoms), dtype=int)

    def _classify_atoms(self):
        """The class names and the class index of every particle. The
        reacting atoms get one class per topology residue,
        'fragment:<chain>:<resid>:<resname>', so that e.g. the substrate and
        the reacting part of a cut residue can be told apart."""

        fragment = set(self.fragment_atoms.tolist())
        base = []
        fragment_residues = []
        for residue in self.topology.residues():
            names = {a.name for a in residue.atoms()}
            if residue.name.upper() in _WATER_RESIDUES:
                cls = "water"
            elif residue.name.upper() in _ION_RESIDUES:
                cls = "ions"
            elif {"N", "CA", "C"} <= names:
                cls = "protein"
            else:
                cls = "other"
            for atom in residue.atoms():
                if atom.index in fragment:
                    label = (f"fragment:{residue.chain.id}:{residue.id}:"
                             f"{residue.name}")
                    if label not in fragment_residues:
                        fragment_residues.append(label)
                    base.append((atom.index, label))
                else:
                    base.append((atom.index, cls))

        names = ["protein"] + fragment_residues + ["water", "ions", "other"]
        assert_msg_critical(
            len(set(names)) == len(names),
            f"Atom class names are not unique: {names}")
        position = {name: k for k, name in enumerate(names)}
        classes = np.empty(self.n_atoms, dtype=int)
        for index, name in base:
            classes[index] = position[name]
        return tuple(names), classes

    # ------------------------------------------------------------------ #
    # Columns
    # ------------------------------------------------------------------ #

    def _make_columns(self):
        columns = [f"{t['label']}/{t['kind']}" for t in self.terms]
        for ham in self.hamiltonians:
            for term in self.terms:
                columns += [
                    f"{term['label']}/{ham}/bonded",
                    f"{term['label']}/{ham}/restraint"
                ]
                if term["kind"] == "distance":
                    columns += [
                        f"{term['label']}/{ham}/coul",
                        f"{term['label']}/{ham}/lj"
                    ]
            for residue in self.residues:
                base = f"{residue['label']}/{ham}"
                columns += [f"{base}/bonded/{c}" for c in BONDED_CATEGORIES]
                columns += [f"{base}/restraint"]
                columns += [
                    f"{base}/nb_internal/coul", f"{base}/nb_internal/lj",
                    f"{base}/nb_rest/coul", f"{base}/nb_rest/lj"
                ]
                for cls in self.class_names:
                    columns += [
                        f"{base}/nb_rest/{cls}/coul",
                        f"{base}/nb_rest/{cls}/lj"
                    ]
            for cut in self.cuts:
                columns += [
                    f"{cut['label']}/{ham}/bonded/{c}"
                    for c in BONDED_CATEGORIES
                ]
            if self.class_pairs:
                for a, b in self.class_pair_list():
                    base = f"{CLASS_PAIR_LABEL}/{ham}/{a}{CLASS_PAIR_SEPARATOR}{b}"
                    columns += [f"{base}/coul", f"{base}/lj"]
        return columns

    def class_pair_list(self):
        """The class pairs (a, b), a <= b in class order, of the class-pair
        columns."""

        names = self.class_names
        return [(names[i], names[j])
                for i in range(len(names))
                for j in range(i, len(names))]

    # ------------------------------------------------------------------ #
    # Frames
    # ------------------------------------------------------------------ #

    def frame_geometry(self, positions, box):
        """
        The Hamiltonian-independent part of one frame: tracked geometries and
        the distances from every tracked residue to every other atom.

        :param positions:
            (N, 3) positions in nm.
        :param box:
            (3, 3) box vectors in nm, or None.

        :return:
            A dict with 'values' (column -> value) and per-residue distances.
        """

        positions = np.asarray(positions, dtype=np.float64)
        if not self.periodic:
            box = None
        assert_msg_critical(
            not self.periodic or box is not None,
            "The system is periodic but a frame has no box")

        values = {}
        for term in self.terms:
            x = positions[term["atoms"]]
            column = f"{term['label']}/{term['kind']}"
            if term["kind"] == "distance":
                values[column] = float(
                    np.linalg.norm(_minimum_image(x[1] - x[0], box)))
            elif term["kind"] == "angle":
                values[column] = _angle(x, box)
            else:
                values[column] = _dihedral(x, box)

        residues = []
        for residue in self.residues:
            idx = residue["atoms"]
            others = residue["others"]
            if box is not None:
                ref = positions[self._group_ref]
                d = ref - positions[idx[0]]
                shift = _minimum_image(d, box) - d
                shifted = positions[others] + shift[self._atom_group[others]]
                d_int = _minimum_image(
                    positions[idx][None, :, :] - positions[idx][:, None, :],
                    box)
            else:
                shifted = positions[others]
                d_int = positions[idx][None, :, :] - positions[idx][:, None, :]
            r_cross = np.linalg.norm(shifted[None, :, :] -
                                     positions[idx][:, None, :],
                                     axis=2)
            residues.append({
                "r_cross": r_cross,
                "r_internal": np.linalg.norm(d_int, axis=2),
            })

        pairs = {}
        for term in self.terms:
            if term["kind"] == "distance":
                pairs[term["label"]] = values[f"{term['label']}/distance"]

        return {"values": values, "residues": residues, "pairs": pairs}

    def evaluator(self, system, hamiltonian):
        """
        An evaluator of one Hamiltonian's system. Close it when done.

        :param system:
            The OpenMM system of this Hamiltonian.
        :param hamiltonian:
            'E1', 'E2' or 'E_int': the column label it writes.
        """

        return _SystemEvaluator(self, system, hamiltonian)

    # ------------------------------------------------------------------ #
    # Output
    # ------------------------------------------------------------------ #

    def metadata(self):
        """Description of the selections and the conventions, for json."""

        def atom_info(i):
            atom = self._atoms[i]
            return {
                "index": int(i),
                "chain": atom.residue.chain.id,
                "resid": str(atom.residue.id),
                "resname": atom.residue.name,
                "name": atom.name,
            }

        return {
            "hamiltonians": list(self.hamiltonians),
            "terms": [{
                "label": t["label"],
                "kind": t["kind"],
                "atoms": [atom_info(i) for i in t["atoms"]],
            } for t in self.terms],
            "residues": [{
                "label": r["label"],
                "kind": r["kind"],
                "atoms": [atom_info(i) for i in r["atoms"]],
            } for r in self.residues],
            "cuts": [{
                "label": c["label"],
                "residue": c["residue"],
                "fragment": [atom_info(i) for i in c["fragment"]],
                "rest": [atom_info(i) for i in c["rest"]],
            } for c in self.cuts],
            "atom_classes": {
                cls: int(np.sum(self.atom_classes == k))
                for k, cls in enumerate(self.class_names)
            },
            "class_pairs": self.class_pairs,
            "units": {
                "energy": "kJ/mol",
                "distance": "nm",
                "angle": "degree",
                "dihedral": "degree",
            },
            "conventions": [
                "E1/E2: reactant/product PES systems; E_int: the integration "
                "system of the frame's own lambda.",
                "bonded/<category>: terms with every atom in the residue; "
                "terms shared with another residue are not counted.",
                "nb_rest/coul and nb_rest/lj: E(all) - E(residue off) - "
                "E(residue alone) with the system's NonbondedForce (PME "
                "included), plus the soft-core pair forces crossing the "
                "residue boundary.",
                "nb_rest/<class>: pairwise, plain Coulomb at the minimum "
                "image of each residue without cutoff, Lennard-Jones with "
                "the system's cutoff and switching; no dispersion correction.",
                "nb_internal: pairwise inside the residue, as nb_rest/<class>.",
                "A pair term's coul/lj: the pair as the Hamiltonian counts it "
                "(exceptions, soft-core pair forces), plain Coulomb at the "
                "minimum image distance.",
                "nb_rest of two tracked residues both contain their mutual "
                "interaction: do not sum nb_rest over residues.",
                "Atom classes: protein, fragment:<chain>:<resid>:<resname> "
                "(the reacting atoms of that residue), water, ions, other.",
                "A residue holding reacting atoms and other atoms is tracked "
                "as <residue>[fragment] and <residue>[rest]; "
                "<residue>[cut]/<H>/bonded/<category> holds the bonded terms "
                "with atoms on both sides (every atom in the residue).",
                f"{CLASS_PAIR_LABEL}/<H>/<a>{CLASS_PAIR_SEPARATOR}<b>/"
                "{coul,lj}: NonbondedForce energy between atom classes a and "
                "b (a before b in atom_classes), E(a, b on) - E(a on) - "
                f"E(b on); <a>{CLASS_PAIR_SEPARATOR}<a>: E(a on) alone. PME "
                "included, exceptions between classes pairwise, long-range "
                "dispersion correction split analytically. The sum over all "
                "pairs is the NonbondedForce energy; soft-core pair forces "
                "are not included.",
            ],
        }

    @staticmethod
    def read_csv(path, rows=None, n_rows=None):
        """
        Reads a TrajectoryAnalysis.csv.

        :param path:
            The csv file.
        :param rows:
            Energies.csv row indices to return, in order (None: all analysed
            rows). A row that was not analysed (frame_stride) comes back NaN.
        :param n_rows:
            The number of rows of Energies.csv, checked against the file.

        :return:
            Dict with 'names' (comma-joined column names), 'values' (rows x
            columns), 'frame', 'file' and 'metadata_file'.
        """

        with Path(path).open("r") as file:
            header = [cell.strip() for cell in file.readline().split(",")]
        data = np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2)
        assert_msg_critical(
            data.shape[1] == len(header),
            f"{path} has {data.shape[1]} columns but its header names "
            f"{len(header)}")
        assert_msg_critical(header[:4] == ["frame", "Lambda", "replica",
                                           "direction"],
                            f"{path} does not start with frame, Lambda, "
                            "replica, direction")
        frames = data[:, 0].astype(int)
        names = header[4:]
        values = data[:, 4:]
        metadata = Path(path).with_suffix(".json")
        if n_rows is not None and metadata.is_file():
            with metadata.open("r") as file:
                expected = json.load(file).get("n_frames")
            assert_msg_critical(
                expected is None or expected == n_rows,
                f"{path} was computed for {expected} frames, but "
                f"Energies.csv has {n_rows} rows")
        if rows is None:
            rows = frames
        rows = np.asarray(rows, dtype=int)
        out = np.full((len(rows), len(names)), np.nan)
        position = {int(f): k for k, f in enumerate(frames)}
        for k, row in enumerate(rows):
            if int(row) in position:
                out[k] = values[position[int(row)]]
        return {
            "names": ",".join(names),
            "values": out,
            "frame": rows,
            "file": str(Path(path).resolve()),
            "metadata_file": str(metadata.resolve()),
        }

    @staticmethod
    def column(analysis, name):
        """One column of a loaded analysis (result['trajectory_analysis'])."""

        names = analysis["names"].split(",")
        assert_msg_critical(name in names,
                            f"No column {name!r} in the trajectory analysis")
        return analysis["values"][:, names.index(name)]

    @staticmethod
    def lambda_averages(result):
        """
        Per-lambda averages of a loaded analysis.

        :param result:
            One configuration's results (with 'trajectory_analysis', and
            'Lambda_frame' from the analysis itself).

        :return:
            Dict with 'Lambda', 'names' and, per column, 'mean', 'sem' and
            'n' arrays (lambda x columns). The SEM is over the (replica,
            direction) block means and is NaN with fewer than 2 blocks;
            frames within a block are correlated, so a frame-level SEM would
            be too small.
        """

        analysis = result["trajectory_analysis"]
        lam = np.round(np.asarray(analysis["Lambda_frame"], dtype=float), 3)
        blocks = (np.asarray(analysis["Replica_frame"]) * 2 +
                  np.asarray(analysis["Direction_frame"]))
        values = analysis["values"]
        lambdas = np.unique(lam)
        n_cols = values.shape[1]
        mean = np.full((len(lambdas), n_cols), np.nan)
        sem = np.full((len(lambdas), n_cols), np.nan)
        count = np.zeros((len(lambdas), n_cols), dtype=int)
        for k, L in enumerate(lambdas):
            sel = lam == L
            v = values[sel]
            finite = np.isfinite(v)
            count[k] = finite.sum(axis=0)
            with np.errstate(invalid="ignore"):
                mean[k] = np.nanmean(np.where(finite, v, np.nan), axis=0)
            block_means = []
            for b in np.unique(blocks[sel]):
                vb = v[blocks[sel] == b]
                if np.isfinite(vb).any():
                    with np.errstate(invalid="ignore"):
                        block_means.append(np.nanmean(vb, axis=0))
            if len(block_means) >= 2:
                bm = np.array(block_means)
                n_b = np.isfinite(bm).sum(axis=0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    sem[k] = np.nanstd(bm, axis=0, ddof=1) / np.sqrt(n_b)
        return {
            "Lambda": lambdas,
            "names": analysis["names"],
            "mean": mean,
            "sem": sem,
            "n": count,
        }


class _SystemEvaluator:
    """
    Evaluates the tracked quantities of one Hamiltonian's system for frames
    handed to it. Holds OpenMM contexts until close() is called.
    """

    def __init__(self, analysis, system, hamiltonian):

        self.analysis = analysis
        self.hamiltonian = hamiltonian
        self.matched = {}
        self._contexts = []
        self._group_columns = []
        self._nb_layers = []
        self._class_layer = None

        assert_msg_critical(
            system.getNumParticles() == analysis.n_atoms,
            f"The {hamiltonian} system has {system.getNumParticles()} "
            f"particles, the topology {analysis.n_atoms}")

        self._read_nonbonded(system)
        self._build_term_contexts(system)
        if analysis.residues and self._nb is not None:
            self._build_nonbonded_contexts(system)
        if analysis.class_pairs and self._nb is not None:
            self._build_class_contexts(system)

    # --- setup ---------------------------------------------------------- #

    def _read_nonbonded(self, system):
        forces = [
            f for f in system.getForces() if isinstance(f, mm.NonbondedForce)
        ]
        assert_msg_critical(
            len(forces) <= 1,
            f"The {self.hamiltonian} system has {len(forces)} "
            "NonbondedForces; at most one is supported")
        if not forces:
            self._nb = None
            return
        nb = forces[0]
        assert_msg_critical(
            nb.getNumParticleParameterOffsets() == 0
            and nb.getNumExceptionParameterOffsets() == 0,
            "NonbondedForce parameter offsets are not supported")
        self._nb = nb

        # Parameters from the XML (in nm, e, kJ/mol): reading them one
        # particle at a time through OpenMM's unit objects is slow
        self._nb_xml = mm.XmlSerializer.serialize(nb)
        root = ET.fromstring(self._nb_xml)
        q, sig, eps = _xml_parameters(root.find("Particles")).T
        self.q, self.sig, self.eps = q, sig, eps
        n = len(q)
        exceptions = root.find("Exceptions")
        ex = np.array([(int(e.get("p1")), int(e.get("p2")))
                       for e in exceptions],
                      dtype=int).reshape(-1, 2)
        ex_par = _xml_parameters(exceptions)
        self._ex = ex
        self._ex_par = ex_par

        self._method = int(root.get("method"))
        if nb.getNonbondedMethod() == mm.NonbondedForce.NoCutoff:
            self.cutoff = np.inf
        else:
            self.cutoff = float(root.get("cutoff"))
        if root.get("useSwitchingFunction") == "1" and np.isfinite(
                self.cutoff):
            self.switch = float(root.get("switchingDistance"))
        else:
            self.switch = None

        # Exceptions per tracked pair and per tracked residue
        self._pair_exceptions = {}
        pair_keys = {
            tuple(sorted(t["atoms"])): t["label"]
            for t in self.analysis.terms if t["kind"] == "distance"
        }
        for k, (p1, p2) in enumerate(ex):
            key = (min(p1, p2), max(p1, p2))
            if key in pair_keys:
                self._pair_exceptions[pair_keys[key]] = ex_par[k]

        self._residue_exceptions = []
        for residue in self.analysis.residues:
            idx = residue["atoms"]
            in_r = np.zeros(n, dtype=bool)
            in_r[idx] = True
            a, b = in_r[ex[:, 0]], in_r[ex[:, 1]]
            internal = {}
            for k in np.where(a & b)[0]:
                i, j = ex[k]
                internal[(min(i, j), max(i, j))] = ex_par[k]
            # Pairs inside the residue: the plain combination, or the
            # exception where there is one
            ia, ib = np.triu_indices(len(idx), 1)
            i_atoms, j_atoms = idx[ia], idx[ib]
            int_par = np.column_stack([
                q[i_atoms] * q[j_atoms], 0.5 * (sig[i_atoms] + sig[j_atoms]),
                np.sqrt(eps[i_atoms] * eps[j_atoms])
            ])
            int_is_ex = np.zeros(len(ia), dtype=bool)
            for p, (i, j) in enumerate(zip(i_atoms, j_atoms)):
                par = internal.get((min(i, j), max(i, j)))
                if par is not None:
                    int_par[p] = par
                    int_is_ex[p] = True
            cross = np.where(a ^ b)[0]
            inside = np.where(a[cross], ex[cross, 0], ex[cross, 1])
            outside = np.where(a[cross], ex[cross, 1], ex[cross, 0])
            local = np.searchsorted(idx, inside)
            other_pos = np.searchsorted(residue["others"], outside)
            self._residue_exceptions.append({
                "ia": ia,
                "ib": ib,
                "int_par": int_par,
                "int_is_ex": int_is_ex,
                "local": local,
                "other_pos": other_pos,
                "par": ex_par[cross],
            })

    def _role(self, force):
        name = self.analysis._group_names.get(force.getForceGroup(), "")
        if name in _COULOMB_GROUPS:
            return "coul"
        if name in _LJ_GROUPS:
            return "lj"
        if name in _RESTRAINT_GROUPS:
            return "restraint"
        return "bonded"

    def _build_term_contexts(self, system):
        """One filtered copy of every force per output column, each column
        in its own force group, packed into as few contexts as needed."""

        analysis = self.analysis
        ham = self.hamiltonian
        # column -> list of (force index, term indices)
        selections = {}
        skipped = []

        def add(column, f_index, term_index):
            selections.setdefault(column, {}).setdefault(f_index,
                                                         []).append(term_index)

        residue_sets = [set(r["atoms"].tolist()) for r in analysis.residues]
        cut_sets = [(set(c["atoms"].tolist()), set(c["fragment"].tolist()))
                    for c in analysis.cuts]
        atom_class = analysis.atom_classes

        parsed_forces = {}
        for f_index, force in enumerate(system.getForces()):
            ftype = force.__class__.__name__
            if isinstance(force, mm.NonbondedForce) or (ftype
                                                        in _IGNORED_FORCE_TYPES):
                continue
            parsed = _parse_force_terms(force)
            if parsed is None:
                skipped.append(f"{force.getName()} ({ftype})")
                continue
            parsed_forces[f_index] = parsed
            role = self._role(force)
            nb_role = role in ("coul", "lj")

            for t_index, atoms in enumerate(parsed["atoms"]):
                atom_set = set(atoms)
                for term in analysis.terms:
                    target = term["atoms"]
                    if len(atoms) != len(target):
                        continue
                    if not _same_term(atoms, target):
                        continue
                    if nb_role:
                        if len(target) == 2:
                            add(f"{term['label']}/{ham}/{role}", f_index,
                                t_index)
                    else:
                        add(f"{term['label']}/{ham}/{role}", f_index, t_index)

                for residue, rset in zip(analysis.residues, residue_sets):
                    base = f"{residue['label']}/{ham}"
                    n_in = len(atom_set & rset)
                    if n_in == 0:
                        continue
                    if nb_role:
                        if n_in == len(atom_set):
                            add(f"{base}/nb_internal/{role}", f_index, t_index)
                        else:
                            partner = next(a for a in atoms if a not in rset)
                            cls = analysis.class_names[atom_class[partner]]
                            add(f"{base}/nb_rest/{role}", f_index, t_index)
                            add(f"{base}/nb_rest/{cls}/{role}", f_index,
                                t_index)
                    elif n_in == len(atom_set):
                        if role == "restraint":
                            add(f"{base}/restraint", f_index, t_index)
                        else:
                            category = {
                                2: "bond",
                                3: "angle",
                                4: "torsion"
                            }.get(len(atoms), "other")
                            add(f"{base}/bonded/{category}", f_index, t_index)

                if nb_role or role == "restraint":
                    continue
                for cut, (cut_set, frag_set) in zip(analysis.cuts, cut_sets):
                    if (atom_set <= cut_set and atom_set & frag_set
                            and atom_set - frag_set):
                        category = {
                            2: "bond",
                            3: "angle",
                            4: "torsion"
                        }.get(len(atoms), "other")
                        add(f"{cut['label']}/{ham}/bonded/{category}",
                            f_index, t_index)

        if skipped:
            analysis.ostream.print_warning(
                f"{ham}: forces not decomposed (no per-term atoms): " +
                ", ".join(skipped))

        # Record what the tracked terms matched
        for term in analysis.terms:
            for role in ("bonded", "restraint", "coul", "lj"):
                column = f"{term['label']}/{ham}/{role}"
                for f_index, t_indices in selections.get(column, {}).items():
                    parsed = parsed_forces[f_index]
                    for t_index in t_indices:
                        attrib = dict(parsed["terms"][t_index].attrib)
                        self.matched.setdefault(column, []).append({
                            "force": parsed["name"],
                            "type": parsed["type"],
                            "atoms": list(parsed["atoms"][t_index]),
                            "parameters": {
                                k: v for k, v in attrib.items()
                                if not re.fullmatch(r"p\d+", k)
                            },
                        })

        # Pack the columns into contexts, one force group per column
        columns = list(selections)
        for start in range(0, len(columns), _MAX_FORCE_GROUPS):
            chunk = columns[start:start + _MAX_FORCE_GROUPS]
            sub = _empty_copy(system)
            group_columns = []
            for group, column in enumerate(chunk):
                for f_index, t_indices in selections[column].items():
                    sub.addForce(
                        _filtered_force(parsed_forces[f_index], t_indices,
                                        group))
                group_columns.append((group, column))
            self._contexts.append(self._context(sub))
            self._group_columns.append(group_columns)

    def _build_nonbonded_contexts(self, system):
        """Two copies of the NonbondedForce, Coulomb only and Lennard-Jones
        only, in which every tracked residue and the rest of the system are
        switched on and off by a global parameter each (parameter offsets).

        An exception inside one group follows that group's switch. An
        exception between two groups would need the product of two switches,
        so its direct term is left out here (the pair stays excluded, which
        keeps the PME exclusion correction) and added pairwise in evaluate().
        """

        self._dispersion = self._dispersion_cross_coefficients()

        # A particle follows one switch, so overlapping groups (the fragment
        # and a residue it reaches into) go into separate layers, each with
        # its own pair of contexts
        layers = []
        for k, residue in enumerate(self.analysis.residues):
            atoms = set(residue["atoms"].tolist())
            for layer in layers:
                if not atoms & layer["atoms"]:
                    layer["members"].append(k)
                    layer["atoms"] |= atoms
                    break
            else:
                layers.append({"members": [k], "atoms": atoms})

        self._nb_layers = []
        for n_layer, layer in enumerate(layers):
            group = np.zeros(len(self.q), dtype=int)
            for g, k in enumerate(layer["members"]):
                group[self.analysis.residues[k]["atoms"]] = g + 1
            names = [
                f"evb_analysis_switch_{n_layer}_{g}"
                for g in range(len(layer["members"]) + 1)
            ]
            self._nb_layers.append({
                "members": layer["members"],
                "switches": names,
                "contexts": self._switch_contexts(system, group, names),
            })

    def _switch_contexts(self, system, group, names):
        """The Coulomb and Lennard-Jones contexts of one layer: group[i] is
        the switch (index into names) of particle i."""

        contexts = {}
        ex_group = group[self._ex[:, 0]]
        same = ex_group == group[self._ex[:, 1]]

        for kind in ("coul", "lj"):
            nb = mm.XmlSerializer.deserialize(self._nb_xml)
            nb.setForceGroup(0)
            if kind == "lj":
                # Without charges PME only adds work; Lennard-Jones is the
                # same with every periodic cutoff method
                if self._method in (mm.NonbondedForce.PME,
                                    mm.NonbondedForce.Ewald):
                    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            # OpenMM computes the dispersion correction once, from the
            # parameters at context creation, not when a switch changes;
            # its residue-rest part is added in evaluate() instead
            nb.setUseDispersionCorrection(False)
            for name in names:
                nb.addGlobalParameter(name, 1.0)

            # plain floats: the API is fast without unit objects
            values = self.q if kind == "coul" else self.eps
            for i, sig in enumerate(self.sig.tolist()):
                nb.setParticleParameters(i, 0.0, sig, 0.0)
            for i in np.nonzero(values)[0].tolist():
                value = float(values[i])
                if kind == "coul":
                    nb.addParticleParameterOffset(names[group[i]], i, value,
                                                  0.0, 0.0)
                else:
                    nb.addParticleParameterOffset(names[group[i]], i, 0.0,
                                                  0.0, value)

            ex_values = self._ex_par[:, 0 if kind == "coul" else 2]
            nonzero = (self._ex_par[:, 0] != 0.0) | (self._ex_par[:, 2] != 0.0)
            for k in np.nonzero(nonzero)[0].tolist():
                p1, p2 = self._ex[k].tolist()
                nb.setExceptionParameters(k, p1, p2, 0.0,
                                          float(self._ex_par[k, 1]), 0.0)
            for k in np.nonzero(same & (ex_values != 0.0))[0].tolist():
                value = float(ex_values[k])
                if kind == "coul":
                    nb.addExceptionParameterOffset(names[ex_group[k]], k,
                                                   value, 0.0, 0.0)
                else:
                    nb.addExceptionParameterOffset(names[ex_group[k]], k, 0.0,
                                                   0.0, value)

            sub = _empty_copy(system)
            sub.addForce(nb)
            contexts[kind] = self._context(sub)
        return contexts

    def _dispersion_cross_coefficients(self):
        """
        Per tracked residue, the residue-rest part of OpenMM's long-range
        dispersion correction times the box volume (kJ/mol nm^3), or zeros
        when the force has none. Follows NonbondedForceImpl::
        calcDispersionCorrection: 8 pi N^2 / (N (N + 1) / 2) times the sum
        over particle pairs of eps sig^12 / (9 rc^9) - eps sig^6 / (3 rc^3)
        + eps X(sig), X the integral of the switched-off part of the
        Lennard-Jones potential, here by Gauss-Legendre quadrature.
        """

        n_res = len(self.analysis.residues)
        if not self._has_dispersion():
            return np.zeros(n_res)

        coefficients = np.zeros(n_res)
        for k, residue in enumerate(self.analysis.residues):
            idx = residue["atoms"]
            types, counts = self._lj_types(residue["others"])
            g = self._dispersion_g(self.sig[idx][:, None],
                                   types[None, :, 0],
                                   self.eps[idx][:, None],
                                   types[None, :, 1])
            coefficients[k] = np.sum(g * counts[None, :])
        return self._dispersion_prefactor() * coefficients

    def _has_dispersion(self):
        return bool(self._nb.getUseDispersionCorrection() and self._method
                    in (mm.NonbondedForce.CutoffPeriodic,
                        mm.NonbondedForce.Ewald, mm.NonbondedForce.PME))

    def _dispersion_prefactor(self):
        # 8 pi N^2 over OpenMM's number of interactions N (N + 1) / 2
        n = len(self.q)
        return 8.0 * np.pi * n * n / (n * (n + 1) / 2.0)

    def _lj_types(self, atoms):
        """The distinct (sig, eps) of the atoms, and their counts."""
        return np.unique(np.column_stack([self.sig[atoms], self.eps[atoms]]),
                         axis=0,
                         return_counts=True)

    def _dispersion_g(self, sig_a, sig_b, eps_a, eps_b):
        """Per pair of Lennard-Jones types (broadcast), the pair term of
        NonbondedForceImpl::calcDispersionCorrection before the prefactor:
        eps sig^12 / (9 rc^9) - eps sig^6 / (3 rc^3) + eps X(sig)."""

        rc = self.cutoff
        sig = 0.5 * (sig_a + sig_b)
        eps = np.sqrt(eps_a * eps_b)
        s6 = sig**6
        g = eps * (s6 * s6 / (9.0 * rc**9) - s6 / (3.0 * rc**3))
        if self.switch is not None:
            nodes, weights = np.polynomial.legendre.leggauss(64)
            rs = self.switch
            r = 0.5 * (rc - rs) * nodes + 0.5 * (rc + rs)
            w = 0.5 * (rc - rs) * weights
            x = (r - rs) / (rc - rs)
            off = 1.0 - (1.0 - 6.0 * x**5 + 15.0 * x**4 - 10.0 * x**3)
            lj = (sig[..., None] / r)**12 - (sig[..., None] / r)**6
            g = g + eps * np.sum(w * lj * off * r * r, axis=-1)
        return g

    def _dispersion_class_coefficients(self):
        """
        Per class pair (a <= b, as class_pair_list), its part of OpenMM's
        long-range dispersion correction times the box volume. OpenMM counts
        n (n + 1) / 2 pairs within a particle type and n m between two
        types; split over the classes this gives, for a != b, the sum over
        i in a, j in b, and for a == b the sum over i <= j in a (i == j
        included). The parts add up to OpenMM's correction.
        """

        analysis = self.analysis
        n_cls = len(analysis.class_names)
        coefficients = np.zeros((n_cls, n_cls))
        if not self._has_dispersion():
            return coefficients
        types = []
        for c in range(n_cls):
            atoms = np.where(analysis.atom_classes == c)[0]
            types.append(self._lj_types(atoms) if len(atoms) else None)
        for a in range(n_cls):
            if types[a] is None:
                continue
            ta, na = types[a]
            for b in range(a, n_cls):
                if types[b] is None:
                    continue
                tb, nb = types[b]
                g = self._dispersion_g(ta[:, None, 0], tb[None, :, 0],
                                       ta[:, None, 1], tb[None, :, 1])
                total = np.sum(g * na[:, None] * nb[None, :])
                if a == b:
                    total = 0.5 * (total + np.sum(np.diag(g) * na))
                coefficients[a, b] = total
        return self._dispersion_prefactor() * coefficients

    def _build_class_contexts(self, system):
        """The Coulomb and Lennard-Jones contexts with one switch per atom
        class, for the class-pair columns."""

        analysis = self.analysis
        group = analysis.atom_classes
        names = [
            f"evb_analysis_class_{c}" for c in range(len(analysis.class_names))
        ]
        present = [
            c for c in range(len(names)) if np.any(group == c)
        ]
        # Exceptions between two classes: direct term computed pairwise
        ex_a, ex_b = group[self._ex[:, 0]], group[self._ex[:, 1]]
        nonzero = (self._ex_par[:, 0] != 0.0) | (self._ex_par[:, 2] != 0.0)
        cross = np.where((ex_a != ex_b) & nonzero)[0]
        self._class_layer = {
            "switches": names,
            "present": present,
            "contexts": self._switch_contexts(system, group, names),
            "dispersion": self._dispersion_class_coefficients(),
            "cross_ex": cross,
            "cross_a": np.minimum(ex_a[cross], ex_b[cross]),
            "cross_b": np.maximum(ex_a[cross], ex_b[cross]),
        }

    def _context(self, system):
        integrator = mm.VerletIntegrator(0.001)
        platform = self.analysis.platform
        if platform is None:
            return mm.Context(system, integrator)
        platform = mm.Platform.getPlatformByName(platform)
        if self.analysis.platform_properties:
            return mm.Context(system, integrator, platform,
                              self.analysis.platform_properties)
        return mm.Context(system, integrator, platform)

    # --- evaluation ----------------------------------------------------- #

    @staticmethod
    def _set_frame(context, positions, box):
        if box is not None:
            context.setPeriodicBoxVectors(mm.Vec3(*box[0]), mm.Vec3(*box[1]),
                                          mm.Vec3(*box[2]))
        context.setPositions(positions)

    @staticmethod
    def _energy(context, groups=-1):
        state = context.getState(getEnergy=True, groups=groups)
        return state.getPotentialEnergy().value_in_unit(
            mmunit.kilojoule_per_mole)

    def _lj(self, sig, eps, r, cutoff=True):
        """Lennard-Jones of pairs, with the force's cutoff and switching
        (cutoff=False for exceptions, which have neither)."""

        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            s6 = (sig / r)**6
            e = 4.0 * eps * (s6 * s6 - s6)
        e = np.where(eps == 0.0, 0.0, e)
        if not cutoff or not np.isfinite(self.cutoff):
            return e
        e = np.where(r < self.cutoff, e, 0.0)
        if self.switch is not None:
            x = np.clip((r - self.switch) / (self.cutoff - self.switch), 0.0,
                        1.0)
            e = e * (1.0 - 6.0 * x**5 + 15.0 * x**4 - 10.0 * x**3)
        return e

    def evaluate(self, frame, positions, box, values):
        """
        Adds this Hamiltonian's columns for one frame to values.

        :param frame:
            The frame_geometry() of the frame.
        :param positions:
            (N, 3) positions in nm.
        :param box:
            (3, 3) box vectors in nm, or None.
        :param values:
            Dict column -> value to add to.
        """

        analysis = self.analysis
        ham = self.hamiltonian
        if not analysis.periodic:
            box = None

        for column in analysis.columns:
            if f"/{ham}/" in column:
                values[column] = 0.0

        for context, group_columns in zip(self._contexts, self._group_columns):
            self._set_frame(context, positions, box)
            for group, column in group_columns:
                values[column] += self._energy(context, 1 << group)

        if self._nb is None:
            return

        # Tracked pairs in the NonbondedForce
        for term in analysis.terms:
            if term["kind"] != "distance":
                continue
            i, j = term["atoms"]
            r = frame["pairs"][term["label"]]
            ex = self._pair_exceptions.get(term["label"])
            if ex is not None:
                coul = ONE_4PI_EPS0 * ex[0] / r
                lj = float(self._lj(ex[1], ex[2], r, cutoff=False))
            else:
                coul = ONE_4PI_EPS0 * self.q[i] * self.q[j] / r
                lj = float(
                    self._lj(0.5 * (self.sig[i] + self.sig[j]),
                             np.sqrt(self.eps[i] * self.eps[j]), r))
            values[f"{term['label']}/{ham}/coul"] += coul
            values[f"{term['label']}/{ham}/lj"] += lj

        if self._class_layer is not None:
            self._evaluate_class_pairs(positions, box, values)

        if not analysis.residues:
            return

        # Residue totals with the system's own electrostatics:
        # E(all) - E(residue off) - E(residue alone)
        n_res = len(analysis.residues)
        totals = np.zeros((n_res, 2))
        for layer in self._nb_layers:
            switches = layer["switches"]
            for c, kind in enumerate(("coul", "lj")):
                context = layer["contexts"][kind]
                self._set_frame(context, positions, box)
                for name in switches:
                    context.setParameter(name, 1.0)
                e_all = self._energy(context)
                for g, k in enumerate(layer["members"]):
                    context.setParameter(switches[g + 1], 0.0)
                    totals[k, c] += e_all - self._energy(context)
                    context.setParameter(switches[g + 1], 1.0)
                for name in switches:
                    context.setParameter(name, 0.0)
                for g, k in enumerate(layer["members"]):
                    context.setParameter(switches[g + 1], 1.0)
                    totals[k, c] -= self._energy(context)
                    context.setParameter(switches[g + 1], 0.0)
                for name in switches:
                    context.setParameter(name, 1.0)

        n_classes = len(analysis.class_names)
        for k, residue in enumerate(analysis.residues):
            base = f"{residue['label']}/{ham}"

            idx = residue["atoms"]
            others = residue["others"]
            classes = residue["other_classes"]
            geo = frame["residues"][k]
            r = geo["r_cross"]
            q_r, q_o = self.q[idx], self.q[others]

            # Coulomb per class, no cutoff
            per_atom = ONE_4PI_EPS0 * (q_r @ (q_o[None, :] / r))
            coul = np.bincount(classes, weights=per_atom, minlength=n_classes)

            # Lennard-Jones per class within the cutoff
            if np.isfinite(self.cutoff):
                ii, jj = np.nonzero(r < self.cutoff)
            else:
                ii, jj = np.indices(r.shape).reshape(2, -1)
            o_atoms = others[jj]
            lj_pairs = self._lj(0.5 * (self.sig[idx[ii]] + self.sig[o_atoms]),
                                np.sqrt(self.eps[idx[ii]] * self.eps[o_atoms]),
                                r[ii, jj])
            lj = np.bincount(classes[jj], weights=lj_pairs,
                             minlength=n_classes)

            # Exceptions replace the plain pair
            ex = self._residue_exceptions[k]
            if len(ex["local"]):
                li, oj = ex["local"], ex["other_pos"]
                r_ex = r[li, oj]
                a_i, a_j = idx[li], others[oj]
                plain_c = ONE_4PI_EPS0 * self.q[a_i] * self.q[a_j] / r_ex
                plain_lj = self._lj(0.5 * (self.sig[a_i] + self.sig[a_j]),
                                    np.sqrt(self.eps[a_i] * self.eps[a_j]),
                                    r_ex)
                par = ex["par"]
                ex_c = ONE_4PI_EPS0 * par[:, 0] / r_ex
                ex_lj = self._lj(par[:, 1], par[:, 2], r_ex, cutoff=False)
                coul += np.bincount(classes[oj],
                                    weights=ex_c - plain_c,
                                    minlength=n_classes)
                lj += np.bincount(classes[oj],
                                  weights=ex_lj - plain_lj,
                                  minlength=n_classes)
                # their direct term is not in the switched contexts
                totals[k] += (np.sum(ex_c), np.sum(ex_lj))

            values[f"{base}/nb_rest/coul"] += totals[k, 0]
            if box is not None:
                totals[k, 1] += self._dispersion[k] / (box[0, 0] * box[1, 1] *
                                                       box[2, 2])
            values[f"{base}/nb_rest/lj"] += totals[k, 1]

            for c, cls in enumerate(analysis.class_names):
                values[f"{base}/nb_rest/{cls}/coul"] += coul[c]
                values[f"{base}/nb_rest/{cls}/lj"] += lj[c]

            # Inside the residue; exceptions have no cutoff
            d = geo["r_internal"][ex["ia"], ex["ib"]]
            par, is_ex = ex["int_par"], ex["int_is_ex"]
            lj_int = np.where(
                is_ex, self._lj(par[:, 1], par[:, 2], d, cutoff=False),
                self._lj(par[:, 1], par[:, 2], d))
            values[f"{base}/nb_internal/coul"] += float(
                np.sum(ONE_4PI_EPS0 * par[:, 0] / d))
            values[f"{base}/nb_internal/lj"] += float(np.sum(lj_int))

    def _evaluate_class_pairs(self, positions, box, values):
        """Adds the class-pair columns: E(a, b on) - E(a on) - E(b on) and
        E(a on) in the class-switched contexts, plus the exceptions between
        classes and the split dispersion correction."""

        analysis = self.analysis
        layer = self._class_layer
        switches, present = layer["switches"], layer["present"]
        names = analysis.class_names
        n_cls = len(names)
        energy = np.zeros((2, n_cls, n_cls))

        for c, kind in enumerate(("coul", "lj")):
            context = layer["contexts"][kind]
            self._set_frame(context, positions, box)
            for name in switches:
                context.setParameter(name, 0.0)
            alone = {}
            for a in present:
                context.setParameter(switches[a], 1.0)
                alone[a] = self._energy(context)
                context.setParameter(switches[a], 0.0)
                energy[c, a, a] = alone[a]
            for i, a in enumerate(present):
                context.setParameter(switches[a], 1.0)
                for b in present[i + 1:]:
                    context.setParameter(switches[b], 1.0)
                    energy[c, a, b] = (self._energy(context) - alone[a] -
                                       alone[b])
                    context.setParameter(switches[b], 0.0)
                context.setParameter(switches[a], 0.0)

        # Exceptions between classes: their direct term is not in the
        # switched contexts (the pair stays excluded there)
        cross = layer["cross_ex"]
        if len(cross):
            p1, p2 = self._ex[cross, 0], self._ex[cross, 1]
            r = np.linalg.norm(_minimum_image(positions[p2] - positions[p1],
                                              box),
                               axis=1)
            par = self._ex_par[cross]
            ex_c = ONE_4PI_EPS0 * par[:, 0] / r
            ex_lj = self._lj(par[:, 1], par[:, 2], r, cutoff=False)
            pair = layer["cross_a"] * n_cls + layer["cross_b"]
            energy[0] += np.bincount(pair, weights=ex_c,
                                     minlength=n_cls * n_cls).reshape(
                                         n_cls, n_cls)
            energy[1] += np.bincount(pair, weights=ex_lj,
                                     minlength=n_cls * n_cls).reshape(
                                         n_cls, n_cls)

        if box is not None:
            energy[1] += layer["dispersion"] / (box[0, 0] * box[1, 1] *
                                                box[2, 2])

        ham = self.hamiltonian
        position = {name: k for k, name in enumerate(names)}
        for a, b in analysis.class_pair_list():
            base = (f"{CLASS_PAIR_LABEL}/{ham}/{a}{CLASS_PAIR_SEPARATOR}{b}")
            values[f"{base}/coul"] += energy[0, position[a], position[b]]
            values[f"{base}/lj"] += energy[1, position[a], position[b]]

    def close(self):
        """Releases the OpenMM contexts."""

        self._contexts = []
        self._group_columns = []
        self._nb_layers = []
        self._class_layer = None


# ---------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------- #


_PARSE_CACHE = {}


def _parse_force_terms(force):
    """
    The per-term structure of a force from its XML: the root element, the
    term container, the term elements and the atoms of each term. None for a
    force whose terms are not defined by particle indices.

    Cached by the serialized force: the lambda systems of a run share most
    of their forces (the whole protein), and parsing those is the bulk of
    the setup. The cache holds the last 64 forces.
    """

    xml = mm.XmlSerializer.serialize(force)
    if xml in _PARSE_CACHE:
        return _PARSE_CACHE[xml]
    parsed = _parse_force_xml(xml, force.getName())
    if len(_PARSE_CACHE) >= 64:
        _PARSE_CACHE.pop(next(iter(_PARSE_CACHE)))
    _PARSE_CACHE[xml] = parsed
    return parsed


def _parse_force_xml(xml, name):
    """_parse_force_terms without the cache."""

    root = ET.fromstring(xml)
    is_cmap = root.get("type") == "CMAPTorsionForce"
    for child in root:
        if child.tag not in _TERM_CONTAINERS:
            continue
        terms = list(child)
        atoms = []
        if is_cmap:
            keys = ["a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4"]
            for element in terms:
                atoms.append(
                    tuple(dict.fromkeys(int(element.get(k)) for k in keys)))
        elif terms:
            # every term of a force carries the same attributes
            keys = sorted(
                (k for k in terms[0].attrib if re.fullmatch(r"p\d+", k)),
                key=lambda k: int(k[1:]))
            if not keys:
                return None
            atoms = [
                tuple(int(element.get(k)) for k in keys) for element in terms
            ]
        return {
            "root": root,
            "container": child,
            "terms": terms,
            "atoms": atoms,
            "name": name,
            "type": root.get("type"),
        }
    return None


def _filtered_force(parsed, keep, group):
    """A copy of a parsed force holding only the terms in keep, in the given
    force group."""

    root = parsed["root"]
    new_root = ET.Element(root.tag, dict(root.attrib))
    new_root.text, new_root.tail = root.text, root.tail
    for child in root:
        if child is parsed["container"]:
            container = ET.SubElement(new_root, child.tag, dict(child.attrib))
            container.extend(parsed["terms"][k] for k in sorted(set(keep)))
        else:
            new_root.append(child)
    new_root.set("forceGroup", str(group))
    return mm.XmlSerializer.deserialize(ET.tostring(new_root,
                                                    encoding="unicode"))


def _xml_parameters(elements):
    """(n, 3) array of the q, sig and eps attributes of NonbondedForce
    particle or exception elements."""

    values = [
        (float(e.get("q")), float(e.get("sig")), float(e.get("eps")))
        for e in elements
    ]
    return np.array(values, dtype=np.float64).reshape(-1, 3)


def _empty_copy(system):
    """A system with the same particles and box and no forces."""

    sub = mm.System()
    for i in range(system.getNumParticles()):
        sub.addParticle(system.getParticleMass(i))
    sub.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    return sub


def _same_term(atoms, target):
    """Whether a force term acts on the tracked atoms: the same pair, the
    same angle (either direction) or the same four atoms of a torsion."""

    if len(atoms) == 3:
        return list(atoms) == list(target) or list(atoms) == list(target)[::-1]
    return set(atoms) == set(target)


def _minimum_image(d, box):
    """Minimum image of displacement(s) d (..., 3) in an OpenMM box (rows a,
    b, c in reduced form); unchanged without a box."""

    if box is None:
        return d
    box = np.asarray(box, dtype=np.float64)
    d = d - np.round(d[..., 2:3] / box[2, 2]) * box[2]
    d = d - np.round(d[..., 1:2] / box[1, 1]) * box[1]
    d = d - np.round(d[..., 0:1] / box[0, 0]) * box[0]
    return d


def _angle(x, box):
    """Angle at the middle atom, in degrees."""

    u = _minimum_image(x[0] - x[1], box)
    v = _minimum_image(x[2] - x[1], box)
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _dihedral(x, box):
    """Dihedral angle of four atoms, in degrees (IUPAC sign)."""

    b1 = _minimum_image(x[1] - x[0], box)
    b2 = _minimum_image(x[2] - x[1], box)
    b3 = _minimum_image(x[3] - x[2], box)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    y = np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2))
    return float(np.degrees(np.arctan2(y, np.dot(n1, n2))))

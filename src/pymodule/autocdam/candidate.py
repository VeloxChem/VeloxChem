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
"""``CDAMCandidate``: one metal + dummy-site model, in angstrom and kcal/mol.

The candidate is the source of truth for a CDAM point: geometry, charges, LJ
and C4 terms, masses, and the metal-dummy restraints. ``with_parameters``
returns a copy with fitted values applied; ``to_generated_batch`` turns it into
the Python-side ``GeneratedBatch`` that the generator and calculator use.
Sections, in order:

* geometries: dummy-site unit vectors per coordination geometry. Add one to
  ``GEOMETRIES``; its coordination number is the number of vectors;
* files: PDB and OpenMM ForceField XML text. These exports are optional;
  the XML is also what ``CoordinateGenerator.create_openmm_template`` feeds to
  OpenMM's ``ForceField`` through StringIO, so nothing touches disk;
* the candidate itself.
"""

import math
from copy import deepcopy
from pathlib import Path

from .models import GeneratedBatch, ParameterSet, SystemData, site_parameters_openmm

# ==========================================================================
# geometries
# ==========================================================================


def _unit_directions(vectors):
    """Normalize template vertices so every dummy has the same M--D radius."""
    directions = []
    for vector in vectors:
        length = math.sqrt(sum(value * value for value in vector))
        directions.append(tuple(value / length for value in vector))
    return tuple(directions)


def _square_antiprism_directions():
    """Return vertices of a regular square antiprism on the unit sphere."""
    height = math.sqrt(1.0 / (1.0 + 2.0 * math.sqrt(2.0)))
    radius = math.sqrt(1.0 - height * height)
    diagonal = radius / math.sqrt(2.0)
    return (
        (radius, 0.0, height),
        (0.0, radius, height),
        (-radius, 0.0, height),
        (0.0, -radius, height),
        (diagonal, diagonal, -height),
        (-diagonal, diagonal, -height),
        (-diagonal, -diagonal, -height),
        (diagonal, -diagonal, -height),
    )


def _dodecahedral_directions():
    """Return spherized D2d trigonal-dodecahedron vertices on the unit sphere."""
    q = 0.169022229424176
    r = math.sqrt(q)
    s = math.sqrt((1.0 - q) / (2.0 * q))
    t = math.sqrt(2.0 - 2.0 * q)
    return _unit_directions((
        (t, r, 0.0),
        (-t, r, 0.0),
        (0.0, -r, t),
        (0.0, -r, -t),
        (1.0, -s, 0.0),
        (-1.0, -s, 0.0),
        (0.0, s, 1.0),
        (0.0, s, -1.0),
    ))


GEOMETRIES = {
    "octahedral": (
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, -1.0),
    ),
    "tetrahedral": tuple(
        tuple(value / math.sqrt(3.0)
              for value in direction)
        for direction in ((1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0),
                          (-1.0, -1.0, 1.0))),
    "square_planar": (
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
    ),
    "square_antiprismatic": _square_antiprism_directions(),
    "dodecahedral": _dodecahedral_directions(),
}


def normalize_geometry_name(geometry: str) -> str:
    """``"Square-Planar"`` -> ``"square_planar"``."""
    return geometry.strip().lower().replace("-", "_")


def dummy_angles(directions):
    """Every (i, j, angle_radians) between template directions, i < j."""
    angles = []
    for left in range(len(directions)):
        for right in range(left + 1, len(directions)):
            cosine = sum(
                a * b for a, b in zip(directions[left], directions[right]))
            angles.append((left, right, math.acos(max(-1.0, min(1.0, cosine)))))
    return angles


# ==========================================================================
# candidate files
# ==========================================================================

# Water models with OpenMM's bundled XML: (force-field file, Modeller model name).
WATER_MODELS = {
    "tip3p": ("tip3p.xml", "tip3p"),
    "spce": ("spce.xml", "spce"),
    "spc/e": ("spce.xml", "spce"),
    "tip4pew": ("tip4pew.xml", "tip4pew"),
    "tip4p-ew": ("tip4pew.xml", "tip4pew"),
}


def candidate_pdb(candidate):
    """Return the PDB geometry (angstrom)."""
    names = [candidate.metal.upper(), *candidate.dummy_names()]
    elements = [candidate.metal, *("EP" for _ in candidate.dummy_names())]
    lines = [f"REMARK   1 AutoCDAM candidate {candidate.candidate_id}"]
    for serial, (name, element,
                 position) in enumerate(zip(names, elements,
                                            candidate.coordinates()),
                                        start=1):
        x, y, z = position
        lines.append(
            f"HETATM{serial:5d} {name:<4s} CDM A{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
        )
    lines.extend(f"CONECT{1:5d}{serial:5d}"
                 for serial in range(2, candidate.coordination_number + 2))
    return "\n".join([*lines, "TER", "END", ""])


def candidate_forcefield_xml(c):
    """Return an OpenMM ForceField XML file as plain text (nm, kJ/mol, amu).

    OpenMM's standard ``NonbondedForce`` represents Coulomb plus 12-6 with
    built-in Lorentz-Berthelot mixing, which needs no XML selector. LJ sits on
    the metal only; dummies are charge-only. C4 metadata identifies the
    metal-water pair that the generator adds as a ``CustomNonbondedForce``.
    """
    site = c.openmm_parameters()
    metal, dummies = c.metal_type(), c.dummy_types()
    lines = ["<ForceField>", f"  <!-- AutoCDAM candidate {c.candidate_id} -->"]

    lines.append("  <AtomTypes>")
    lines.append(f'    <Type name="{metal}" class="{metal}" '
                 f'element="{c.metal}" mass="{c.central_mass:.12g}"/>')
    lines.extend(
        f'    <Type name="{t}" class="{t}" mass="{c.dummy_mass:.12g}"/>'
        for t in dummies)
    lines.append("  </AtomTypes>")

    # One CDM residue with metal-dummy bonds.
    lines += [
        "  <Residues>", '    <Residue name="CDM">',
        f'      <Atom name="{c.metal.upper()}" type="{metal}"/>'
    ]
    for name, atom_type in zip(c.dummy_names(), dummies):
        lines.append(f'      <Atom name="{name}" type="{atom_type}"/>')
        lines.append(
            f'      <Bond atomName1="{c.metal.upper()}" atomName2="{name}"/>')
    lines += ["    </Residue>", "  </Residues>"]

    lines.append("  <HarmonicBondForce>")
    lines.extend(f'    <Bond type1="{metal}" type2="{t}" '
                 f'length="{site["distance"]:.12g}" k="{site["bond_k"]:.12g}"/>'
                 for t in dummies)
    lines.append("  </HarmonicBondForce>")

    lines.append("  <HarmonicAngleForce>")
    for left, right, angle in dummy_angles(c.geometries[c.geometry]):
        lines.append(
            f'    <Angle type1="{c.dummy_type(left + 1)}" type2="{metal}" '
            f'type3="{c.dummy_type(right + 1)}" angle="{angle:.12g}" k="{site["angle_k"]:.12g}"/>'
        )
    lines.append("  </HarmonicAngleForce>")

    # OpenMM requires every combined NonbondedForce tag to use the same 1-4
    # scales. The bundled water XML uses these values; the CDAM star has no
    # 1-4 bonded pairs of its own.
    if c.water_model and c.water_model.lower() in WATER_MODELS:
        coulomb14scale, lj14scale = "0.833333", "0.5"
    else:
        coulomb14scale, lj14scale = "0", "0"
    lines += [
        "  <!-- OpenMM NonbondedForce built-in default: Lorentz-Berthelot LJ mixing. -->",
        "  <!-- sigma_ij=(sigma_i+sigma_j)/2; epsilon_ij=sqrt(epsilon_i*epsilon_j). -->",
        "  <!-- No XML mixing-rule selector is required; AutoCDAM metadata only documents the rule. -->",
        f'  <NonbondedForce coulomb14scale="{coulomb14scale}" lj14scale="{lj14scale}">',
        f'    <Atom type="{metal}" charge="{c.central_charge:.12g}" '
        f'sigma="{site["sigma"]:.12g}" epsilon="{site["epsilon"]:.12g}"/>',
    ]
    lines.extend(
        f'    <Atom type="{t}" charge="{c.dummy_charge:.12g}" sigma="0" epsilon="0"/>'
        for t in dummies)
    lines.append("  </NonbondedForce>")

    lines.append(f'  <AutoCDAM potential_form="{c.potential_form}" '
                 f'water_model="{c.water_model or "unspecified"}" '
                 f'c4_metal_water="{site["c4"]:.12g}" '
                 'c4_unit="kJ mol^-1 nm^4" partner_role="water_oxygen" '
                 'lj_mixing_rule="lorentz_berthelot"/>')
    if c.potential_form == "lj12_6_4":
        lines.append(
            "  <!-- Add -C4_metal_water/r^4 for central_metal-water_oxygen pairs. -->"
        )
    return "\n".join([*lines, "</ForceField>", ""])


# ==========================================================================
# candidate
# ==========================================================================

# Physical parameters a proposal may set; see CDAMCandidate.with_parameters for units.
PARAMETER_NAMES = frozenset({
    "dummy_charge",
    "metal_dummy_distance",
    "metal_lj_sigma",
    "metal_lj_epsilon",
    "c4_metal_water",
    "bond_force_constant",
    "angle_force_constant",
})
# Short user-facing names accepted by AutoCDAM(parameters=...).
PARAMETER_ALIASES = {
    "sigma": "metal_lj_sigma",
    "epsilon": "metal_lj_epsilon",
    "c4": "c4_metal_water",
}
# TODO: use the VeloxChem element library for metal masses instead of this table.
METAL_MASSES = {
    "Al": 26.9815385,
    "Ca": 40.078,
    "Co": 58.933194,
    "Cr": 51.9961,
    "Fe": 55.845,
    "Mg": 24.305,
    "Mn": 54.938044,
    "Ni": 58.6934,
    "Zr": 91.224,
    "Zn": 65.38,
}


class CDAMCandidate:
    """One CDAM parameter point: geometry, charges, and potential settings.

    Dummy sites carry charge only. The central charge compensates their total
    charge so the complete ion retains ``formal_charge``. Positive dummy
    charges have no hard-coded upper limit; fitting bounds belong to the study.
    TIP3P, sigma 2 Å and epsilon 0.2 kcal/mol are convenient starting
    defaults, not fitted or transferable literature parameters.
    """

    geometries = GEOMETRIES
    parameter_names = PARAMETER_NAMES
    parameter_aliases = PARAMETER_ALIASES
    metal_masses = METAL_MASSES

    def __init__(
            self,
            metal,
            formal_charge,
            geometry="octahedral",
            coordination_number=None,
            dummy_charge=0.50,
            metal_dummy_distance=0.90,
            potential_form="lj12_6",
            metal_lj_sigma=2.0,
            metal_lj_epsilon=0.2,
            c4_metal_water=0.0,
            water_model="tip3p",
            metal_mass=None,
            dummy_mass=3.0,
            bond_force_constant=800.0,  # TODO: check against VeloxChem defaults
            angle_force_constant=250.0,  # TODO: check against VeloxChem defaults
    ):
        """Store the template settings and validate them.

        Distances and sigma are in angstrom, epsilon in kcal/mol, C4 in
        kcal/mol*angstrom^4, charges in e and masses in amu. The bond force
        constant is in kcal/mol/angstrom^2 and the angle constant in
        kcal/mol/rad^2. ``coordination_number`` defaults to the number of
        directions of ``geometry``; ``metal_mass`` defaults to the
        ``METAL_MASSES`` table and is the mass of the whole ion (dummy masses
        are subtracted from the central site). ``metal_lj_sigma`` and
        ``metal_lj_epsilon`` may be None until forces are built. Raises
        ValueError for an invalid combination (see ``_validate``).
        """
        self.metal = self._normalise_metal(metal)
        self.formal_charge = formal_charge
        self.geometry = normalize_geometry_name(geometry)
        expected = len(self.geometries.get(self.geometry, ()))
        self._coordination_number = expected if coordination_number is None else coordination_number
        self.dummy_charge = dummy_charge
        self.metal_dummy_distance = metal_dummy_distance
        self.potential_form = potential_form
        self.metal_lj_sigma = metal_lj_sigma
        self.metal_lj_epsilon = metal_lj_epsilon
        self.c4_metal_water = c4_metal_water
        self.water_model = water_model
        self.metal_mass = metal_mass if metal_mass is not None else self.metal_masses.get(
            self.metal)
        self.dummy_mass = dummy_mass
        self.bond_force_constant = bond_force_constant
        self.angle_force_constant = angle_force_constant
        self._validate()

    # ------------------------------------------------------------ derived values
    @property
    def coordination_number(self):
        """Number of dummy sites."""
        return self._coordination_number

    @property
    def central_charge(self):
        """Charge on the central metal site in e, so the ion totals ``formal_charge``."""
        return float(
            self.formal_charge) - self.coordination_number * self.dummy_charge

    @property
    def central_mass(self):
        """Mass of the central site in amu: the ion mass minus all dummy masses."""
        return self.metal_mass - self.coordination_number * self.dummy_mass

    @property
    def candidate_id(self):
        """Readable identifier built from metal, charge, geometry, q, r and potential form."""
        return (
            f"{self.metal}{self.formal_charge}_{self.geometry}_"
            f"q{self.dummy_charge:.3f}_r{self.metal_dummy_distance:.3f}_{self.potential_form}"
        )

    def openmm_parameters(self) -> dict[str, float]:
        """Metal LJ, C4, metal-dummy distance and force constants in OpenMM units."""
        return site_parameters_openmm(sigma=self.metal_lj_sigma,
                                      epsilon=self.metal_lj_epsilon,
                                      c4=self.c4_metal_water,
                                      distance=self.metal_dummy_distance,
                                      bond_k=self.bond_force_constant,
                                      angle_k=self.angle_force_constant)

    def coordinates(self):
        """Return metal and dummy coordinates in angstrom."""
        coordinates = [[0.0, 0.0, 0.0]]
        for direction in self.geometries[self.geometry]:
            coordinates.append(
                [self.metal_dummy_distance * value for value in direction])
        return coordinates

    def dummy_names(self):
        """Return the dummy atom names ``D1`` ... ``Dn``."""
        return [f"D{index}" for index in range(1, self.coordination_number + 1)]

    def dummy_types(self):
        """Return the force-field atom types of all dummy sites, in site order."""
        return [
            self.dummy_type(index)
            for index in range(1, self.coordination_number + 1)
        ]

    @staticmethod
    def dummy_type(index):
        """Return the atom type of the 1-based dummy site ``index``."""
        return f"CDAM_D{index}"

    def metal_type(self):
        """Return the force-field atom type of the central metal."""
        return f"CDAM_{self.metal}"

    # ------------------------------------------------------------------- copies
    def with_parameters(self, parameters):
        """Copy this template and apply named physical parameters.

        LJ sigma and M--D distance are in angstrom; epsilon and C4 are in
        kcal/mol and kcal/mol*angstrom^4. Charge is in elementary charges.
        Unknown names are rejected so a scan cannot silently be ineffective.
        """
        values = parameters.values if isinstance(parameters,
                                                 ParameterSet) else parameters
        unknown = set(values) - self.parameter_names
        if unknown:
            raise ValueError(f"unsupported CDAM parameters: {sorted(unknown)}")
        candidate = deepcopy(self)
        for name, value in values.items():
            setattr(candidate, name, float(value))
        candidate._validate()
        candidate.require_lj_parameters()
        return candidate

    def with_geometry(self, geometry, metal_dummy_distance=None):
        """Copy this template with another dummy geometry (and its coordination number)."""
        return type(self)(metal=self.metal,
                          formal_charge=self.formal_charge,
                          geometry=geometry,
                          dummy_charge=self.dummy_charge,
                          metal_dummy_distance=(self.metal_dummy_distance
                                                if metal_dummy_distance is None
                                                else metal_dummy_distance),
                          potential_form=self.potential_form,
                          metal_lj_sigma=self.metal_lj_sigma,
                          metal_lj_epsilon=self.metal_lj_epsilon,
                          c4_metal_water=self.c4_metal_water,
                          water_model=self.water_model,
                          metal_mass=self.metal_mass,
                          dummy_mass=self.dummy_mass,
                          bond_force_constant=self.bond_force_constant,
                          angle_force_constant=self.angle_force_constant)

    # ------------------------------------------------------- batch round trip
    @classmethod
    def from_generated_batch(cls, batch):
        """Recover a CDAM template from authoritative Python-side batch data."""
        data, meta = batch.system, batch.metadata
        nb, bonded = data.nonbonded_parameters, data.bonded_parameters
        return cls(metal=meta["metal"],
                   formal_charge=meta["formal_charge"],
                   geometry=meta["geometry"],
                   dummy_charge=data.charges[1],
                   metal_dummy_distance=bonded["metal_dummy_distance_angstrom"],
                   potential_form=nb["potential_form"],
                   metal_lj_sigma=nb["metal_lj_sigma_angstrom"],
                   metal_lj_epsilon=nb["metal_lj_epsilon_kcal_per_mol"],
                   c4_metal_water=nb["c4_metal_water_kcal_angstrom4_per_mol"],
                   water_model=meta["water_model"],
                   metal_mass=sum(a["mass_amu"] for a in data.atoms),
                   dummy_mass=data.atoms[1]["mass_amu"],
                   bond_force_constant=bonded[
                       "bond_force_constant_kcal_per_mol_angstrom2"],
                   angle_force_constant=bonded[
                       "angle_force_constant_kcal_per_mol_rad2"])

    def to_generated_batch(self):
        """Return this candidate as a one-frame generated batch."""
        metal_atom = {
            "name": self.metal.upper(),
            "element": self.metal,
            "role": "metal",
            "mass_amu": self.central_mass
        }
        dummy_atoms = [{
            "name": name,
            "element": "EP",
            "role": "dummy",
            "mass_amu": self.dummy_mass
        } for name in self.dummy_names()]
        system = SystemData(
            atoms=[metal_atom, *dummy_atoms],
            residues=[{
                "name": "CDM",
                "chain_id": "A",
                "residue_id": 1
            }],
            bonds=[
                (0, index) for index in range(1, self.coordination_number + 1)
            ],
            atom_types=[self.metal_type(), *self.dummy_types()],
            charges=[
                self.central_charge,
                *([self.dummy_charge] * self.coordination_number)
            ],
            bonded_parameters={
                "metal_dummy_distance_angstrom": self.metal_dummy_distance,
                "bond_force_constant_kcal_per_mol_angstrom2":
                    self.bond_force_constant,
                "angle_force_constant_kcal_per_mol_rad2":
                    self.angle_force_constant,
            },
            nonbonded_parameters={
                "potential_form": self.potential_form,
                "metal_lj_sigma_angstrom": self.metal_lj_sigma,
                "metal_lj_epsilon_kcal_per_mol": self.metal_lj_epsilon,
                "c4_metal_water_kcal_angstrom4_per_mol": self.c4_metal_water,
            },
            metadata=self.metadata(),
        )
        return GeneratedBatch(system=system,
                              configurations=[self.coordinates()],
                              metadata=self.metadata())

    def metadata(self):
        """Return a JSON-friendly summary of the candidate (units in the key names)."""
        return {
            "mode": "cdam_candidate",
            "candidate_id": self.candidate_id,
            "metal": self.metal,
            "formal_charge": self.formal_charge,
            "geometry": self.geometry,
            "coordination_number": self.coordination_number,
            "metal_dummy_distance_angstrom": self.metal_dummy_distance,
            "dummy_charge_e": self.dummy_charge,
            "central_charge_e": self.central_charge,
            "potential_form": self.potential_form,
            "water_model": self.water_model,
        }

    # ------------------------------------------------------------------- files
    def write_files(self, directory, stem=None):
        """Write matching PDB and OpenMM XML files."""
        self.require_lj_parameters()
        output_directory = Path(directory)
        output_directory.mkdir(parents=True, exist_ok=True)
        file_stem = stem or self.candidate_id
        pdb_path = output_directory / f"{file_stem}.pdb"
        xml_path = output_directory / f"{file_stem}.xml"
        pdb_path.write_text(self.to_pdb(), encoding="utf-8")
        xml_path.write_text(self.to_openmm_xml(), encoding="utf-8")
        return pdb_path, xml_path

    def to_pdb(self):  # TODO: align with VeloxChem PDB output
        """Return the PDB geometry."""
        return candidate_pdb(self)

    def to_openmm_xml(self):  # TODO: align with VeloxChem XML output
        """Return an OpenMM ForceField XML file as plain text."""
        self.require_lj_parameters()
        return candidate_forcefield_xml(self)

    # -------------------------------------------------------------- validation
    def _validate(self):
        """Raise ValueError if geometry, charges, masses or force terms are inconsistent.

        LJ sigma/epsilon are checked separately by ``require_lj_parameters``
        because a template may leave them unset.
        """
        for name in ("dummy_charge", "metal_dummy_distance", "c4_metal_water",
                     "dummy_mass", "bond_force_constant",
                     "angle_force_constant"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if self.geometry not in self.geometries:
            raise ValueError(f"unsupported geometry {self.geometry!r}")
        if self.coordination_number != len(self.geometries[self.geometry]):
            required = len(self.geometries[self.geometry])
            raise ValueError(
                f"{self.geometry} requires coordination_number={required}")
        if not isinstance(self.formal_charge, int) or self.formal_charge <= 0:
            raise ValueError("formal_charge must be a positive integer")
        if self.potential_form not in ("lj12_6", "lj12_6_4"):
            raise ValueError("potential_form must be 'lj12_6' or 'lj12_6_4'")
        if self.dummy_charge < 0.0:
            raise ValueError("dummy_charge must be nonnegative")
        if self.metal_dummy_distance <= 0.0:
            raise ValueError("metal_dummy_distance must be positive")
        if (self.dummy_mass <= 0 or self.bond_force_constant < 0 or
                self.angle_force_constant < 0 or self.c4_metal_water < 0):
            raise ValueError(
                "masses must be positive and force coefficients nonnegative")
        if self.metal_mass is None:
            raise ValueError(
                f"provide metal_mass for unsupported metal {self.metal}")
        if self.metal_mass <= self.coordination_number * self.dummy_mass:
            raise ValueError("metal mass must exceed total dummy mass")
        if self.potential_form == "lj12_6" and self.c4_metal_water != 0.0:
            raise ValueError("lj12_6 requires c4_metal_water=0.0")

    def require_lj_parameters(self):
        """Refuse to build forces without a positive sigma and nonnegative epsilon."""
        if self.metal_lj_sigma is None or self.metal_lj_epsilon is None:
            raise ValueError(
                "writing OpenMM XML requires metal_lj_sigma and metal_lj_epsilon"
            )
        if (not math.isfinite(self.metal_lj_sigma) or
                not math.isfinite(self.metal_lj_epsilon) or
                self.metal_lj_sigma <= 0.0 or self.metal_lj_epsilon < 0.0):
            raise ValueError(
                "metal_lj_sigma must be positive and metal_lj_epsilon cannot be negative"
            )

    @staticmethod
    def _normalise_metal(metal):
        """Return the element symbol capitalised (``zn`` -> ``Zn``); raise ValueError if not alphabetic."""
        metal = metal.strip()
        if not metal or not metal.isalpha():
            raise ValueError("metal must be an element symbol, such as 'Zn'")
        return metal[0].upper() + metal[1:].lower()

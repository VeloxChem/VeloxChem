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
"""Turn CDAM parameter proposals into coordinate batches and OpenMM systems.

``CoordinateGenerator`` builds vacuum and solvated OpenMM systems once per
candidate template, then clones and updates them for each proposal. Units are
converted from the candidate's angstrom/kcal/mol to OpenMM's nm/kJ/mol here.
PDB/XML export is optional; fitting uses Python objects directly.

The candidate model itself lives in ``candidate.py``.
"""

import io
import math
import sys

import numpy as np
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass

from ..errorhandler import assert_msg_critical
from .geometry import WATER_RESIDUES
from .candidate import WATER_MODELS
from .models import (OpenMMSystems, ParameterSet, default_ostream,
                     ANGSTROM_TO_NM, site_parameters_openmm)

__all__ = ["CoordinateGenerator"]

C4_FORCE_NAME = "CDAM_1264"
TEMPLATE_CUTOFF_NM = 1.0
MODEL_METADATA_KEYS = ("metal", "formal_charge", "geometry", "potential_form",
                       "water_model")


class CoordinateGenerator:
    """Create candidate coordinate batches and update their OpenMM systems.

    Batch generation updates candidate coordinates and cloned OpenMM systems
    from a common template. File export remains an optional candidate utility.
    """

    def __init__(self, ostream=None):
        """Use ``ostream`` for progress messages (the package default when None)."""
        assert_msg_critical("openmm" in sys.modules,
                            f"{type(self).__name__}: OpenMM is required")
        self.ostream = default_ostream(ostream)

    def generate_parameter_batch(self,
                                 candidate,
                                 proposals,
                                 *,
                                 openmm_template=None,
                                 sample_ids=None):
        """Generate one independent coordinate/system batch per parameter set.

        This does not write PDB/XML. A fixed solvent box can be reused: its
        coordinates are copied, solute coordinates and forces are updated, and
        the calculator minimizes each candidate before collecting properties.
        """
        proposals = list(proposals)
        ids = (list(sample_ids) if sample_ids is not None else
               [f"sample-{i:06d}" for i in range(len(proposals))])
        if len(ids) != len(proposals):
            raise ValueError("sample_ids must match the number of proposals")
        batches = []
        for parameters, sample_id in zip(proposals, ids):
            if not isinstance(parameters, ParameterSet):
                parameters = ParameterSet(parameters)
            batch = candidate.with_parameters(parameters).to_generated_batch()
            batch.parameters = parameters
            batch.metadata["sample_id"] = sample_id
            if openmm_template is not None:
                batch.openmm_systems = self.update_openmm_systems(
                    openmm_template, batch)
            batches.append(batch)
        return batches

    def create_openmm_template(self, candidate, *, padding_nm=1.2):
        """Build vacuum and solvated OpenMM objects once, entirely in memory.

        OpenMM's ForceField reads the candidate XML from StringIO; no candidate
        files are generated. Supported water models use OpenMM's bundled XML.
        More complex solvent/framework systems can be supplied as OpenMMSystems.
        """
        candidate.require_lj_parameters()
        if candidate.water_model is None or candidate.water_model.lower(
        ) not in WATER_MODELS:
            raise ValueError(
                "provide a supported water_model or a prepared OpenMMSystems template"
            )
        if not math.isfinite(padding_nm) or padding_nm <= TEMPLATE_CUTOFF_NM:
            raise ValueError(
                "padding_nm must exceed the 1.0 nm nonbonded cutoff")
        topology = app.Topology()
        residue = topology.addResidue("CDM", topology.addChain("A"))
        metal = topology.addAtom(candidate.metal.upper(),
                                 app.element.get_by_symbol(candidate.metal),
                                 residue)
        for name in candidate.dummy_names():
            topology.addBond(metal, topology.addAtom(name, None, residue))
        positions = unit.Quantity(
            np.asarray(candidate.coordinates()) * ANGSTROM_TO_NM,
            unit.nanometer)
        xml = candidate.to_openmm_xml()
        vacuum = app.ForceField(io.StringIO(xml)).createSystem(
            topology, nonbondedMethod=app.NoCutoff)
        water_xml, model = WATER_MODELS[candidate.water_model.lower()]
        forcefield = app.ForceField(io.StringIO(xml), water_xml)
        modeller = app.Modeller(topology, positions)
        modeller.addSolvent(forcefield,
                            model=model,
                            padding=padding_nm * unit.nanometer,
                            neutralize=False)
        solvated = forcefield.createSystem(modeller.topology,
                                           nonbondedMethod=app.PME,
                                           nonbondedCutoff=TEMPLATE_CUTOFF_NM *
                                           unit.nanometer,
                                           constraints=app.HBonds)
        template = OpenMMSystems(vacuum, topology, positions, solvated,
                                 modeller.topology, modeller.positions)
        return self.update_openmm_systems(template,
                                          candidate.to_generated_batch())

    @staticmethod
    def update_openmm_systems(template, batch):
        """Clone and update both systems from a generated candidate's data.

        Updates occur before creating a Context. Changing constraints/geometry
        requires a fresh Context, rather than updateParametersInContext alone.
        The caller's template and every other proposal remain untouched.
        """
        if batch.metadata.get("mode") != "cdam_candidate":
            raise ValueError(
                "OpenMM updates require a generated CDAM candidate")
        model_metadata = {
            key: batch.metadata[key] for key in MODEL_METADATA_KEYS
        }
        if template.model_metadata and template.model_metadata != model_metadata:
            raise ValueError(
                "template ion/geometry/potential/water model differs from the generated candidate"
            )
        forces = _site_forces_openmm_units(batch.system)
        output = []
        for original, topology, positions in (
            (template.vacuum_system, template.vacuum_topology,
             template.vacuum_positions),
            (template.solvated_system, template.solvated_topology,
             template.solvated_positions),
        ):
            output.extend(
                _updated_system(original,
                                topology,
                                positions,
                                batch,
                                forces,
                                require_water=original
                                is template.solvated_system))
        return OpenMMSystems(*output, model_metadata=model_metadata)


# ------------------------------------------------------------ system updates


def _site_forces_openmm_units(data):
    """The candidate's force parameters converted to nm and kJ/mol."""
    nb, bonded = data.nonbonded_parameters, data.bonded_parameters
    return {
        "potential_form": nb["potential_form"],
        **site_parameters_openmm(sigma=nb["metal_lj_sigma_angstrom"],
                                 epsilon=nb["metal_lj_epsilon_kcal_per_mol"],
                                 c4=nb["c4_metal_water_kcal_angstrom4_per_mol"],
                                 distance=bonded["metal_dummy_distance_angstrom"],
                                 bond_k=bonded["bond_force_constant_kcal_per_mol_angstrom2"],
                                 angle_k=bonded["angle_force_constant_kcal_per_mol_rad2"])
    }


def _updated_system(original, topology, original_positions, batch, forces, *,
                    require_water):
    """One cloned (system, topology, positions) with the candidate's site applied."""
    data = batch.system
    system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(original))
    atoms = list(topology.atoms())
    indices = _site_indices(atoms, data)
    if system.getNumParticles() != len(atoms):
        raise ValueError("OpenMM system and topology particle counts differ")
    metal, dummies = indices[0], indices[1:]
    positions = _place_site(original_positions, indices,
                            batch.configurations[0])
    for index, atom in zip(indices, data.atoms):
        system.setParticleMass(index, float(atom["mass_amu"]))
    nonbonded = _single_nonbonded_force(system)
    _set_site_nonbonded(nonbonded, indices, data.charges, forces)
    _set_site_restraints(system, metal, dummies, forces)
    _set_c4(system,
            nonbonded,
            atoms,
            metal,
            forces,
            require_water=require_water)
    return system, topology, unit.Quantity(positions, unit.nanometer)


def _site_indices(atoms, data):
    """Topology indices of the CDAM atoms, ordered like the candidate (metal first)."""
    solute = [a for a in atoms if a.residue.name == "CDM"]
    names = [a["name"] for a in data.atoms]
    if len(solute) != len(names) or {a.name for a in solute} != set(names):
        raise ValueError(
            "template CDM atom names/count differ from the generated candidate")
    return [next(a.index for a in solute if a.name == name) for name in names]


def _place_site(original_positions, indices, site_coordinates_angstrom):
    """Copy the template positions and move the site atoms onto the candidate geometry."""
    positions = np.array(original_positions.value_in_unit(unit.nanometer),
                         dtype=float,
                         copy=True)
    center = positions[indices[0]].copy()
    for index, coordinates in zip(indices, site_coordinates_angstrom):
        positions[index] = center + np.array(coordinates) * ANGSTROM_TO_NM
    return positions


def _single_nonbonded_force(system):
    """Return the system's only ``NonbondedForce``; raise ValueError if there is not exactly one."""
    nonbonded = [
        f for f in system.getForces() if isinstance(f, mm.NonbondedForce)
    ]
    if len(nonbonded) != 1:
        raise ValueError(
            "template requires exactly one standard NonbondedForce")
    return nonbonded[0]


def _set_site_nonbonded(force, indices, charges, forces):
    """Metal charge + LJ, charge-only dummies, and no interactions inside the site."""
    metal, dummies = indices[0], indices[1:]
    force.setParticleParameters(metal, charges[0], forces["sigma"],
                                forces["epsilon"])
    for index, charge in zip(dummies, charges[1:]):
        force.setParticleParameters(index, charge, 0.0, 0.0)
    # CDAM internal Coulomb/LJ interactions are excluded (star topology).
    # Refuse exceptions to the environment rather than retaining stale charges.
    for i in range(force.getNumExceptions()):
        left, right, _, sigma, _ = force.getExceptionParameters(i)
        if left in indices and right in indices:
            force.setExceptionParameters(i, left, right, 0.0, sigma, 0.0)
        elif left in indices or right in indices:
            raise ValueError(
                "template contains unsupported CDAM-environment bonded exceptions"
            )


def _set_site_restraints(system, metal, dummies, forces):
    """Metal-dummy constraints at the candidate distance; bond and D-M-D angle constants."""
    existing = {}
    for i in range(system.getNumConstraints()):
        left, right, _ = system.getConstraintParameters(i)
        existing[frozenset((left, right))] = i
    for dummy in dummies:
        index = existing.get(frozenset((metal, dummy)))
        if index is None:
            system.addConstraint(metal, dummy, forces["distance"])
        else:
            system.setConstraintParameters(index, metal, dummy,
                                           forces["distance"])
    for bonded in system.getForces():
        if isinstance(bonded, mm.HarmonicBondForce):
            for i in range(bonded.getNumBonds()):
                left, right, _, _ = bonded.getBondParameters(i)
                if metal in (left, right) and (left in dummies or
                                               right in dummies):
                    bonded.setBondParameters(i, left, right, forces["distance"],
                                             forces["bond_k"])
        elif isinstance(bonded, mm.HarmonicAngleForce):
            for i in range(bonded.getNumAngles()):
                left, center, right, angle, _ = bonded.getAngleParameters(i)
                if center == metal and left in dummies and right in dummies:
                    bonded.setAngleParameters(i, left, center, right, angle,
                                              forces["angle_k"])


def _set_c4(system, nonbonded, atoms, metal, forces, *, require_water):
    """Add or update the -C4/r^4 metal-water-oxygen term (12-6-4 only)."""
    c4_forces = [
        f for f in system.getForces() if
        isinstance(f, mm.CustomNonbondedForce) and f.getName() == C4_FORCE_NAME
    ]
    if len(c4_forces) > 1:
        raise ValueError("duplicate CDAM C4 forces in template")
    if forces["potential_form"] != "lj12_6_4":
        if c4_forces:
            raise ValueError("12-6 template unexpectedly contains a C4 force")
        return
    if c4_forces:
        c4_force = c4_forces[0]
        if (c4_force.getNumPerParticleParameters() != 1 or
                c4_force.getPerParticleParameterName(0) != "C4"):
            raise ValueError("incompatible CDAM C4 force in template")
    else:
        c4_force = new_c4_force(system, nonbonded, len(atoms))
    oxygens = [
        a.index for a in atoms if a.residue.name in WATER_RESIDUES and
        a.element is not None and a.element.symbol == "O"
    ]
    if not oxygens and require_water:
        raise ValueError("C4 template has no recognized water oxygens")
    set_c4_pairs(c4_force, len(atoms), metal, oxygens, forces["c4"])
    for i in range(c4_force.getNumGlobalParameters()):
        if c4_force.getGlobalParameterName(i) == "c4_lambda":
            c4_force.setGlobalParameterDefaultValue(i, 0.0)


def set_c4_pairs(force, n_particles, metal, partners, c4) -> None:
    """Give the metal ``c4`` and each partner oxygen 1, and restrict pairs to metal-partner.

    Without the interaction group, partner-partner pairs (e.g. water O-O) would
    also feel -1*1/r^4 and the solvent would collapse.
    """
    partners = sorted(partners)
    members = set(partners)
    for index in range(n_particles):
        force.setParticleParameters(
            index, [c4 if index == metal else float(index in members)])
    if force.getNumInteractionGroups() == 0:
        force.addInteractionGroup([metal], partners)
    else:
        force.setInteractionGroupParameters(0, [metal], partners)


def new_c4_force(system, nonbonded, n_particles):
    """Add an empty C4 force (``c4_lambda`` for FEP decoupling, the nonbonded exclusions).

    Shared by the aqueous template and framework systems; callers set the
    per-particle C4 values and the metal-partner interaction group.
    """
    force = mm.CustomNonbondedForce("-(1-c4_lambda)*C41*C42/r^4")
    force.setName(C4_FORCE_NAME)
    force.addGlobalParameter("c4_lambda", 0.0)
    force.addPerParticleParameter("C4")
    for _ in range(n_particles):
        force.addParticle([0.0])
    for i in range(nonbonded.getNumExceptions()):
        left, right, *_ = nonbonded.getExceptionParameters(i)
        force.addExclusion(left, right)
    if nonbonded.getNonbondedMethod() == mm.NonbondedForce.NoCutoff:
        force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    else:
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(nonbonded.getCutoffDistance())
    system.addForce(force)
    return force

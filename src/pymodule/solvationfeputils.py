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

from collections import Counter
from io import StringIO
from contextlib import redirect_stderr
from xml.etree import ElementTree
import numpy as np

with redirect_stderr(StringIO()) as fg_err:
    try:
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError:
        pass

# Force groups used for the per-term energy comparison of the solute
_TERM_GROUPS = {'bonds': 1, 'angles': 2, 'torsions': 3, 'custom bonded terms': 4, 'nonbonded terms': 5}

_CUSTOM_BONDED = ('CustomBondForce', 'CustomAngleForce', 'CustomTorsionForce', 'CustomCompoundBondForce')

# Standard bonded forces and their term class in the energy comparison
_STANDARD_BONDED = {'HarmonicBondForce': 'bonds', 'HarmonicAngleForce': 'angles', 'PeriodicTorsionForce': 'torsions',
                    'RBTorsionForce': 'torsions', 'CMAPTorsionForce': 'torsions'}


def _value(quantity, target_unit=None):
    """
    Returns a plain float from an OpenMM quantity (MD units if no unit is given).
    """

    if not hasattr(quantity, 'value_in_unit'):
        return float(quantity)
    if target_unit is None:
        return quantity.value_in_unit_system(unit.md_unit_system)
    return quantity.value_in_unit(target_unit)


def _get_nonbonded_force(system):
    """
    Returns the only NonbondedForce of a system.
    """

    forces = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError(f'Expected one NonbondedForce, found {len(forces)}.')
    return forces[0]


def _custom_bonded_terms(force):
    """
    Returns (term getter, number of terms, per-term parameter names, adder name) of a custom
    bonded force. The getter returns (particles, parameters) for a term.
    """

    name = type(force).__name__
    if name == 'CustomBondForce':
        get = lambda i: (force.getBondParameters(i)[:2], force.getBondParameters(i)[2])
        params = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
        return get, force.getNumBonds(), params, 'addBond'
    if name == 'CustomAngleForce':
        get = lambda i: (force.getAngleParameters(i)[:3], force.getAngleParameters(i)[3])
        params = [force.getPerAngleParameterName(k) for k in range(force.getNumPerAngleParameters())]
        return get, force.getNumAngles(), params, 'addAngle'
    if name == 'CustomTorsionForce':
        get = lambda i: (force.getTorsionParameters(i)[:4], force.getTorsionParameters(i)[4])
        params = [force.getPerTorsionParameterName(k) for k in range(force.getNumPerTorsionParameters())]
        return get, force.getNumTorsions(), params, 'addTorsion'
    get = lambda i: tuple(force.getBondParameters(i))
    params = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
    return get, force.getNumBonds(), params, 'addBond'


def _restrict_custom_bonded_force(force, local):
    """
    Returns a copy of a custom bonded force with only the terms whose particles all belong to
    the solute, renumbered to solute indices. The copy is made through the XML serialization, so
    tabulated functions, global and per-term parameters are kept.
    """

    root = ElementTree.fromstring(mm.XmlSerializer.serialize(force))
    for container in ('Bonds', 'Angles', 'Torsions'):
        node = root.find(container)
        if node is None:
            continue
        for term in list(node):
            keys = [k for k in term.attrib if k[0] == 'p' and k[1:].isdigit()]
            particles = [int(term.attrib[k]) for k in keys]
            if all(p in local for p in particles):
                for k, p in zip(keys, particles):
                    term.set(k, str(local[p]))
            else:
                node.remove(term)
    root.set('usesPeriodic', '0')
    return mm.XmlSerializer.deserialize(ElementTree.tostring(root, encoding='unicode'))


def _is_restrictable(force):
    """
    Returns True if _restrict_force can copy a force into the solute-only system.
    """

    if type(force).__name__ in _STANDARD_BONDED or type(force).__name__ in _CUSTOM_BONDED:
        return True
    if isinstance(force, mm.RMSDForce):
        return True
    if isinstance(force, mm.CustomCVForce):
        return all(_is_restrictable(force.getCollectiveVariable(k)) for k in range(force.getNumCollectiveVariables()))
    return False


def _restrict_force(force, local):
    """
    Returns a copy of a bonded-type force with only the terms whose particles all belong to the
    solute, renumbered to solute indices (local: solvated index -> solute index, in solute order).
    RMSDForce and CustomCVForce must act on solute atoms only; a CustomCVForce is copied through
    its XML serialization with every collective variable restricted, so its expression, parameters
    and tabulated functions are kept.
    """

    name = type(force).__name__

    if name in _CUSTOM_BONDED:
        return _restrict_custom_bonded_force(force, local)

    if isinstance(force, mm.HarmonicBondForce):
        copy_force = mm.HarmonicBondForce()
        for i in range(force.getNumBonds()):
            a, b, r0, k = force.getBondParameters(i)
            if a in local and b in local:
                copy_force.addBond(local[a], local[b], r0, k)
        return copy_force

    if isinstance(force, mm.HarmonicAngleForce):
        copy_force = mm.HarmonicAngleForce()
        for i in range(force.getNumAngles()):
            a, b, c, theta0, k = force.getAngleParameters(i)
            if a in local and b in local and c in local:
                copy_force.addAngle(local[a], local[b], local[c], theta0, k)
        return copy_force

    if isinstance(force, mm.PeriodicTorsionForce):
        copy_force = mm.PeriodicTorsionForce()
        for i in range(force.getNumTorsions()):
            a, b, c, d, n, phase, k = force.getTorsionParameters(i)
            if all(x in local for x in (a, b, c, d)):
                copy_force.addTorsion(local[a], local[b], local[c], local[d], n, phase, k)
        return copy_force

    if isinstance(force, mm.RBTorsionForce):
        copy_force = mm.RBTorsionForce()
        for i in range(force.getNumTorsions()):
            a, b, c, d, *coefficients = force.getTorsionParameters(i)
            if all(x in local for x in (a, b, c, d)):
                copy_force.addTorsion(local[a], local[b], local[c], local[d], *coefficients)
        return copy_force

    if isinstance(force, mm.CMAPTorsionForce):
        copy_force = mm.CMAPTorsionForce()
        for i in range(force.getNumMaps()):
            copy_force.addMap(*force.getMapParameters(i))
        for i in range(force.getNumTorsions()):
            m, *atoms = force.getTorsionParameters(i)
            if all(x in local for x in atoms):
                copy_force.addTorsion(m, *[local[x] for x in atoms])
        return copy_force

    if isinstance(force, mm.RMSDForce):
        reference = force.getReferencePositions().value_in_unit(unit.nanometer)
        particles = list(force.getParticles()) or list(range(len(reference)))
        if not all(p in local for p in particles):
            raise ValueError(f'RMSDForce "{force.getName()}" acts on atoms outside the solute.')
        return mm.RMSDForce([reference[a] for a in local], [local[p] for p in particles])

    if isinstance(force, mm.CustomCVForce):
        root = ElementTree.fromstring(mm.XmlSerializer.serialize(force))
        variables = root.find('CollectiveVariables')
        for k, node in enumerate(variables):
            restricted = _restrict_force(force.getCollectiveVariable(k), local)
            node.remove(node.find('Force'))
            node.append(ElementTree.fromstring(mm.XmlSerializer.serialize(restricted)))
        return mm.XmlSerializer.deserialize(ElementTree.tostring(root, encoding='unicode'))

    raise ValueError(f'Force "{force.getName()}" ({name}) cannot be copied into the solute-only system.')


def _restrict_to_solute(system, solute_indices):
    """
    Returns the solute's own terms of a system as a non-periodic system of the solute atoms:
    masses, constraints between solute atoms, bonded terms (standard and custom) whose atoms all
    belong to the solute, RMSDForce and CustomCVForce terms that act on solute atoms only, and the
    solute's NonbondedForce particles and exceptions (NoCutoff). One force group per term class
    (_TERM_GROUPS).

    Raises ValueError for virtual sites in the solute and for NonbondedForce parameter offsets.
    """

    local = {a: k for k, a in enumerate(solute_indices)}
    restricted = mm.System()

    for a in solute_indices:
        if system.isVirtualSite(a):
            raise ValueError(f'Solute particle {a} is a virtual site; this is not supported.')
        restricted.addParticle(system.getParticleMass(a))
    for i in range(system.getNumConstraints()):
        a, b, distance = system.getConstraintParameters(i)
        if a in local and b in local:
            restricted.addConstraint(local[a], local[b], distance)

    new_forces = []
    for force in system.getForces():
        name = type(force).__name__
        if isinstance(force, mm.NonbondedForce):
            if force.getNumParticleParameterOffsets() > 0 or force.getNumExceptionParameterOffsets() > 0:
                raise ValueError('NonbondedForce parameter offsets are not supported.')
            nonbonded = mm.NonbondedForce()
            nonbonded.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            for a in solute_indices:
                nonbonded.addParticle(*force.getParticleParameters(a))
            for i in range(force.getNumExceptions()):
                a, b, charge_prod, sigma, epsilon = force.getExceptionParameters(i)
                if a in local and b in local:
                    nonbonded.addException(local[a], local[b], charge_prod, sigma, epsilon)
            new_forces.append((nonbonded, 'nonbonded terms'))
        elif name in _STANDARD_BONDED:
            new_forces.append((_restrict_force(force, local), _STANDARD_BONDED[name]))
        elif name in _CUSTOM_BONDED:
            new_forces.append((_restrict_force(force, local), 'custom bonded terms'))
        elif isinstance(force, (mm.RMSDForce, mm.CustomCVForce)):
            particles = _force_particles(force, system.getNumParticles())
            if particles and particles <= set(local) and _is_restrictable(force):
                new_forces.append((_restrict_force(force, local), 'custom bonded terms'))

    for force, group in new_forces:
        force.setForceGroup(_TERM_GROUPS[group])
        restricted.addForce(force)

    return restricted


def extract_solute_system(system, topology, positions, solute_indices):
    """
    Cuts the solute out of a pre-built solvated system, as the system of the vacuum leg.

    The solute keeps exactly the terms it has in the solvated system: masses, constraints, bonded
    terms (standard and custom) within the solute, charges, Lennard-Jones parameters and
    exceptions; nonbonded interactions without cutoff. Terms between the solute and its environment
    (restraints, applied fields, solvent) have no counterpart in vacuum and are not included.

    :param system:
        The solvated OpenMM system.
    :param topology:
        The topology of the solvated system.
    :param positions:
        The positions of the solvated system.
    :param solute_indices:
        The solute's particle indices.

    :return:
        A tuple (system, topology, positions) of the solute alone.
    """

    vacuum_system = _restrict_to_solute(system, list(solute_indices))
    for force in vacuum_system.getForces():
        force.setForceGroup(0)
    vacuum_system.addForce(mm.CMMotionRemover())

    solute = set(solute_indices)
    modeller = app.Modeller(topology, positions)
    modeller.delete([atom for atom in topology.atoms() if atom.index not in solute])
    vacuum_topology = modeller.getTopology()
    vacuum_topology.setPeriodicBoxVectors(None)

    return vacuum_system, vacuum_topology, modeller.getPositions()


def _parameter_tables(system):
    """
    Returns the parameters of a solute-only system in a canonical, comparable form.
    """

    def r(x):
        return round(float(x), 6)

    tables = {
        'masses': [_value(system.getParticleMass(i), unit.dalton) for i in range(system.getNumParticles())],
        'constraints': Counter(), 'bonds': Counter(), 'angles': Counter(), 'torsions': Counter(),
        'charges': [], 'lj': [], 'exceptions': {},
    }
    for i in range(system.getNumConstraints()):
        a, b, distance = system.getConstraintParameters(i)
        tables['constraints'][(min(a, b), max(a, b), r(_value(distance, unit.nanometer)))] += 1
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            for i in range(force.getNumBonds()):
                a, b, r0, k = force.getBondParameters(i)
                tables['bonds'][(min(a, b), max(a, b), r(_value(r0, unit.nanometer)), r(_value(k)))] += 1
        elif isinstance(force, mm.HarmonicAngleForce):
            for i in range(force.getNumAngles()):
                a, b, c, theta0, k = force.getAngleParameters(i)
                key = (a, b, c) if a < c else (c, b, a)
                tables['angles'][key + (r(_value(theta0, unit.radian)), r(_value(k)))] += 1
        elif isinstance(force, mm.PeriodicTorsionForce):
            for i in range(force.getNumTorsions()):
                a, b, c, d, n, phase, k = force.getTorsionParameters(i)
                key = min((a, b, c, d), (d, c, b, a))
                tables['torsions'][key + (n, r(_value(phase, unit.radian)), r(_value(k)))] += 1
        elif isinstance(force, mm.NonbondedForce) and force.getNumParticles() == system.getNumParticles():
            for i in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                tables['charges'].append(_value(charge, unit.elementary_charge))
                tables['lj'].append((_value(sigma, unit.nanometer), _value(epsilon, unit.kilojoule_per_mole)))
            for i in range(force.getNumExceptions()):
                a, b, charge_prod, sigma, epsilon = force.getExceptionParameters(i)
                tables['exceptions'][(min(a, b), max(a, b))] = (_value(charge_prod), _value(sigma, unit.nanometer),
                                                               _value(epsilon, unit.kilojoule_per_mole))
    return tables


def _term_energies(system, coordinates):
    """
    Returns the energy of every term class of a solute-only system (Reference platform).
    """

    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    context.setPositions(coordinates * unit.nanometer)
    return {name: context.getState(getEnergy=True, groups={group}).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole) for name, group in _TERM_GROUPS.items()}


def compare_solute(system, solute_system, solute_indices, solute_coordinates, num_conformations=4,
                   displacement=0.01, tolerance=1e-4, seed=1):
    """
    Compares the solute of a pre-built solvated system with the solute-only (vacuum) system.

    The free energy cycle closes only if the solute Hamiltonian is the same in both legs: charges,
    Lennard-Jones parameters, excluded and 1-4 pairs, bonded terms and constraints. Readable messages
    come from a comparison of the parameters; the verdict for the energy terms comes from their
    energies at the given solute geometry and at randomly displaced geometries, so that equivalent
    terms written differently (e.g. reversed torsions) are accepted.

    :param system:
        The solvated OpenMM system (e.g. from system.xml).
    :param solute_system:
        The OpenMM system of the solute alone (vacuum leg).
    :param solute_indices:
        The solute's particle indices in the solvated system, in the order of the solute system.
    :param solute_coordinates:
        The solute coordinates in nm, shape (number of solute atoms, 3).
    :param num_conformations:
        The number of geometries for the energy comparison (the first is the given one).
    :param displacement:
        The standard deviation of the random displacements in nm.
    :param tolerance:
        The largest accepted energy difference in kJ/mol.
    :param seed:
        The seed of the random displacements.

    :return:
        A tuple (errors, notes) of lists of strings. Errors make the two legs' solute Hamiltonians
        differ. Notes are differences without effect on the free energy (masses change only the
        dynamics; the sigma of atoms without Lennard-Jones interaction only enters a soft-core term).
    """

    errors, notes = [], []
    num_atoms = len(solute_indices)
    if solute_system.getNumParticles() != num_atoms:
        errors.append(f'The solute system has {solute_system.getNumParticles()} particles, '
                      f'the solute in the solvated system has {num_atoms}.')
        return errors, notes

    solvated = _restrict_to_solute(system, list(solute_indices))
    vacuum = _restrict_to_solute(solute_system, list(range(num_atoms)))
    tables_s, tables_v = _parameter_tables(solvated), _parameter_tables(vacuum)

    mass_diff = np.abs(np.array(tables_s['masses']) - np.array(tables_v['masses']))
    if mass_diff.max() > 1e-6:
        notes.append(f'Masses differ for {int((mass_diff > 1e-6).sum())} solute atom(s) (largest difference '
                     f'{mass_diff.max():.4f} Da); this only affects the dynamics.')

    charge_diff = np.abs(np.array(tables_s['charges']) - np.array(tables_v['charges']))
    if charge_diff.max() > 1e-6:
        errors.append(f'Charges differ for {int((charge_diff > 1e-6).sum())} solute atom(s) '
                      f'(largest difference {charge_diff.max():.2e} e).')

    for k, ((sigma_s, eps_s), (sigma_v, eps_v)) in enumerate(zip(tables_s['lj'], tables_v['lj'])):
        if abs(eps_s - eps_v) > 1e-8 * max(1.0, abs(eps_v)) or (eps_v != 0.0 and abs(sigma_s - sigma_v) > 1e-6):
            errors.append(f'Lennard-Jones parameters of solute atom {k} differ: sigma {sigma_s:.6f} and '
                          f'{sigma_v:.6f} nm, epsilon {eps_s:.6f} and {eps_v:.6f} kJ/mol.')
        elif eps_v == 0.0 and abs(sigma_s - sigma_v) > 1e-6:
            notes.append(f'Solute atom {k} has epsilon = 0 and different placeholder sigma values '
                         f'({sigma_s:g} and {sigma_v:g} nm).')

    pairs_s, pairs_v = set(tables_s['exceptions']), set(tables_v['exceptions'])
    if pairs_s != pairs_v:
        errors.append(f'Excluded and 1-4 pairs differ: {len(pairs_s - pairs_v)} pairs only in the solvated '
                      f'system, {len(pairs_v - pairs_s)} only in the solute system.')
    else:
        differing = [p for p in pairs_s
                     if not np.allclose(tables_s['exceptions'][p], tables_v['exceptions'][p], rtol=1e-6, atol=1e-8)
                     and not (tables_s['exceptions'][p][0] == tables_v['exceptions'][p][0] == 0.0 and
                              tables_s['exceptions'][p][2] == tables_v['exceptions'][p][2] == 0.0)]
        if differing:
            errors.append(f'Parameters of {len(differing)} excluded or 1-4 pairs differ (e.g. the 1-4 scaling '
                          'factors).')

    if tables_s['constraints'] != tables_v['constraints']:
        errors.append(f'Constraints differ ({sum(tables_s["constraints"].values())} and '
                      f'{sum(tables_v["constraints"].values())} constrained solute bonds).')
    for name in ('bonds', 'angles', 'torsions'):
        if tables_s[name] != tables_v[name]:
            notes.append(f'The {name} are written differently ({sum(tables_s[name].values())} and '
                         f'{sum(tables_v[name].values())} terms); see the energy comparison.')

    rng = np.random.default_rng(seed)
    coordinates = np.asarray(solute_coordinates, dtype=float)
    largest = {name: 0.0 for name in _TERM_GROUPS}
    for c in range(num_conformations):
        xyz = coordinates if c == 0 else coordinates + rng.normal(scale=displacement, size=coordinates.shape)
        energies_s, energies_v = _term_energies(solvated, xyz), _term_energies(vacuum, xyz)
        for name in _TERM_GROUPS:
            largest[name] = max(largest[name], abs(energies_s[name] - energies_v[name]))
    for name, difference in largest.items():
        if difference > tolerance:
            errors.append(f'The solute {name} differ in energy by up to {difference:.4g} kJ/mol.')

    return errors, notes


def _force_particles(force, num_particles):
    """
    Returns the particles a force acts on, or None if it can act on every particle.
    """

    if isinstance(force, mm.CustomExternalForce):
        return {force.getParticleParameters(i)[0] for i in range(force.getNumParticles())}
    if isinstance(force, mm.CustomCentroidBondForce):
        particles = set()
        for g in range(force.getNumGroups()):
            particles.update(force.getGroupParameters(g)[0])
        return particles
    if type(force).__name__ in _CUSTOM_BONDED:
        get, num_terms, _, _ = _custom_bonded_terms(force)
        return {x for i in range(num_terms) for x in get(i)[0]}
    if isinstance(force, mm.HarmonicBondForce):
        return {x for i in range(force.getNumBonds()) for x in force.getBondParameters(i)[:2]}
    if isinstance(force, mm.HarmonicAngleForce):
        return {x for i in range(force.getNumAngles()) for x in force.getAngleParameters(i)[:3]}
    if isinstance(force, (mm.PeriodicTorsionForce, mm.RBTorsionForce)):
        return {x for i in range(force.getNumTorsions()) for x in force.getTorsionParameters(i)[:4]}
    if isinstance(force, mm.CMAPTorsionForce):
        return {x for i in range(force.getNumTorsions()) for x in force.getTorsionParameters(i)[1:]}
    if isinstance(force, mm.RMSDForce):
        particles = list(force.getParticles())
        return set(particles) if particles else set(range(num_particles))
    if isinstance(force, mm.CustomCVForce):
        particles = set()
        for k in range(force.getNumCollectiveVariables()):
            variable_particles = _force_particles(force.getCollectiveVariable(k), num_particles)
            if variable_particles is None:
                return None
            particles |= variable_particles
        return particles
    return None


def nonbonded_force_is_solute_inert(system, force, solute_indices, positions, num_conformations=3,
                                    displacement=0.02, seed=2):
    """
    Checks whether a CustomNonbondedForce has no interaction involving the solute.

    True if its interaction groups never contain a solute atom, or if a copy restricted to the
    pairs with at least one solute atom gives zero energy and forces at the given configuration and
    at randomly displaced solute geometries (Reference platform).

    :param system:
        The OpenMM system.
    :param force:
        The CustomNonbondedForce (or its index in the system).
    :param solute_indices:
        The solute's particle indices.
    :param positions:
        The positions of all particles in nm, shape (number of particles, 3).

    :return:
        A tuple (inert, energy, largest force component) of the first failing evaluation.
    """

    if isinstance(force, int):
        force = system.getForce(force)
    force = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force))
    solute = set(solute_indices)
    everyone = set(range(system.getNumParticles()))

    groups = [force.getInteractionGroupParameters(k) for k in range(force.getNumInteractionGroups())]
    if not groups:
        new_groups = [(solute, everyone)]
    else:
        new_groups = []
        for set1, set2 in groups:
            set1, set2 = set(set1), set(set2)
            if set1 & solute and set2:
                new_groups.append((set1 & solute, set2))
            if set2 & solute and set1 - solute:
                new_groups.append((set1 - solute, set2 & solute))
    if not new_groups:
        return True, 0.0, 0.0

    # Interaction groups cannot be removed: reuse the existing ones (emptied) and add more if needed
    for k in range(max(len(new_groups), force.getNumInteractionGroups())):
        set1, set2 = new_groups[k] if k < len(new_groups) else (set(), set())
        if k < force.getNumInteractionGroups():
            force.setInteractionGroupParameters(k, sorted(set1), sorted(set2))
        else:
            force.addInteractionGroup(sorted(set1), sorted(set2))
    force.setForceGroup(0)

    test_system = mm.System()
    for i in range(system.getNumParticles()):
        test_system.addParticle(system.getParticleMass(i))
    test_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    test_system.addForce(force)
    context = mm.Context(test_system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))

    rng = np.random.default_rng(seed)
    coordinates = np.asarray(positions, dtype=float)
    for c in range(num_conformations):
        xyz = coordinates.copy()
        if c > 0:
            xyz[list(solute_indices)] += rng.normal(scale=displacement, size=(len(solute_indices), 3))
        context.setPositions(xyz * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        largest_force = float(np.abs(state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)).max())
        if abs(energy) > 1e-9 or largest_force > 1e-9:
            return False, energy, largest_force
    return True, 0.0, 0.0


def _classify_force(force, system, solute_indices, positions, charges, in_variable=False):
    """
    Classifies one force (or one collective variable of a CustomCVForce) with respect to the
    solute. Returns a dictionary with category, detail and, where relevant, parameter (index of the
    charge parameter of a charge-coupled force) or translation_only (centroid restraints).
    """

    solute = set(solute_indices)
    name = type(force).__name__
    particles = _force_particles(force, system.getNumParticles())

    if particles is not None and not particles & solute:
        return {'category': 'environment', 'detail': 'no solute atom involved'}

    if isinstance(force, mm.CustomNonbondedForce):
        inert, energy, largest_force = nonbonded_force_is_solute_inert(system, force, solute_indices, positions)
        if inert:
            return {'category': 'environment' if in_variable else 'solute-inert',
                    'detail': 'no interaction involving the solute'}
        return {'category': 'unsupported',
                'detail': f'interacts with the solute (energy {energy:.4g} kJ/mol, largest force '
                          f'{largest_force:.4g} kJ/mol/nm) and cannot be switched off by the driver'}

    if isinstance(force, mm.CustomCVForce):
        parts = [(force.getCollectiveVariableName(k),
                  _classify_force(force.getCollectiveVariable(k), system, solute_indices, positions, charges,
                                  in_variable=True))
                 for k in range(force.getNumCollectiveVariables())]
        for variable, part in parts:
            if part['category'] == 'unsupported':
                return {'category': 'unsupported', 'detail': f'collective variable "{variable}": {part["detail"]}'}
        categories = {part['category'] for _, part in parts}
        if categories == {'environment'}:
            return {'category': 'environment', 'detail': 'no solute atom involved'}
        if categories == {'intramolecular'}:
            return {'category': 'intramolecular',
                    'detail': 'collective variables within the solute; copied into the vacuum leg'}
        if particles is not None and particles <= solute:
            variable = [v for v, part in parts if part['category'] != 'intramolecular'][0]
            return {'category': 'unsupported',
                    'detail': f'acts on the solute only, but collective variable "{variable}" cannot be copied '
                              'into the vacuum leg'}
        return {'category': 'restraint',
                'detail': 'collective variables involving solute and environment atoms, kept in every state'}

    if isinstance(force, mm.RMSDForce):
        if particles <= solute:
            return {'category': 'intramolecular',
                    'detail': 'RMSD of solute atoms after superposition; copied into the vacuum leg'}
        return {'category': 'restraint', 'detail': 'RMSD over solute and environment atoms, kept in every state'}

    if name in _STANDARD_BONDED or name in _CUSTOM_BONDED:
        if in_variable:
            # The collective variable is the sum of all its terms: it acts within the solute only if
            # every term does
            if particles <= solute:
                return {'category': 'intramolecular', 'detail': 'terms within the solute; copied into the vacuum leg'}
            return {'category': 'restraint',
                    'detail': 'terms involving solute and environment atoms, kept in every state'}
        get, num_terms, _, _ = _custom_bonded_terms(force)
        mixed = [i for i in range(num_terms) if 0 < len(set(get(i)[0]) & solute) < len(get(i)[0])]
        if mixed:
            return {'category': 'restraint',
                    'detail': f'{len(mixed)} term(s) between solute and environment atoms, kept in every state'}
        return {'category': 'intramolecular', 'detail': 'terms within the solute; copied into the vacuum leg'}

    if isinstance(force, mm.CustomExternalForce):
        if not in_variable:
            particle_terms = {force.getParticleParameters(k)[0]: k for k in range(force.getNumParticles())}
            if solute <= set(particle_terms):
                for p in range(force.getNumPerParticleParameters()):
                    values = np.array([force.getParticleParameters(particle_terms[a])[1][p] for a in solute_indices])
                    if np.allclose(values, charges, rtol=0.0, atol=1e-6):
                        return {'category': 'charge-coupled', 'parameter': p,
                                'detail': f'per-particle parameter "{force.getPerParticleParameterName(p)}" equals '
                                          f'the solute charges ({force.getEnergyFunction()}); scaled with the charges'}
        return {'category': 'restraint', 'detail': 'external potential on solute atoms, kept in every state'}

    if isinstance(force, mm.CustomCentroidBondForce):
        groups = [set(force.getGroupParameters(g)[0]) for g in range(force.getNumGroups())]
        whole = all(group == solute or not group & solute for group in groups)
        return {'category': 'restraint', 'translation_only': whole,
                'detail': 'centroid of the whole solute relative to its environment, kept in every state' if whole
                          else 'involves part of the solute (orientation or conformation), kept in every state'}

    if isinstance(force, mm.NonbondedForce):
        return {'category': 'unsupported',
                'detail': 'a NonbondedForce inside a collective variable is not reached by the alchemical scaling'}

    return {'category': 'unsupported', 'detail': 'acts on solute atoms and is not handled by the driver'}


def classify_forces(system, solute_indices, positions):
    """
    Classifies every force of a pre-built system with respect to the alchemical solute.

    Categories:
        - 'handled': NonbondedForce (scaled by the driver) and standard bonded forces.
        - 'global': barostats, thermostats, center-of-mass motion removal.
        - 'environment': no solute atom involved.
        - 'intramolecular': custom bonded terms, RMSD restraints and collective-variable forces
          within the solute only; copied into the vacuum leg, so they cancel.
        - 'solute-inert': CustomNonbondedForce without any interaction involving the solute.
        - 'charge-coupled': CustomExternalForce with a per-particle parameter equal to the solute's
          NonbondedForce charges (an applied field); it must be scaled with the charges.
        - 'restraint': other terms between the solute and its environment (restraints, external
          potentials, collective variables); kept unchanged in every lambda state of the solvated leg.
        - 'unsupported': interactions between the solute and its environment that the driver cannot
          switch off (e.g. Lennard-Jones terms in a CustomNonbondedForce, also inside a CustomCVForce),
          collective-variable forces on the solute alone that cannot be copied into the vacuum leg,
          and unknown forces acting on the solute.

    :param system:
        The OpenMM system.
    :param solute_indices:
        The solute's particle indices.
    :param positions:
        The positions of all particles in nm, shape (number of particles, 3).

    :return:
        A list with one dictionary per force: index, type, name, category, detail and, for
        'charge-coupled' forces, the index of the charge parameter ('parameter').
    """

    handled = (mm.NonbondedForce, mm.HarmonicBondForce, mm.HarmonicAngleForce, mm.PeriodicTorsionForce,
               mm.RBTorsionForce, mm.CMAPTorsionForce)
    global_forces = (mm.CMMotionRemover, mm.MonteCarloBarostat, mm.MonteCarloAnisotropicBarostat,
                     mm.MonteCarloMembraneBarostat, mm.AndersenThermostat)
    if hasattr(mm, 'MonteCarloFlexibleBarostat'):
        global_forces += (mm.MonteCarloFlexibleBarostat,)

    nonbonded = _get_nonbonded_force(system)
    charges = np.array([_value(nonbonded.getParticleParameters(a)[0], unit.elementary_charge)
                        for a in solute_indices])
    rows = []

    for index, force in enumerate(system.getForces()):
        row = {'index': index, 'type': type(force).__name__, 'name': force.getName()}
        rows.append(row)

        if isinstance(force, global_forces):
            row.update(category='global', detail='')
        elif isinstance(force, handled):
            row.update(category='handled', detail='scaled by the driver' if isinstance(force, mm.NonbondedForce)
                       else 'bonded terms; the solute terms are copied into the vacuum leg')
        else:
            row.update(_classify_force(force, system, solute_indices, positions, charges))

    return rows


def scale_charge_coupled_force(system, solute_indices, force_index, parameter_index, factor):
    """
    Multiplies the charge parameter of the solute atoms in a charge-coupled CustomExternalForce.

    :param system:
        The OpenMM system.
    :param solute_indices:
        The solute's particle indices.
    :param force_index:
        The index of the CustomExternalForce.
    :param parameter_index:
        The index of its per-particle charge parameter.
    :param factor:
        The scaling factor (the same as for the solute's NonbondedForce charges).
    """

    force = system.getForce(force_index)
    solute = set(solute_indices)
    for k in range(force.getNumParticles()):
        particle, parameters = force.getParticleParameters(k)
        if particle in solute:
            parameters = list(parameters)
            parameters[parameter_index] *= factor
            force.setParticleParameters(k, particle, parameters)

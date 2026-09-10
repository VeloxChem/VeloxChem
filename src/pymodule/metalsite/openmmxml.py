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

"""
Writing the fitted metal site out as an OpenMM force field XML.

create_enzyme_system puts the fitted terms onto one System, by index. That
system describes one topology and nothing else: it cannot be solvated,
extended or handed on, because everything downstream of it builds its own
system and throws the edits away. The same parameters written as a force
field XML can be loaded beside the protein force field and used to build
whatever system the caller wants.

What that costs is a restructured topology. A metal-ligand bond has to be
a real bond for OpenMM to write a bonded term across it, and once it is,
the coordinating residues no longer match the protein force field's own
templates. So the sidechains the active site was cut from -- from CB
outward, exactly where extract_active_site truncates -- are moved into one
residue of their own together with the metals, and that residue is what
this module writes a template for. The backbone stays where it is, which
is what keeps the chain intact, the residue numbering readable and a PDB
round trip lossless.

Two templates would then be ambiguous, since OpenMM matches a residue to a
template on the bond graph alone and the stubs left behind are all
isomorphic. The XML therefore carries an InitializationScript that
registers a template matcher keyed on chain and residue id, which is
consulted before the graph search: the file stays self-contained and no
caller has to pass a residueTemplates mapping. Those ids are what ties
the file to a topology, so a PDB of it has to be written with
keepIds=True -- PDBFile numbers the residues from one otherwise, and the
file can then no longer load the topology it was written for.

The fitted metal terms are added by a Script rather than written as
<Bond> and <Angle> entries, and that is forced by how amber14 is keyed.
Telling two metal-ligand bonds of the same kind apart needs an atom type
per coordinating atom -- five of the metal-oxygen bonds of a two-zinc site
share one pair of atom classes while their fitted equilibria span 0.196 to
0.279 nm -- but every bonded term in amber14 names its atoms by type
rather than by class, so an atom given a type of its own inherits none of
the protein force field's covalent terms and all of them would have to be
written out again. Leaving every atom on the type the protein force field
gave it keeps those terms exactly as they were, and the metal terms are
put on by name instead. The script runs on every createSystem, so a
solvated or rebuilt system gets them too.
"""

import xml.etree.ElementTree as ET
import string
import sys

import numpy as np

from ..errorhandler import assert_msg_critical
from . import core
from .printing import stream

try:
    import openmm.app as mmapp
    from openmm.app.internal import compiled
except ImportError:
    pass

# The residue the metals and the coordinating sidechains are moved into.
SITE_RESIDUE_NAME = 'MSI'
# Not a standard PDB residue name, which is what makes PDBFile write a
# CONECT record for every bond the residue takes part in and read it back:
# a metal-to-protein bond is only read back when the protein atom sits in a
# residue the reader does not think it already knows the bonding of.
SITE_RESIDUE_ID = '1'

# What a run writes beside the system.
SITE_XML_FILE = 'metal_site.xml'
SITE_TOPOLOGY_FILE = 'metal_site.pdb'

# One letter per source residue, appended to the atom name it came in with
# so that the merged residue holds no two atoms of the same name. Kept to a
# single character because a PDB atom name is four columns wide.
SLOT_LETTERS = string.ascii_uppercase


# ----------------------------------------------------------------------
# restructuring
# ----------------------------------------------------------------------


def restructure_topology(topology,
                         positions,
                         active_site,
                         forcefield,
                         site_residue_name=SITE_RESIDUE_NAME,
                         ostream=None):
    """
    Moves the active site into a residue of its own and bonds the metals.

    The atoms the active site was cut from -- everything it holds except
    the capping hydrogens, which stand for alpha carbons that stay behind
    -- are taken out of their residues and put into one new residue at the
    end of a chain of its own. The backbone of every coordinating residue
    stays exactly where it was, so the chain keeps its residue list intact
    and nothing downstream of this sees a gap where a residue used to be.

    The metal-ligand bonds are added as real topology bonds, which is what
    makes OpenMM write bonded terms across them -- and, with them, the
    1-2 and 1-3 exclusions and the 1-4 scaling a bonded metal model should
    have and the injected system never had.

    :param topology:
        The protonated OpenMM topology the active site was extracted from.
    :param positions:
        Its positions as an (N, 3) array in Angstrom.
    :param active_site:
        The active site, for the map back to the topology.
    :param forcefield:
        The force field generator, for the metal bonds to add.

    :return:
        A dictionary holding the new topology and positions, the map from
        old atom index to new and back, where each active site atom ended
        up and under what name, and the new residues the sidechains were
        taken out of, keyed by the index they had in the old topology.
    """

    ostream = stream(ostream)

    assert_msg_critical('openmm.app' in sys.modules,
                        'restructure_topology: openmm is required')

    atom_map = active_site['atom_map']
    caps = set(active_site['cap_indices'])
    site_indices = [index for index in sorted(atom_map) if index not in caps]
    moved = {atom_map[index]: index for index in site_indices}

    assert_msg_critical(
        len(moved) == len(site_indices), 'restructure_topology: two active '
        'site atoms map to one topology atom')

    atoms = list(topology.atoms())

    assert_msg_critical(
        max(moved) < len(atoms), 'restructure_topology: the atom map does '
        'not index this topology. Extract the active site from the topology '
        'the system is built from; prepare_protein renumbers the atoms, so '
        'it must run first.')

    names = _site_atom_names(atoms, atom_map, site_indices)

    new_topology = mmapp.Topology()
    box_vectors = topology.getPeriodicBoxVectors()
    if box_vectors is not None:
        new_topology.setPeriodicBoxVectors(box_vectors)

    index_map = {}
    stub_residues = {}

    for chain in topology.chains():
        kept = []
        for residue in chain.residues():
            original = list(residue.atoms())
            remaining = [atom for atom in original if atom.index not in moved]
            if remaining:
                kept.append((residue, remaining, len(remaining) != len(original)))

        if not kept:
            continue

        new_chain = new_topology.addChain(chain.id)
        for residue, remaining, is_stub in kept:
            new_residue = new_topology.addResidue(residue.name, new_chain,
                                                  residue.id,
                                                  residue.insertionCode)
            if is_stub:
                stub_residues[residue.index] = new_residue
            for atom in remaining:
                index_map[atom.index] = new_topology.addAtom(
                    atom.name, atom.element, new_residue, atom.id).index

    site_chain = new_topology.addChain(_free_chain_id(topology))
    site_residue = new_topology.addResidue(site_residue_name, site_chain,
                                           SITE_RESIDUE_ID)

    site_atom_index = {}
    for index in site_indices:
        atom = atoms[atom_map[index]]
        new_atom = new_topology.addAtom(names[index], atom.element,
                                        site_residue)
        index_map[atom.index] = new_atom.index
        site_atom_index[index] = new_atom.index

    new_atoms = list(new_topology.atoms())
    bonded = set()

    for first, second in topology.bonds():
        pair = (index_map[first.index], index_map[second.index])
        new_topology.addBond(new_atoms[pair[0]], new_atoms[pair[1]])
        bonded.add(frozenset(pair))

    metal_bonds = metal_bond_keys(forcefield, active_site)
    for key in metal_bonds:
        pair = (site_atom_index[key[0]], site_atom_index[key[1]])
        if frozenset(pair) not in bonded:
            new_topology.addBond(new_atoms[pair[0]], new_atoms[pair[1]])
            bonded.add(frozenset(pair))

    new_positions = np.zeros((new_topology.getNumAtoms(), 3))
    for old_index, new_index in index_map.items():
        new_positions[new_index] = positions[old_index]

    ostream.print_info(
        f'Moved {len(site_indices)} atoms out of {len(stub_residues)} '
        f'residues into residue {site_residue_name}, and bonded '
        f'{len(metal_bonds)} metal-ligand contacts.')
    ostream.flush()

    return {
        'topology': new_topology,
        'positions': new_positions,
        'atom_index_map': index_map,
        'site_atom_index': site_atom_index,
        'site_atom_name': names,
        'site_indices': site_indices,
        'atom_source': {new: old for old, new in index_map.items()},
        'site_residue': site_residue,
        'stub_residues': stub_residues,
    }


def _site_atom_names(atoms, atom_map, site_indices):
    """
    Names the atoms of the merged residue, uniquely and in four characters.

    An atom keeps the name it came in with, with one letter appended for
    the residue it came out of, so that a name still says what the atom was.
    Anything that would not fit a PDB atom name field, or that would repeat
    a name already taken, falls back to the element and a running number.

    :return:
        The name for each active site index.
    """

    slots = {}
    for index in site_indices:
        residue = atoms[atom_map[index]].residue
        if residue.index not in slots:
            assert_msg_critical(
                len(slots) < len(SLOT_LETTERS),
                '_site_atom_names: the active site holds more than '
                f'{len(SLOT_LETTERS)} residues')
            slots[residue.index] = SLOT_LETTERS[len(slots)]

    names = {}
    taken = set()
    for index in site_indices:
        atom = atoms[atom_map[index]]
        name = f'{atom.name}{slots[atom.residue.index]}'
        if len(name) > 4 or name in taken:
            symbol = atom.element.symbol if atom.element is not None else 'X'
            counter = 1
            name = f'{symbol}{counter}'
            while name in taken or len(name) > 4:
                counter += 1
                name = f'{symbol}{counter}'
        taken.add(name)
        names[index] = name

    return names


def _free_chain_id(topology):
    """
    Returns a chain id the topology does not already use.
    """

    used = {chain.id for chain in topology.chains()}
    for candidate in SLOT_LETTERS + string.digits:
        if candidate not in used:
            return candidate

    return 'MS'


def metal_bond_keys(forcefield, active_site):
    """
    The metal bonds that become bonds of the restructured topology.

    A term that names a capping hydrogen is left out the way
    create_enzyme_system leaves it out: the cap stands for an alpha carbon
    the protein force field parameterizes itself.

    :return:
        The bond keys, as active site index pairs.
    """

    caps = set(active_site['cap_indices'])
    bonds, _ = core.get_metal_keys(forcefield, active_site)

    return [key for key in bonds if not caps & set(key)]


def metal_term_keys(forcefield, active_site):
    """
    Every metal term that reaches the restructured topology, by kind.

    :return:
        The tuple of the bond, angle and improper keys, with the terms
        naming a capping hydrogen already dropped.
    """

    caps = set(active_site['cap_indices'])
    bonds, angles = core.get_metal_keys(forcefield, active_site)
    impropers = core.get_metal_impropers(forcefield, active_site)

    return ([key for key in bonds if not caps & set(key)],
            [key for key in angles if not caps & set(key)],
            [key for key in impropers if not caps & set(key)])


# ----------------------------------------------------------------------
# what the protein force field says
# ----------------------------------------------------------------------


def protein_atom_parameters(topology, forcefield_files, ostream=None):
    """
    Reads the atom type and charge the protein force field gives every atom.

    Asked of the topology as it stands, before any metal bond is added, so
    that every residue still matches a template of its own. Every atom
    keeps the type read here, which is what leaves the covalent terms of
    the coordinating residues exactly as the protein force field wrote
    them; only the charges change, and only where the fit covers them.

    :param topology:
        The protonated topology, without the metal bonds.
    :param forcefield_files:
        The OpenMM force field files for the protein.

    :return:
        The tuple of the loaded ForceField and, per topology atom index,
        the atom type and charge it was given.
    """

    ostream = stream(ostream)

    assert_msg_critical('openmm.app' in sys.modules,
                        'protein_atom_parameters: openmm is required')

    protein_ff = mmapp.ForceField(*forcefield_files)
    bonded_to_atom = protein_ff._buildBondedToAtomList(topology)
    templates = protein_ff.getMatchingTemplates(topology)

    parameters = {}
    for residue, template in zip(topology.residues(), templates):
        assert_msg_critical(
            template is not None, 'protein_atom_parameters: the protein '
            f'force field has no template for {residue.name}{residue.id}. '
            'The topology must be the one the active site was extracted '
            'from, before any metal bond was added to it.')

        matches = compiled.matchResidueToTemplate(residue, template,
                                                  bonded_to_atom, False, False)

        assert_msg_critical(
            matches is not None, 'protein_atom_parameters: the template '
            f'{template.name} stopped matching {residue.name}{residue.id}')

        for atom, index in zip(residue.atoms(), matches):
            template_atom = template.atoms[index]
            parameters[atom.index] = {
                'type': template_atom.type,
                'charge': template_atom.parameters['charge'],
            }

    return protein_ff, parameters


# ----------------------------------------------------------------------
# templates
# ----------------------------------------------------------------------


def build_templates(topology,
                    active_site,
                    restructured,
                    partial_charges,
                    protein_parameters,
                    ostream=None):
    """
    Builds the residue templates the restructured topology needs.

    One for the merged site residue, carrying the fitted charges, and one
    for every backbone stub the sidechains were taken out of, carrying the
    protein force field's own charges with the shift that gives the
    coordination region its charge back. The shift is the one
    redistribute_backbone_charges applies to a built system, computed
    through the same function, so the two paths cannot drift apart on it.

    :param topology:
        The original protonated topology.
    :param active_site:
        The active site.
    :param restructured:
        What restructure_topology returned.
    :param partial_charges:
        The active site charges, capping hydrogens included.
    :param protein_parameters:
        What protein_atom_parameters returned for the original topology.

    :return:
        The tuple of the template list and the shift each atom the active
        site does not cover took on.
    """

    ostream = stream(ostream)

    charges = core.redistribute_cap_charges(active_site, partial_charges)
    correction = core.backbone_charge_shift(
        lambda index: protein_parameters[index]['charge'], topology,
        active_site, charges)
    shift = correction['shift']

    atom_map = active_site['atom_map']
    site_atom_index = restructured['site_atom_index']
    source = restructured['atom_source']

    types = {}
    new_charges = {}

    for index in restructured['site_indices']:
        new_index = site_atom_index[index]
        types[new_index] = protein_parameters[atom_map[index]]['type']
        new_charges[new_index] = float(charges[index])

    for stub in restructured['stub_residues'].values():
        for atom in stub.atoms():
            protein = protein_parameters[source[atom.index]]
            types[atom.index] = protein['type']
            new_charges[atom.index] = protein['charge'] + shift

    bonded = {}
    for first, second in restructured['topology'].bonds():
        bonded.setdefault(first.index, set()).add(second.index)
        bonded.setdefault(second.index, set()).add(first.index)

    site_residue = restructured['site_residue']
    templates = [
        _template_of_residue(
            site_residue, bonded, types, new_charges,
            f'{site_residue.name}-{site_residue.chain.id}{site_residue.id}')
    ]

    for stub in restructured['stub_residues'].values():
        templates.append(
            _template_of_residue(
                stub, bonded, types, new_charges,
                f'{stub.name}-{stub.chain.id}{stub.id}-{site_residue.name}'))

    ostream.print_info(
        f'Built {len(templates)} residue templates; the coordination region '
        f'gives each of the {len(correction["uncovered"])} atoms the active '
        f'site does not cover {shift:+.4f} e.')
    ostream.flush()

    return templates, shift


def _template_of_residue(residue, bonded, types, charges, name):
    """
    Turns one residue of the restructured topology into a template.

    The bonds and the external bonds are read off the topology rather than
    worked out from what the residue is meant to be, so a template says
    exactly what OpenMM will be matching it against.

    :return:
        The template, as its name, the residue it is for, its atoms with
        their types and charges, its internal bonds and its external ones.
    """

    members = list(residue.atoms())
    member_indices = {atom.index for atom in members}
    name_of = {atom.index: atom.name for atom in members}

    atoms = []
    for atom in members:
        assert_msg_critical(
            atom.index in types, '_template_of_residue: no atom type for '
            f'{atom.name} of {residue.name}{residue.id}')
        atoms.append((atom.name, types[atom.index], charges[atom.index]))

    internal = []
    external = []
    for atom in members:
        for other in sorted(bonded.get(atom.index, ())):
            if other in member_indices:
                if other > atom.index:
                    internal.append((atom.name, name_of[other]))
            else:
                external.append(atom.name)

    return {
        'name': name,
        'chain': residue.chain.id,
        'residue_id': residue.id,
        'residue_name': residue.name,
        'atoms': atoms,
        'bonds': internal,
        'external': external,
    }


# ----------------------------------------------------------------------
# the xml
# ----------------------------------------------------------------------

# Registered while the file is loaded and consulted before the graph
# search, so the file needs no residueTemplates mapping from its caller.
# OpenMM matches a residue to a template on the bond graph alone, and the
# backbone stubs are all isomorphic to one another, so without this the
# first two of them are a "Multiple non-identical matching templates"
# error. Keyed on chain and residue id, which a rebuilt topology keeps and
# an object-keyed mapping would not.
MATCHER_SCRIPT = '\n'.join([
    'class _MetalSiteTemplateMatcher(object):',
    '    def __init__(self, assignment):',
    '        self.assignment = assignment',
    '',
    '    def __call__(self, forcefield, residue, bondedToAtom,',
    '                 ignoreExternalBonds, ignoreExtraParticles):',
    '        name = self.assignment.get((residue.chain.id, residue.id))',
    '        if name is None:',
    '            return None',
    '        return forcefield._templates[name]',
    '',
    'self.registerTemplateMatcher(_MetalSiteTemplateMatcher({assignment}))',
    '',
])

# The fitted metal terms. Written here rather than as <Bond>, <Angle> and
# <Improper> entries because every bonded term in amber14 names its atoms
# by type, so telling one metal-ligand bond from another -- which needs a
# type per coordinating atom -- would cost every covalent term those atoms
# take part in as well. See the module docstring.
#
# Impropers would have to be here in any case: OpenMM builds one out of a
# residue template only for a central atom and three atoms bonded to it,
# and the improper a bidentate carboxylate gets names the metal, which its
# central carbon is not bonded to.
#
# The proper torsions the protein force field writes across a metal bond
# are zeroed first, before anything is added. They come from its wildcard
# entries -- X-CT-S-X for a bound cysteine, X-CR-NB-X and X-CV-NB-X for a
# bound histidine -- and reach the metal only because the bond is now a
# real one. Two of the three keep the metal in the plane of an imidazole,
# which is what the fitted improper on the coordinating nitrogen already
# does, at four times the barrier and without anything having fitted it.
# Nothing about the metal should come from a wildcard, so they go, and the
# fit is left to describe the coordination on its own.
#
# Executed with the locals of createSystem as its namespace, so `sys` is
# the System being built and `topology` the topology it is being built
# for. Nothing imported by forcefield.py is in scope, which is why the
# forces are found by class name.
TERM_SCRIPT = '\n'.join([
    '_ms_site = {site}',
    '_ms_expected = {expected}',
    '_ms_drop_torsions = {drop_torsions}',
    '_ms_bonds = {bonds}',
    '_ms_angles = {angles}',
    '_ms_impropers = {impropers}',
    '',
    '_ms_index = {{}}',
    'for _ms_residue in topology.residues():',
    '    if (_ms_residue.chain.id, _ms_residue.id) == _ms_site:',
    '        for _ms_atom in _ms_residue.atoms():',
    '            _ms_index[_ms_atom.name] = _ms_atom.index',
    '',
    'if len(_ms_index) != _ms_expected:',
    '    raise ValueError(',
    '        "metal site force field: expected %d atoms in residue %s of "',
    '        "chain %s, found %d. The metal site is keyed on chain and "',
    '        "residue id, so it cannot be renamed or renumbered after the "',
    '        "force field was written."',
    '        % (_ms_expected, _ms_site[1], _ms_site[0], len(_ms_index)))',
    '',
    '_ms_forces = {{}}',
    'for _ms_force in sys.getForces():',
    '    _ms_forces[_ms_force.__class__.__name__] = _ms_force',
    '',
    'for _ms_kind, _ms_terms in (("HarmonicBondForce", _ms_bonds),',
    '                            ("HarmonicAngleForce", _ms_angles),',
    '                            ("PeriodicTorsionForce", _ms_impropers)):',
    '    if _ms_terms and _ms_kind not in _ms_forces:',
    '        raise ValueError(',
    '            "metal site force field: the system has no %s to add the "',
    '            "fitted metal terms to" % _ms_kind)',
    '',
    '_ms_pairs = set()',
    'for _ms_names, _ms_length, _ms_k in _ms_bonds:',
    '    _ms_pairs.add(frozenset((_ms_index[_ms_names[0]],',
    '                             _ms_index[_ms_names[1]])))',
    '',
    'if _ms_drop_torsions and "PeriodicTorsionForce" in _ms_forces:',
    '    _ms_force = _ms_forces["PeriodicTorsionForce"]',
    '    for _ms_i in range(_ms_force.getNumTorsions()):',
    '        _ms_p = _ms_force.getTorsionParameters(_ms_i)',
    '        if (frozenset(_ms_p[0:2]) in _ms_pairs',
    '                or frozenset(_ms_p[1:3]) in _ms_pairs',
    '                or frozenset(_ms_p[2:4]) in _ms_pairs):',
    '            _ms_force.setTorsionParameters(',
    '                _ms_i, _ms_p[0], _ms_p[1], _ms_p[2], _ms_p[3],',
    '                _ms_p[4], _ms_p[5], 0.0)',
    '',
    'for _ms_names, _ms_length, _ms_k in _ms_bonds:',
    '    _ms_forces["HarmonicBondForce"].addBond(',
    '        _ms_index[_ms_names[0]], _ms_index[_ms_names[1]],',
    '        _ms_length, _ms_k)',
    '',
    'for _ms_names, _ms_angle, _ms_k in _ms_angles:',
    '    _ms_forces["HarmonicAngleForce"].addAngle(',
    '        _ms_index[_ms_names[0]], _ms_index[_ms_names[1]],',
    '        _ms_index[_ms_names[2]], _ms_angle, _ms_k)',
    '',
    'for _ms_names, _ms_period, _ms_phase, _ms_barrier in _ms_impropers:',
    '    _ms_forces["PeriodicTorsionForce"].addTorsion(',
    '        _ms_index[_ms_names[0]], _ms_index[_ms_names[1]],',
    '        _ms_index[_ms_names[2]], _ms_index[_ms_names[3]],',
    '        _ms_period, _ms_phase, _ms_barrier)',
    '',
])


def forcefield_xml(active_site,
                   forcefield,
                   restructured,
                   templates,
                   forcefield_files=(),
                   drop_torsions_across_metal_bonds=True,
                   ostream=None):
    """
    Writes the force field XML for the restructured topology.

    :param active_site:
        The active site.
    :param forcefield:
        The force field generator carrying the fitted metal terms.
    :param restructured:
        What restructure_topology returned.
    :param templates:
        The residue templates, from build_templates.
    :param forcefield_files:
        The protein force field files the templates name atom types of,
        recorded in a comment so that a reader knows what the file has to
        be loaded beside.
    :param drop_torsions_across_metal_bonds:
        Whether to zero the wildcard proper torsions the protein force
        field writes across a metal bond once it is a real one. See the
        comment on TERM_SCRIPT for why they should go.

    :return:
        The XML as a string.
    """

    ostream = stream(ostream)

    root = ET.Element('ForceField')

    # every atom keeps the type the protein force field gave it, so the
    # file says nothing on its own and has to be loaded beside the same
    # files it was written against
    if forcefield_files:
        root.append(
            ET.Comment(' Load beside: ' + ', '.join(forcefield_files) + ' '))

    residues = ET.SubElement(root, 'Residues')
    for template in templates:
        element = ET.SubElement(residues, 'Residue', name=template['name'])
        for name, atom_type, charge in template['atoms']:
            ET.SubElement(element,
                          'Atom',
                          name=name,
                          type=atom_type,
                          charge=repr(float(charge)))
        for first, second in template['bonds']:
            ET.SubElement(element, 'Bond', atomName1=first, atomName2=second)
        for name in template['external']:
            ET.SubElement(element, 'ExternalBond', atomName=name)

    assignment = {(template['chain'], template['residue_id']): template['name']
                  for template in templates}
    initialization = ET.SubElement(root, 'InitializationScript')
    initialization.text = '\n' + MATCHER_SCRIPT.format(
        assignment=repr(assignment))

    bonds, angles, impropers = metal_term_keys(forcefield, active_site)
    names = restructured['site_atom_name']
    site_residue = restructured['site_residue']

    bond_terms = [((names[key[0]], names[key[1]]),
                   float(forcefield.bonds[key]['equilibrium']),
                   float(forcefield.bonds[key]['force_constant']))
                  for key in bonds]

    angle_terms = [((names[key[0]], names[key[1]], names[key[2]]),
                    float(np.radians(forcefield.angles[key]['equilibrium'])),
                    float(forcefield.angles[key]['force_constant']))
                   for key in angles]

    # the central atom of an improper key goes third, which is where
    # OpenMM puts it when it builds one out of a residue template itself
    improper_terms = [
        ((names[key[1]], names[key[2]], names[key[0]], names[key[3]]),
         int(forcefield.impropers[key]['periodicity']),
         float(np.radians(forcefield.impropers[key]['phase'])),
         float(forcefield.impropers[key]['barrier'])) for key in impropers
    ]

    script = ET.SubElement(root, 'Script')
    script.text = '\n' + TERM_SCRIPT.format(
        site=repr((site_residue.chain.id, site_residue.id)),
        expected=len(restructured['site_indices']),
        drop_torsions=repr(bool(drop_torsions_across_metal_bonds)),
        bonds=repr(bond_terms),
        angles=repr(angle_terms),
        impropers=repr(improper_terms))

    ET.indent(root, space=' ')

    ostream.print_info(
        f'Wrote {len(templates)} residue templates, {len(bond_terms)} metal '
        f'bonds, {len(angle_terms)} metal angles and {len(improper_terms)} '
        'metal impropers.')
    ostream.flush()

    return ET.tostring(root, encoding='unicode') + '\n'


# ----------------------------------------------------------------------
# the whole of it
# ----------------------------------------------------------------------


def create_enzyme_forcefield(topology,
                             positions,
                             active_site,
                             forcefield,
                             partial_charges=None,
                             forcefield_files=('amber14-all.xml',
                                               'amber14/tip3pfb.xml'),
                             site_residue_name=SITE_RESIDUE_NAME,
                             drop_torsions_across_metal_bonds=True,
                             ostream=None):
    """
    Writes the fitted metal site as an OpenMM force field XML.

    The counterpart of create_enzyme_system, and what it should be
    preferred to: that one puts the fitted terms onto one System and can
    describe nothing else, while this returns a force field and the
    topology it is for, which can be solvated, extended and rebuilt as
    often as the caller likes.

    The metal-ligand bonds become real bonds of the topology it returns,
    so the system built from it carries the 1-2 and 1-3 exclusions and the
    1-4 scaling a bonded metal model should have. The injected system has
    none of those, and that is the one place the two are meant to differ.

    :param topology:
        The protonated OpenMM topology of the whole enzyme.
    :param positions:
        Its positions as an (N, 3) array in Angstrom.
    :param active_site:
        The active site, for the map back to the topology.
    :param forcefield:
        The force field generator carrying the fitted metal parameters.
    :param partial_charges:
        The charges fitted on the active site. D4 charges are used when
        none are given, the way create_enzyme_system falls back, so that
        the file carries the same charges as the force field built beside
        it rather than the protein force field's own.
    :param forcefield_files:
        The OpenMM force field files for the protein. Every atom keeps the
        type these give it, so the file that comes back has to be loaded
        beside the same ones.
    :param site_residue_name:
        The name of the residue the site is moved into.
    :param drop_torsions_across_metal_bonds:
        Whether to zero the wildcard proper torsions the protein force
        field writes across a metal bond once it is a real one.

    :return:
        A dictionary holding the XML, the restructured topology and its
        positions, the residue templates and what restructure_topology
        returned.
    """

    ostream = stream(ostream)

    assert_msg_critical('openmm.app' in sys.modules,
                        'create_enzyme_forcefield: openmm is required')

    if partial_charges is None:
        partial_charges = core.d4_charges(active_site, ostream=ostream)

    protein_ff_not_used, protein_parameters = protein_atom_parameters(
        topology, forcefield_files, ostream=ostream)

    restructured = restructure_topology(topology,
                                        positions,
                                        active_site,
                                        forcefield,
                                        site_residue_name=site_residue_name,
                                        ostream=ostream)

    templates, shift = build_templates(topology,
                                       active_site,
                                       restructured,
                                       partial_charges,
                                       protein_parameters,
                                       ostream=ostream)

    xml = forcefield_xml(
        active_site,
        forcefield,
        restructured,
        templates,
        forcefield_files=forcefield_files,
        drop_torsions_across_metal_bonds=drop_torsions_across_metal_bonds,
        ostream=ostream)

    return {
        'xml': xml,
        'topology': restructured['topology'],
        'positions': restructured['positions'],
        'templates': templates,
        'backbone_shift': shift,
        'restructured': restructured,
    }

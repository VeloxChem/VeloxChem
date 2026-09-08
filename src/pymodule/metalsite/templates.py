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
Turning what a run left in a folder into a template.

A template is an active site with the two files it was loaded from added to
it, rather than a dictionary describing one: every matching step then takes
two dictionaries of the same shape, and a template is simply an active site
with no topology behind it.

The rules that live here: which geometry a folder is allowed to be built
from and how the two are told apart, that the bonding comes from the force
field's own bonds rather than a third file, that the caps and beta carbons
are read off the comments annotate_atoms wrote, and that a residue
coordinating no metal is discarded with its charge put back on what stays.
"""

import numpy as np

from ..molecule import Molecule
from ..errorhandler import assert_msg_critical
from . import core
from . import matching
from .printing import stream

# The geometry a run leaves behind, in the order it is looked for. Which of
# them a template is allowed to be built from is what the fallback argument
# of load_template_from_folder decides.
GEOMETRY_KINDS = ('qm_opt', 'mm_opt')


def load_geometry(folder, fallback):
    """
    Reads the active site geometry of a template folder and works out what
    kind of geometry it is.

    :param folder:
        The folder of an earlier run.
    :param fallback:
        The lowest acceptable kind, as given to load_template_from_folder.

    :return:
        The tuple of the molecule and its kind.
    """

    opt_path = folder / core.GEOMETRY_FILE
    mm_path = folder / core.MM_GEOMETRY_FILE

    allowed = GEOMETRY_KINDS[:1 if fallback is None else 2]

    if opt_path.is_file():
        geometry = Molecule.read_xyz_file(str(opt_path))

        if mm_path.is_file() and _same_geometry(
                geometry, Molecule.read_xyz_file(str(mm_path))):
            # build_forcefield writes whatever the active site ended
            # up as under the name of the optimized geometry, so an
            # untouched copy of the MM geometry means the QM optimization
            # never ran
            kind = 'mm_opt'
            assert_msg_critical(
                kind in allowed,
                f'MetalForceFieldManager: {opt_path} is identical to '
                f'{mm_path.name}, so it is the crude MM geometry and no '
                'QM optimization ran. Pass fallback="mm_opt" to use it '
                'anyway.')
        else:
            kind = 'qm_opt'

        return geometry, kind

    assert_msg_critical(
        mm_path.is_file(),
        f'MetalForceFieldManager: neither {opt_path.name} nor '
        f'{mm_path.name} found in {folder}, so there is no geometry to '
        'build a template from.')

    assert_msg_critical(
        'mm_opt' in allowed,
        f'MetalForceFieldManager: {folder} holds only {mm_path.name}, a '
        'geometry relaxed on the crude force field. Pass '
        'fallback="mm_opt" to use it anyway.')

    return Molecule.read_xyz_file(str(mm_path)), 'mm_opt'


def _same_geometry(first, second):
    """
    Whether two molecules hold the same coordinates.

    :param first:
        The first molecule.
    :param second:
        The second molecule.

    :return:
        True when they match to within the precision of an xyz file.
    """

    if first.number_of_atoms() != second.number_of_atoms():
        return False

    return np.allclose(first.get_coordinates_in_angstrom(),
                       second.get_coordinates_in_angstrom(),
                       atol=1.0e-6)


def build(name, forcefield, molecule, kind, folder, metal_elements=None, ostream=None):
    """
    Assembles a template from a loaded force field and geometry.

    The bonding topology is not a separate file: the bonds of the force
    field are the edges of the active site, metal-ligand bonds included,
    since build_forcefield builds them from the connectivity matrix. The
    capping hydrogens and the beta carbons come from the comments
    annotate_atoms wrote, which is how a cap is told apart from the
    hydrogens it is symmetric with.

    What the folder holds is not taken as it stands: the residues that
    coordinate no metal are discarded by _prune_unconnected_residues
    before the template is handed back.

    :param name:
        The name of the template.
    :param forcefield:
        The force field generator loaded from the folder.
    :param molecule:
        The active site geometry.
    :param kind:
        The geometry kind.
    :param folder:
        The folder the template came from.

    :return:
        The template dictionary.
    """

    # the same check a run makes of its own force field, which is why it
    # is the core's and not a second copy of it here. A template is an
    # active site with no topology behind it, so {'molecule': ...} is all
    # the adapter it needs.
    core._check_forcefield(forcefield, {'molecule': molecule}, source=folder)

    labels = list(molecule.get_labels())

    metal_indices = [
        index for index, label in enumerate(labels)
        if label in metal_elements
    ]

    assert_msg_critical(
        len(metal_indices) > 0,
        f'MetalForceFieldManager: the template of {folder} holds no metal '
        f'center. Recognized elements: {metal_elements}')

    charges = np.array(
        [atom['charge'] for atom in forcefield.atoms.values()])

    # The truncation records what it made in the comments, and those are
    # what a template reads. The capping hydrogens do sit at exactly zero
    # charge once redistribute_cap_charges has run, but that is the
    # correction doing its job rather than a marker, and it says nothing
    # in a force field whose charges are all zero.
    marked = {
        role: [
            index for index, atom in forcefield.atoms.items()
            if role in (atom.get('comment', '') or '')
        ]
        for role in (core.BETA_CARBON_COMMENT, core.CAP_COMMENT)
    }

    beta_carbon_indices = marked[core.BETA_CARBON_COMMENT]
    cap_indices = marked[core.CAP_COMMENT]

    assert_msg_critical(
        len(cap_indices) > 0,
        f'MetalForceFieldManager: the force field of {folder} does not '
        'say which of its atoms are capping hydrogens. It was written '
        'before annotate_atoms recorded them, so rebuild it.')

    # A template is an active site with the two files it was loaded from
    # added to it, rather than a dictionary holding a description of one:
    # every matching step then takes two dictionaries of one shape, and a
    # template is simply an active site that has no topology behind it.
    template = {
        'molecule': molecule,
        'metal_indices': metal_indices,
        'cap_indices': cap_indices,
        'beta_carbon_indices': beta_carbon_indices,
        'name': name,
        'forcefield': forcefield,
        'geometry_kind': kind,
        'charges': charges,
        'folder': str(folder),
    }

    described = matching.describe(template, forcefield.bonds.keys())

    return _prune_unconnected_residues(described, ostream=ostream)


def _prune_unconnected_residues(template, ostream=None):
    """
    Discards the residues of a template that coordinate no metal.

    A template exists for the coordination it was fitted for, and a
    residue reaching none of the metals says nothing about that one: it
    is a neighbour the truncation happened to keep. Nothing is ever
    transferred from it -- the metal bonds and angles that are, are keyed
    on the residues that stay -- while every structure compared against
    the template has to hold a residue of its kind, in the right place,
    or map onto nothing at all. Dropping it therefore widens what the
    template matches and changes nothing it is used for.

    The atoms of a discarded residue take their charge with them, which
    leaves what is left off a whole number of electrons.
    _compensate_charges puts that back.

    :param template:
        The described template.

    :return:
        The template with those residues dropped and described again, or
        the object it was given when every residue coordinates a metal.
    """

    coarse = template['coarse_topology']
    unconnected = [
        node for node in matching.residue_nodes(coarse)
        if coarse.degree(node) == 0
    ]

    if not unconnected:
        return template

    discard = sorted(atom for node in unconnected
                     for atom in coarse.nodes[node]['atoms'])

    pruned, shift = _drop_atoms(template, discard, ostream=ostream)
    _print_discarded_residues(template,
                              unconnected,
                              discard,
                              shift,
                              ostream=ostream)

    return matching.describe(pruned, pruned['forcefield'].bonds.keys())


def _drop_atoms(template, discard, ostream=None):
    """
    Removes atoms from a template and renumbers everything left.

    A template is indexed by atom throughout -- the geometry, the
    charges, the marked capping hydrogens and beta carbons, and every key
    of the force field -- so all of them are rewritten together against
    one map from the old index to the new. The force field generator is
    rewritten in place: it is loaded for this template alone and is
    shared with nothing.

    :param template:
        The described template.
    :param discard:
        The indices of the atoms to remove.

    :return:
        The tuple of a new template dictionary, without the topologies
        _describe adds to one, and the charge every remaining atom was
        shifted by.
    """

    molecule = template['molecule']
    labels = molecule.get_labels()
    dropped = set(discard)
    keep = [index for index in range(len(labels)) if index not in dropped]
    index_of = {old: new for new, old in enumerate(keep)}

    def renumbered(indices):
        return [index_of[index] for index in indices if index in index_of]

    forcefield = template['forcefield']
    forcefield.atoms = {
        index_of[index]: atom
        for index, atom in forcefield.atoms.items() if index in index_of
    }

    # a term of a discarded residue is internal to it, since a residue
    # that coordinates no metal shares no bond with the rest of the site;
    # the keys are filtered rather than assumed so anything else that
    # crosses the cut goes with it
    for name in ('bonds', 'angles', 'dihedrals', 'impropers'):
        table = getattr(forcefield, name)
        setattr(
            forcefield, name, {
                tuple(index_of[index] for index in key): parameters
                for key, parameters in table.items()
                if all(index in index_of for index in key)
            })

    geometry = Molecule([labels[index] for index in keep],
                        molecule.get_coordinates_in_angstrom()[keep],
                        'angstrom')
    geometry.set_charge(molecule.get_charge())
    geometry.set_multiplicity(molecule.get_multiplicity())
    forcefield.molecule = geometry

    cap_indices = renumbered(template['cap_indices'])
    charges, shift = compensate_charges(template['charges'][keep],
                                        cap_indices)

    # one array means the charges: what the template reports is what the
    # force field carries, here as everywhere else
    for index, atom in forcefield.atoms.items():
        atom['charge'] = float(charges[index])
    forcefield.partial_charges = charges

    pruned = {
        key: value
        for key, value in template.items()
        if key not in ('fine_topology', 'coarse_topology', 'composition')
    }
    pruned.update({
        'molecule': geometry,
        'metal_indices': renumbered(template['metal_indices']),
        'cap_indices': cap_indices,
        'beta_carbon_indices':
        renumbered(template['beta_carbon_indices']),
        'charges': charges,
    })

    return pruned, shift


def compensate_charges(charges, cap_indices):
    """
    Puts the charge the discarded residues took with them back onto the
    atoms that stay.

    What is left of a template is a whole number of electrons, and the
    atoms that went carried a share of the charge that is not one, so the
    difference from the nearest integer is added in equal parts to every
    remaining atom. The capping hydrogens are left out of it, for the
    reason redistribute_cap_charges empties them in the first place: a
    cap stands for a bond to a protein the site does not hold, and charge
    put on one is taken off again by the next redistribution.

    :param charges:
        The charges of the atoms that stay.
    :param cap_indices:
        Which of those are capping hydrogens.

    :return:
        The tuple of the compensated charges and the shift each atom was
        given.
    """

    charges = np.array(charges, dtype=float)
    caps = set(cap_indices)
    rest = [index for index in range(charges.size) if index not in caps]

    carriers = len(rest) > 0
    assert_msg_critical(
        carriers, 'MetalForceFieldManager: discarding the residues that '
        'coordinate no metal leaves nothing but capping hydrogens to '
        'carry the charge')

    total = float(np.sum(charges))
    shift = (round(total) - total) / len(rest)
    charges[rest] += shift

    return charges, shift


def _print_discarded_residues(template, nodes, discard, shift, ostream=None):
    """
    Reports what loading a template threw away.

    :param template:
        The template as it stood before the atoms were dropped.
    :param nodes:
        The coarse nodes of the discarded residues.
    :param discard:
        The indices of the atoms that were dropped.
    :param shift:
        The charge every remaining atom was given to make up for them.
    """

    ostream = stream(ostream)

    coarse = template['coarse_topology']
    charges = template['charges']

    named = ', '.join(
        sorted(f'{coarse.nodes[node]["formula"]}/'
               f'{coarse.nodes[node]["key"][:6]}' for node in nodes))
    carried = float(sum(charges[index] for index in discard))

    ostream.print_info(
        f'Template {template["name"]}: discarding {len(nodes)} residue(s) '
        f'coordinating no metal ({named}), {len(discard)} atom(s) '
        f'carrying {carried:+.3f} e.')
    ostream.print_info(
        f'Shifting every remaining atom but the capping hydrogens by '
        f'{shift:+.4f} e to make the charge of the template whole again.')

    stray = carried - round(carried)
    if abs(stray) > 0.25:
        ostream.print_warning(
            f'Those residues carry {carried:+.3f} e, which is '
            f'{stray:+.3f} e away from a whole number of electrons, so '
            'what is left of the template was rounded onto the nearest '
            'integer and that may not be the charge of the site it '
            'describes.')

    ostream.flush()

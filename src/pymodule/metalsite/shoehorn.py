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
Walking an active site onto a template by editing it.

compare_active_site refuses a template unless the site is made of the same
residues, coordinated the same way and protonated the same way, and each of
those three is something the cutoffs can get wrong on an unrelaxed
structure rather than something about the chemistry. This walks the site
onto the template's terms instead of refusing it.

Every edit goes through the builder's own public edit methods, so what
comes out is a site the builder could have been walked to by hand, and a
failure part way through is undone by putting the request back and
rebuilding from it.

Nothing here holds state. Which template a site was walked onto is the
manager's to remember, so run() returns whether it worked and the manager
records it.
"""

from collections import Counter
from copy import deepcopy
from itertools import permutations
import numpy as np
import math
import sys

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    pass

from ..outputstream import OutputStream
from ..errorhandler import assert_msg_critical
from . import core
from . import matching
from . import printing
from .printing import stream


def run(builder,
        template,
        max_include_radius,
        max_mappings=matching.DEFAULT_MAX_MAPPINGS,
        ostream=None):
    """
    Edits a builder's active site until it is built the way a template is.

    :param builder:
        The builder holding the site. It is edited in place, and left as it
        was found if anything stands in the way.
    :param template:
        The template to walk it onto.
    :param max_include_radius:
        How far out from a metal center, in Angstrom, a residue may be
        picked up from.
    :param max_mappings:
        The limit on how many atom mappings are built.
    :param ostream:
        The output stream, or None to report nothing.

    :return:
        True when the site was walked onto the template, False when
        something the structure does not hold stood in the way.
    """

    ostream = stream(ostream)

    # everything the builder is is a function of the record of the edits
    # made on it, which is what its own edit methods write and what
    # build_active_site() rebuilds from. Keeping a copy of it is therefore
    # the whole of undoing a run that fails part way through.
    snapshot = deepcopy(builder._request)

    printing.print_shoehorn_header(
        template['name'],
        len(matching.residue_nodes(template['coarse_topology'])),
        builder.active_site['residues'],
        max_include_radius,
        ostream=ostream)

    # Every edit rebuilds the site, and each rebuild reports the whole
    # cluster and relaxes it again. Neither says anything about what to edit
    # next -- the decisions are read off the protein, which the crude
    # relaxation does not touch -- so the builder is quietened and its
    # relaxation held back until there is a site worth relaxing. The stream
    # is swapped rather than muted because mute is reference counted, and
    # the caller's own stream may be this very object.
    builder_stream = builder.ostream
    relaxing = builder._mm_opt

    builder.ostream = OutputStream(None)
    builder._mm_opt = False

    try:
        reason = _shoehorn(builder,
                           template,
                           float(max_include_radius),
                           max_mappings=max_mappings,
                           ostream=ostream)
    except Exception:
        builder._mm_opt = relaxing
        _restore(builder, snapshot)
        builder.ostream = builder_stream
        raise

    builder._mm_opt = relaxing

    if reason is not None:
        _restore(builder, snapshot)
        builder.ostream = builder_stream
        ostream.print_warning(
            f'Could not shoehorn the site into {template["name"]}: {reason}. '
            'The active site is as it was.')
        ostream.flush()
        return False

    if relaxing:
        # what comes out is what a build_active_site would have left
        builder.build_active_site()

    builder.ostream = builder_stream

    _print_shoehorn_summary(builder, template['name'], ostream=ostream)

    return True


def described_site(builder):
    """
    Describes the active site of a builder as it now stands.

    Derived on every call rather than kept: every edit rebuilds the site
    and renumbers its atoms, so a description from before one says
    nothing about the site after it.

    :param builder:
        The builder holding the site.

    :return:
        The active site with its two topologies added.
    """

    active_site = builder.active_site

    return matching.describe(
        active_site,
        core.connectivity_bonds(active_site['connectivity_matrix']))


def _shoehorn(builder,
              template,
              max_include_radius,
              max_mappings=matching.DEFAULT_MAX_MAPPINGS,
              ostream=None):
    """
    The three stages of a shoehorning, in order.

    :param builder:
        The builder holding the site to edit.
    :param template:
        The template to edit it onto.
    :param max_include_radius:
        How far out from a metal center a residue may be picked up from.

    :return:
        What stood in the way, or None when nothing did.
    """

    reason = _shoehorn_composition(builder,
                                   template,
                                   max_include_radius,
                                   ostream=ostream)
    if reason is not None:
        return reason

    # the protonation first: the mapping the denticity is forced through
    # cannot be solved while the hydrogens still differ
    for stage in (_shoehorn_protonation, _shoehorn_denticity):
        reason = stage(builder,
                       template,
                       max_mappings=max_mappings,
                       ostream=ostream)
        if reason is not None:
            return reason

    return None


def _restore(builder, request):
    """
    Puts an active site back the way it was before a shoehorning.

    :param builder:
        The builder holding the site.
    :param request:
        The record of its edits, as it stood beforehand.
    """

    builder._request = deepcopy(request)
    builder.build_active_site()


def _print_shoehorn_summary(builder, name, ostream=None):
    """
    Names what the site was walked onto and prints what it became.

    :param builder:
        The builder holding the site.
    :param name:
        The name of the template it was edited onto.
    """

    # bound once: binding_modes derives a new object per access
    modes = builder.binding_modes
    residues = list(builder.protonated_topology.residues())
    site = set(core.active_site_residues(modes))

    variants = [(core.residue_label(residues[res_index]), variant)
                for res_index, variant in sorted(modes['variants'].items())
                if res_index in site]

    printing.print_shoehorn_summary(name, modes, variants, ostream=ostream)


def _candidate_residues(builder, max_radius):
    """
    Every residue a shoehorning is allowed to build the site out of.

    The site's own residues and everything else within max_radius of a
    metal center, each measured twice: how close its sidechain reaches a
    center at all, and how close its nearest donor atom comes, which is
    what a bond would be made through.

    :param builder:
        The builder holding the site.
    :param max_radius:
        How far out from a metal center, in Angstrom, to look.

    :return:
        The candidates, each naming its residue, its family and its
        distances to every metal center by that center's residue index.
    """

    topology = builder.protonated_topology
    positions = np.asarray(builder.enzyme_positions)
    modes = builder.binding_modes

    members = set(core.active_site_residues(modes))
    metals = {
        entry['res_index']: positions[entry['index']]
        for entry in modes['metals']
    }

    candidates = []

    for residue in topology.residues():
        if residue.name in core.UNTRUNCATABLE_RESIDUES:
            continue

        atoms = matching.sidechain_heavy_atoms(residue)
        if not atoms:
            continue

        reach = {}
        donor = {}

        for res_index, metal in metals.items():
            measured = [(float(np.linalg.norm(positions[atom.index] - metal)),
                         atom) for atom in atoms]
            reach[res_index] = min(distance for distance, _ in measured)
            donors = [(distance, atom) for distance, atom in measured
                      if atom.element.symbol in core.DONOR_ELEMENTS]
            if donors:
                distance, atom = min(donors, key=lambda found: found[0])
                donor[res_index] = (distance, atom.name)

        if residue.index not in members and min(reach.values()) > max_radius:
            continue

        candidates.append({
            'res_index': residue.index,
            'resid': str(residue.id),
            'chain': str(residue.chain.id),
            'name': residue.name,
            'label': core.residue_label(residue),
            'key': matching.residue_family_key(topology, residue),
            'reach': reach,
            'donor': donor,
            'member': residue.index in members,
        })

    return candidates


def _template_slots(template):
    """
    What a template asks the site to be made of: one slot per residue of
    it, naming the kind of residue and the metal centers it is on.

    :param template:
        The template.

    :return:
        The slots, in the order the coarse topology holds them.
    """

    coarse = template['coarse_topology']
    slots = []

    for node in matching.residue_nodes(coarse):
        slots.append({
            'key':
            coarse.nodes[node]['family'],
            'formula':
            coarse.nodes[node]['formula'],
            'metals':
            sorted(image[1] for image in coarse.neighbors(node)),
        })

    return slots


def _metal_pairings(builder, template, described):
    """
    Every way the metal centers of a site can stand for those of a
    template.

    There are only ever a handful, so they are all tried and the one
    whose residues fit best is taken, rather than being guessed at from
    what coordinates them. The element has to agree outright: a zinc site
    is not a template for an iron one.

    The site's centers are named by the index of their residue, which no
    rebuild disturbs, unlike the atom indices every rebuild renumbers.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.
    :param described:
        The described active site.

    :return:
        One mapping of template metal index to site metal residue index
        per pairing.
    """

    template_labels = template['molecule'].get_labels()
    query_labels = described['molecule'].get_labels()

    template_metals = list(template['metal_indices'])
    query_metals = list(described['metal_indices'])

    if len(template_metals) != len(query_metals):
        return []

    pairings = []

    # bound once for the whole search: nothing here edits the site, and
    # every permutation would otherwise re-derive the same coordination
    modes = builder.binding_modes

    for order in permutations(query_metals):
        if any(template_labels[first] != query_labels[second]
               for first, second in zip(template_metals, order)):
            continue
        pairings.append({
            first: _metal_res_index(modes, described, second)
            for first, second in zip(template_metals, order)
        })

    return pairings


def _slot_cost(slot, candidate, pairing, max_radius):
    """
    What it costs to build one of a template's residues out of one of the
    site's, as the total distance the bonds it asks for would span.

    :param slot:
        The template residue.
    :param candidate:
        The residue of the structure.
    :param pairing:
        Which site metal center stands for which of the template's.
    :param max_radius:
        How far from a metal center a residue may be picked up from,
        measured on the residue rather than on the atom it would bind
        through.

    :return:
        The cost, or infinity when this residue cannot fill this slot.
    """

    if candidate['key'] != slot['key']:
        return math.inf

    metals = [pairing[metal] for metal in slot['metals']]

    if len(metals) > 1 and candidate['name'] not in core.BRIDGING_RESIDUES:
        # an imidazole nitrogen has one lone pair in the ring plane, so a
        # second metal on it is an artifact however the template reads
        return math.inf

    if not metals:
        # a residue the template holds without coordinating anything;
        # what matters is only that it is in the same neighbourhood
        distance = min(candidate['reach'].values())
        return distance if distance <= max_radius else math.inf

    total = 0.0

    for metal in metals:
        # Whether the residue is in reach is asked of the residue, not of
        # the atom it would bind through. A sidechain turns about its own
        # bonds, so a histidine sitting against the metal with its ring
        # rotated away has its nitrogens further off than its carbons --
        # which says how it is posed, not whether it belongs to the site,
        # and posing it is what the shoehorning goes on to do. Gating on
        # the donor atom rejected exactly the residues worth reaching for.
        if candidate['reach'][metal] > max_radius:
            return math.inf

        found = candidate['donor'].get(metal)

        # nothing to bind through at all, however close it comes
        if found is None:
            return math.inf

        # the bond as it stands is still what ranks one admissible
        # candidate against another: a sidechain already pointing at the
        # metal is likelier the one the template means
        total += found[0]

    return total


def _assign_residues(slots, candidates, pairing, max_radius):
    """
    Solves which residue of the structure fills which residue of the
    template, for one pairing of the metal centers.

    The whole site is assigned at once rather than one metal center at a
    time, which is what lets a residue move from one center to the other:
    a histidine the cutoffs put on the near metal is exactly the one the
    template wants on the far one, and no amount of filling one center's
    deficit from the outside will produce that.

    :param slots:
        The template's residues.
    :param candidates:
        The residues of the structure.
    :param pairing:
        Which site metal center stands for which of the template's.
    :param max_radius:
        How far from a metal center a residue may be picked up from.

    :return:
        The tuple of the total cost and the slot to candidate mapping, or
        None when some residue of the template cannot be filled at all.
    """

    assert_msg_critical(
        'scipy' in sys.modules,
        'MetalForceFieldManager.shoehorn: scipy is required to solve '
        'which residue of the structure fills which of the template')

    if len(candidates) < len(slots):
        return None

    costs = np.array([[
        _slot_cost(slot, candidate, pairing, max_radius)
        for candidate in candidates
    ] for slot in slots])

    finite = np.isfinite(costs)

    if not finite.any(axis=1).all():
        # a slot nothing at all can fill; no assignment exists
        return None

    # linear_sum_assignment cannot be given infinities, so they are made
    # more expensive than any whole assignment of finite ones
    blocked = costs[finite].sum() + 1.0
    rows, columns = linear_sum_assignment(np.where(finite, costs, blocked))

    assignment = {}
    total = 0.0

    for row, column in zip(rows, columns):
        if not finite[row, column]:
            return None
        assignment[int(row)] = candidates[int(column)]
        total += float(costs[row, column])

    if len(assignment) != len(slots):
        return None

    return total, assignment


def _assignment_failure(template,
                        slots,
                        candidates,
                        pairings,
                        max_radius,
                        ostream=None):
    """
    Says which of a template's residues the structure cannot supply.

    Reported off the pairing that got furthest, since that is the one the
    site came closest to being built like.

    :param template:
        The template that could not be matched.
    :param slots:
        Its residues.
    :param candidates:
        The residues of the structure.
    :param pairings:
        The pairings of the metal centers that were tried.
    :param max_radius:
        How far from a metal center a residue may be picked up from.

    :return:
        What stood in the way.
    """

    best = None

    for pairing in pairings:
        missing = []
        for slot in slots:
            if any(
                    math.isfinite(
                        _slot_cost(slot, candidate, pairing, max_radius))
                    for candidate in candidates):
                continue
            missing.append(slot)
        if best is None or len(missing) < len(best[1]):
            best = (pairing, missing)

    pairing, missing = best

    if missing:
        formulas = ', '.join(sorted(slot['formula'] for slot in missing))
        return (
            f'Template {template["name"]} is made of residues the structure '
            f'does not have within {max_radius:.1f} A of the right '
            f'metal center ({formulas})')

    # every residue of the template can be filled by something, but not
    # by enough different somethings at once: two of them are competing
    # for the one residue the structure has in reach.
    #
    # What is available is counted over the residues a slot of that kind
    # could actually be built from, which is not the same as the residues
    # of that kind the structure holds: a histidine whose ring nitrogen is
    # out of reach of the metal the template wants it on is no use here,
    # however close the rest of it comes. Counting it left this case
    # falling through to the message below, which says a competition was
    # lost without naming the distance that lost it.
    counts = Counter(slot['key'] for slot in slots)
    available = Counter()

    for key in counts:
        wanted = [slot for slot in slots if slot['key'] == key]
        available[key] = sum(
            any(
                math.isfinite(_slot_cost(slot, candidate, pairing, max_radius))
                for slot in wanted) for candidate in candidates)

    short = [
        slot['formula'] for slot in slots
        if available[slot['key']] < counts[slot['key']]
    ]

    if short:
        return (f'{template["name"]} is made of more '
                f'{", ".join(sorted(set(short)))} residues than the '
                f'structure has within {max_radius:.1f} A of its metal '
                'centers')

    return (f'the residues of the structure cannot be shared out over '
            f'those of {template["name"]}')


def _shoehorn_composition(builder, template, max_radius, ostream=None):
    """
    Makes the site hold the template's residues, coordinating the metal
    centers the template coordinates.

    Stages one and two of a shoehorning at once, because they are one
    question: which residue of the structure is to be which residue of
    the template. Answered as a single assignment over the whole site,
    for every way the metal centers can be paired up, and the cheapest
    answer is the one built -- so a residue is moved from one center to
    another where that is what the template says, rather than a center's
    deficit being filled from the outside while the residue that belongs
    in it sits on its neighbour.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.
    :param max_radius:
        How far out from a metal center a residue may be picked up from.

    :return:
        What stood in the way, or None when nothing did.
    """

    described = described_site(builder)

    # on the amino acids alone: the protonation is put right afterwards,
    # and until it is, an ASH is not the ASP a template holds
    if matching.coarse_mappings(template, described, match_protonation=False):
        return None

    slots = _template_slots(template)
    candidates = _candidate_residues(builder, max_radius)
    pairings = _metal_pairings(builder, template, described)

    if not pairings:
        return (f'the metal centers of {template["name"]} are not the '
                'metal centers of the site')

    best = None

    for pairing in pairings:
        found = _assign_residues(slots, candidates, pairing, max_radius)
        if found is None:
            continue
        if best is None or found[0] < best[0]:
            best = (found[0], pairing, found[1])

    if best is None:
        return _assignment_failure(template,
                                   slots,
                                   candidates,
                                   pairings,
                                   max_radius,
                                   ostream=ostream)

    _, pairing, assignment = best

    _apply_assignment(builder,
                      template,
                      slots,
                      pairing,
                      assignment,
                      ostream=ostream)

    if not matching.coarse_mappings(
            template, described_site(builder), match_protonation=False):
        return ('which residue coordinates which metal still differs from '
                f'{template["name"]}')

    return None


def _apply_assignment(builder,
                      template,
                      slots,
                      pairing,
                      assignment,
                      ostream=None):
    """
    Edits the site into the assignment that was solved.

    The bonds the template does not make go first and the ones it makes
    second, so that a residue moving from one metal center to another is
    never bonded to both at once -- which for a histidine add_metal_bond
    refuses outright. What is left holding nothing is dropped last, once
    every center has what it is owed.

    :param builder:
        The builder holding the site.
    :param template:
        The template being matched.
    :param slots:
        Its residues.
    :param pairing:
        Which site metal center stands for which of the template's.
    :param assignment:
        Which residue of the structure fills which slot.
    """

    wanted = {}
    donors = {}

    for index, candidate in assignment.items():
        res_index = candidate['res_index']
        wanted.setdefault(res_index,
                          set()).update(pairing[metal]
                                        for metal in slots[index]['metals'])
        donors[res_index] = candidate

    _include_assigned(builder, template, donors, ostream=ostream)
    _remove_unwanted_bonds(builder, wanted, ostream=ostream)
    _add_wanted_bonds(builder, wanted, donors, ostream=ostream)
    _drop_unassigned(builder, template, donors, ostream=ostream)


def _include_assigned(builder, template, donors, ostream=None):
    """
    Puts every residue the assignment uses into the site before any bond
    is made.

    A residue remove_residue took out is remembered as excluded, and a
    bond to a residue the site excludes is refused by the extraction
    rather than quietly putting it back, so asking for it by name is
    what clears the way. It is also how a residue the template holds
    without coordinating anything gets in at all.

    :param builder:
        The builder holding the site.
    :param template:
        The template being matched.
    :param donors:
        The candidate of each assigned residue.
    """

    members = set(core.active_site_residues(builder.binding_modes))
    included = []

    for res_index, candidate in sorted(donors.items()):
        if res_index in members:
            continue
        builder.include_residue(candidate['resid'], chain=candidate['chain'])
        included.append(candidate['label'])

    if included:
        ostream.print_info('Included ' + ', '.join(included) +
                           f', which {template["name"]} is made of.')
        ostream.flush()


def _remove_unwanted_bonds(builder, wanted, ostream=None):
    """
    Takes out every metal bond the assignment does not ask for.

    :param builder:
        The builder holding the site.
    :param wanted:
        The metal centers each residue is to coordinate, by residue
        index.
    """

    unwanted = []

    modes = builder.binding_modes
    metal_res = {
        entry['index']: entry['res_index']
        for entry in modes['metals']
    }

    for ligand in modes['ligands']:
        for metal in ligand['metals']:
            if metal_res[metal] in wanted.get(ligand['res_index'], set()):
                continue
            unwanted.append({
                'resid': ligand['residue'][3:],
                'chain': str(ligand['chain']),
                'atom': ligand['atom'],
                'metal': metal_res[metal],
                'label': f'{ligand["residue"]} {ligand["atom"]}',
            })

    for entry in unwanted:
        metal = _metal_entry(builder, entry['metal'])
        builder.remove_metal_bond(entry['resid'],
                                  metal=metal['index'],
                                  atom=entry['atom'],
                                  chain=entry['chain'])

    if unwanted:
        ostream.print_info('Unbound ' + ', '.join(entry['label']
                                                  for entry in unwanted) +
                           ', which the template does not coordinate.')
        ostream.flush()


def _add_wanted_bonds(builder, wanted, donors, ostream=None):
    """
    Makes every metal bond the assignment asks for and the site does not
    already have.

    :param builder:
        The builder holding the site.
    :param wanted:
        The metal centers each residue is to coordinate, by residue
        index.
    :param donors:
        The candidate of each of those residues, which names the atom the
        bond is made through.
    """

    added = []

    for res_index in sorted(wanted):
        for metal_res_index in sorted(wanted[res_index]):
            if _bonded(builder, res_index, metal_res_index):
                continue

            candidate = donors[res_index]
            found = candidate['donor'].get(metal_res_index)

            metal = _metal_entry(builder, metal_res_index)
            builder.add_metal_bond(candidate['resid'],
                                   metal['index'],
                                   atom=None if found is None else found[1],
                                   chain=candidate['chain'])
            added.append(f'{candidate["label"]} to '
                         f'{metal["element"]} {metal["index"]}')

    if added:
        ostream.print_info('Bonded ' + ', '.join(added) +
                           ', which the template coordinates.')
        ostream.flush()


def _bonded(builder, res_index, metal_res_index):
    """
    Whether one residue coordinates one metal center as things stand.

    :param builder:
        The builder holding the site.
    :param res_index:
        The residue.
    :param metal_res_index:
        The residue of the metal center.

    :return:
        True when the two are bonded.
    """

    modes = builder.binding_modes
    metal = [
        entry['index'] for entry in modes['metals']
        if entry['res_index'] == metal_res_index
    ]

    if not metal:
        return False

    return any(ligand['res_index'] == res_index and metal[0] in ligand['metals']
               for ligand in modes['ligands'])


def _drop_unassigned(builder, template, donors, ostream=None):
    """
    Drops the residues the assignment had no use for.

    Last of the edits, so that a residue on its way out is never the one
    a metal center is still waiting for, and remove_residue's refusal to
    orphan a center cannot be tripped by the order things are done in.

    :param builder:
        The builder holding the site.
    :param template:
        The template being matched.
    :param donors:
        The candidate of each assigned residue.
    """

    residues = list(builder.protonated_topology.residues())
    dropped = []

    for res_index in sorted(
            set(core.active_site_residues(builder.binding_modes)) -
            set(donors)):
        residue = residues[res_index]
        builder.remove_residue(str(residue.id), chain=str(residue.chain.id))
        dropped.append(core.residue_label(residue))

    if dropped:
        ostream.print_info('Dropped ' + ', '.join(dropped) +
                           f', which {template["name"]} is not made of.')
        ostream.flush()


def _metal_res_index(modes, described, metal):
    """
    The residue index of one of the site's metal centers.

    Takes the coordination rather than the builder, so a caller looking up
    several metals binds it once. binding_modes hands back a new object per
    access, so reading it per lookup rebuilds it per lookup -- which is only
    safe to avoid where no edit happens in between, and every caller here
    is inside a loop that makes none.

    :param modes:
        The binding modes, bound once by the caller.
    :param described:
        The described active site.
    :param metal:
        The index of the metal in that site.

    :return:
        The index of its residue in the topology.
    """

    index = described['atom_map'][metal]
    res_index = {
        entry['index']: entry['res_index']
        for entry in modes['metals']
    }

    known = index in res_index
    assert_msg_critical(
        known, 'MetalForceFieldManager: the active site holds a metal '
        f'center at atom {index} that its binding modes do not')

    return res_index[index]


def _metal_entry(builder, res_index):
    """
    The metal center of one residue, as the coordination now has it.

    Looked up again before every edit rather than kept: an edit
    reprotonates the structure, which renumbers its atoms, so the atom
    index an edit takes is only good until the next one.

    :param builder:
        The builder holding the site.
    :param res_index:
        The index of the metal's residue.

    :return:
        The metal entry of the binding modes.
    """

    entries = {
        entry['res_index']: entry
        for entry in builder.binding_modes['metals']
    }

    known = res_index in entries
    assert_msg_critical(
        known, 'MetalForceFieldManager: the structure no longer holds a '
        f'metal center in residue {res_index}')

    return entries[res_index]


def _best_heavy_mapping(template,
                        described,
                        match_h_count=True,
                        max_mappings=matching.DEFAULT_MAX_MAPPINGS,
                        ostream=None):
    """
    Solves which of the site's atoms is which of the template's.

    Every way the two line up at the coarse level is taken down to the
    atoms and superimposed, and the one that fits best is the answer --
    which is the same thing compare_active_site measures a template by,
    run here to decide what to edit rather than what to report.

    :param template:
        The template to match.
    :param described:
        The described active site.
    :param match_h_count:
        Whether an atom has to carry as many hydrogens as the one it maps
        onto, at both levels. Off for the pass that runs before the
        protonation is put right, where no residue that differs in it
        would map at all.

    :return:
        The mapping over the heavy atoms, or None when there is none.
    """

    maps = []

    for coarse_mapping in matching.coarse_mappings(
            template, described, match_protonation=match_h_count):
        maps.extend(
            matching.heavy_atom_maps(template,
                                     described,
                                     coarse_mapping,
                                     match_h_count=match_h_count,
                                     max_mappings=max_mappings,
                                     ostream=ostream))

    if not maps:
        return None

    coordinates = described['molecule'].get_coordinates_in_angstrom()
    heavy_map, _, _ = matching.best_heavy_map(template, maps, coordinates)

    return heavy_map


def _hydrogen_count(described, index):
    """
    How many hydrogens one atom of a site carries.

    The capping hydrogen is one of them, on both sides of a comparison
    alike, so it cancels rather than having to be told apart here.

    :param described:
        A described active site or a template.
    :param index:
        The atom.

    :return:
        The count.
    """

    labels = described['molecule'].get_labels()

    return sum(1 for other in described['fine_topology'].neighbors(index)
               if labels[other] == 'H')


def _shoehorn_protonation(builder,
                          template,
                          max_mappings=matching.DEFAULT_MAX_MAPPINGS,
                          ostream=None):
    """
    Protonates every residue of the site the way the template has it.

    Run before the denticity: an atom mapping that has to agree on the
    hydrogens cannot be solved while they differ, so this one is solved
    without them and is the last thing that has to be.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.

    :return:
        What stood in the way, or None when nothing did.
    """

    described = described_site(builder)
    heavy_map = _best_heavy_mapping(template,
                                    described,
                                    match_h_count=False,
                                    max_mappings=max_mappings,
                                    ostream=ostream)

    if heavy_map is None:
        return ('the atoms of the site do not map onto those of '
                f'{template["name"]}')

    changes = _protonation_changes(builder,
                                   template,
                                   described,
                                   heavy_map,
                                   max_mappings=max_mappings,
                                   ostream=ostream)

    if isinstance(changes, str):
        return changes

    for change in changes:
        builder.update_protonation_state(change['resid'],
                                         change['variant'],
                                         chain=change['chain'])

    if changes:
        ostream.print_info('Set ' +
                           ', '.join(f'{change["label"]} to {change["variant"]}'
                                     for change in changes) +
                           f', which is how {template["name"]} is protonated.')
        ostream.flush()

    return None


def _protonation_changes(builder,
                         template,
                         described,
                         heavy_map,
                         max_mappings=matching.DEFAULT_MAX_MAPPINGS,
                         ostream=None):
    """
    Works out which residues are protonated unlike the template.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.
    :param described:
        The described active site.
    :param heavy_map:
        The mapping from template index to site index.

    :return:
        What to hand to update_protonation_state, or what stood in the
        way as a string.
    """

    topology = builder.protonated_topology
    atoms = list(topology.atoms())
    residues = list(topology.residues())
    variants = builder.binding_modes['variants']
    metals = set(described['metal_indices'])

    per_residue = {}

    for first, second in heavy_map.items():
        if second in metals:
            continue
        atom = atoms[described['atom_map'][second]]
        per_residue.setdefault(atom.residue.index, []).append(
            (first, atom.name))

    changes = []

    for res_index, pairs in sorted(per_residue.items()):
        residue = residues[res_index]
        current = variants.get(res_index, residue.name)

        wanted = {
            name: _hydrogen_count(template, first)
            for first, name in pairs
        }
        delta = sum(wanted.values()) - sum(
            _hydrogen_count(described, heavy_map[first]) for first, _ in pairs)

        if delta == 0 and not _tautomer_differs(
                template, described, heavy_map, pairs,
                max_mappings=max_mappings):
            continue

        variant = _target_variant(residue, current, delta, wanted)

        if variant is None:
            return (f'{core.residue_label(residue)} would have to be '
                    'protonated the way '
                    f'{template["name"]} has it, and there is no variant '
                    'of it that is')

        if variant == current:
            continue

        changes.append({
            'resid': str(residue.id),
            'chain': str(residue.chain.id),
            'label': core.residue_label(residue),
            'variant': variant,
        })

    return changes


def _tautomer_differs(template,
                      described,
                      heavy_map,
                      pairs,
                      max_mappings=matching.DEFAULT_MAX_MAPPINGS):
    """
    Whether a residue carries its hydrogens on other atoms than the
    template does, with the same number of them.

    The histidines are what this is about: HID and HIE hold one hydrogen
    each and differ only in which ring nitrogen holds it.

    :param template:
        The template to match.
    :param described:
        The described active site.
    :param heavy_map:
        The mapping from template index to site index.
    :param pairs:
        The atoms of the residue, as template index and site atom name.

    :return:
        True when some atom of it disagrees.
    """

    return any(
        _hydrogen_count(template, first) != _hydrogen_count(
            described, heavy_map[first]) for first, _ in pairs)


def _target_variant(residue, current, delta, wanted):
    """
    The protonation variant that gives a residue the template's
    hydrogens.

    A proton is a charge, so what the variant has to be is looked up by
    what the charge has to become: VARIANT_CHARGES is the table the
    active site charge is counted with and known_variants is what OpenMM
    will build. The histidines are the one family where two variants
    share a charge, and there the ring nitrogen that carries the hydrogen
    is what tells HID from HIE -- the same convention the automatic
    choice in _histidine_variant is made on.

    :param residue:
        The residue.
    :param current:
        The variant it is protonated as now.
    :param delta:
        How many hydrogens the template has on it more than the site
        does.
    :param wanted:
        How many hydrogens the template carries on each of its atoms, by
        atom name.

    :return:
        The variant, or None when no variant of the residue is it.
    """

    legal = core.known_variants(residue.name)
    charge = core.VARIANT_CHARGES.get(current)

    if not legal:
        # OpenMM builds it one way only, so as long as it carries as many
        # hydrogens as the template it is already protonated the only way
        # it can be. The relaxed mapping is free to swap two symmetric
        # heavy atoms, which is what puts one here with delta zero.
        return current if delta == 0 else None

    if charge is None:
        return None

    found = [
        variant for variant in legal
        if core.VARIANT_CHARGES.get(variant) == charge + delta
    ]

    if len(found) == 1:
        return found[0]

    if not found:
        return None

    # HID carries its hydrogen on ND1 and HIE on NE2
    tautomers = {'ND1': 'HID', 'NE2': 'HIE'}
    named = {
        tautomers[name]
        for name, count in wanted.items() if count and name in tautomers
    }
    named &= set(found)

    return named.pop() if len(named) == 1 else None


def _shoehorn_denticity(builder,
                        template,
                        max_mappings=matching.DEFAULT_MAX_MAPPINGS,
                        ostream=None):
    """
    Bonds every metal center to exactly the atoms the template bonds it
    to.

    How many atoms of a residue reach a metal is a distance cutoff on an
    unrelaxed structure rather than chemistry, which is why the matching
    refuses to have an opinion about it. What matched on those terms is
    built on the template's terms, and this is where the site is put onto
    them -- through the builder's own edits, so that the site itself says
    what it is bonded like rather than a force field built from it.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.

    :return:
        What stood in the way, or None when nothing did.
    """

    described = described_site(builder)
    heavy_map = _best_heavy_mapping(template,
                                    described,
                                    max_mappings=max_mappings,
                                    ostream=ostream)

    if heavy_map is None:
        return ('the atoms of the site do not map onto those of '
                f'{template["name"]} even once it is protonated like it')

    changes = _denticity_changes(builder,
                                 template,
                                 described,
                                 heavy_map,
                                 max_mappings=max_mappings,
                                 ostream=ostream)

    for change in changes['added']:
        entry = _metal_entry(builder, change['metal'])
        builder.add_metal_bond(change['resid'],
                               entry['index'],
                               atom=change['atom'],
                               chain=change['chain'])

    for change in changes['removed']:
        entry = _metal_entry(builder, change['metal'])
        builder.remove_metal_bond(change['resid'],
                                  metal=entry['index'],
                                  atom=change['atom'],
                                  chain=change['chain'])

    if changes['added'] or changes['removed']:
        ostream.print_info(
            f'Wired the metal centers as {template["name"]} wires its '
            f'own: {len(changes["added"])} metal bond(s) added, '
            f'{len(changes["removed"])} removed.')
        ostream.flush()

    return None


def _denticity_changes(builder,
                       template,
                       described,
                       heavy_map,
                       max_mappings=matching.DEFAULT_MAX_MAPPINGS,
                       ostream=None):
    """
    Works out which metal-ligand bonds the site makes and the template
    does not, and the other way round.

    Everything is named by residue id, atom name and the residue of the
    metal, none of which a rebuild disturbs, since applying the first of
    these changes renumbers the atoms the rest were found on.

    :param builder:
        The builder holding the site.
    :param template:
        The template to match.
    :param described:
        The described active site.
    :param heavy_map:
        The mapping from template index to site index.

    :return:
        What to add and what to remove.
    """

    metals = set(described['metal_indices'])
    bonds, _ = matching.metal_keys(template)

    wanted = {
        frozenset((heavy_map[first], heavy_map[second]))
        for first, second in bonds
    }
    current = {
        frozenset((first, second))
        for first, second in core.connectivity_bonds(
            described['connectivity_matrix']) if metals & {first, second}
    }

    changes = {'added': [], 'removed': []}

    # Bound once for the recording loop. Safe only because this loop
    # applies no edit -- the edits are made by the caller, on what is
    # returned -- and both would otherwise be rebuilt per bond: the atom
    # list walks the whole protein, and binding_modes derives a new object
    # per access by design.
    atoms = list(builder.protonated_topology.atoms())
    modes = builder.binding_modes

    for kind, pairs in (('added', wanted - current), ('removed',
                                                      current - wanted)):
        for pair in sorted(pairs, key=sorted):
            changes[kind].append(_bond_record(atoms, modes, described, pair))

    return changes


def _bond_record(atoms, modes, described, pair):
    """
    Names one metal-ligand bond the way an edit method takes it.

    :param atoms:
        The atoms of the protonated topology, bound once by the caller.
    :param modes:
        The binding modes, bound once by the caller.
    :param described:
        The described active site.
    :param pair:
        The two atoms of the bond, as active site indices.

    :return:
        The residue, its chain, the donor atom and the metal's residue.
    """

    metals = set(described['metal_indices'])
    metal, donor = sorted(pair, key=lambda index: index not in metals)

    atom = atoms[described['atom_map'][donor]]

    return {
        'resid': str(atom.residue.id),
        'chain': str(atom.residue.chain.id),
        'atom': atom.name,
        'metal': _metal_res_index(modes, described, metal),
        'label': f'{core.residue_label(atom.residue)} {atom.name}',
    }
